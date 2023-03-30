/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Nadir Roman Guerrero < nadir.romanguerrero@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include "MorphologyVoxelizer.h"

#include "PointVoxelizer.h"
#include "TriangleVoxelizer.h"

#include <data/sparsevoxeloctree/common/PointToGridPosition.h>

namespace
{
class MorphologySampleVoxelizer
{
private:
    struct SampleBounds
    {
        Ultraliser::Vec3ui_32 min;
        Ultraliser::Vec3ui_32 max;
    };

public:
    MorphologySampleVoxelizer(Ultraliser::SparseOctree &octree)
        : _voxelizer(octree)
        , _oneOverBoundsDims(1.f / octree.getBounds().getDimensions())
        , _minBound(octree.getBounds().getMin())
        , _resolution(1 << octree.getMaxDepth())
        , _voxelSize(octree.getBounds().getDimensions() / _resolution)
    {
    }

    void voxelize(const Ultraliser::Sample &sample)
    {
        auto bb = _getSampleBounds(sample);

        for (uint32_t ix = bb.min.x(); ix <= bb.max.x(); ++ix)
        {
            for (uint32_t iy = bb.min.y(); iy <= bb.max.y(); ++iy)
            {
                for (uint32_t iz = bb.min.z(); iz <= bb.max.z(); ++iz)
                {
                    auto point = Ultraliser::Vec3ui_32(ix, iy, iz);
                    if (!_testSampleCubeIntersection(sample, point))
                    {
                        continue;
                    }
                    _voxelizer.voxelize(point);
                }
            }
        }
    }

    SampleBounds _getSampleBounds(const Ultraliser::Sample &sample)
    {
        auto radius = Ultraliser::Vector3f(sample.getRadius());
        auto pMin = sample.getPosition() - radius;
        auto pMax = sample.getPosition() + radius;
        auto vMin = Ultraliser::PointToGridPosition::convert(pMin, _minBound, _oneOverBoundsDims, _resolution);
        auto vMax = Ultraliser::PointToGridPosition::convert(pMax, _minBound, _oneOverBoundsDims, _resolution);
        return {vMin, vMax};
    }

    bool _testSampleCubeIntersection(const Ultraliser::Sample &sample, const Ultraliser::Vec3ui_32 &voxel)
    {
        auto voxelHalfSize = _voxelSize * 0.5f;

        auto voxelCenter = Ultraliser::Vector3f(
            _minBound[0] + (voxel[0] * _voxelSize[0]) + voxelHalfSize[0],
            _minBound[1] + (voxel[1] * _voxelSize[1]) + voxelHalfSize[1],
            _minBound[2] + (voxel[2] * _voxelSize[2]) + voxelHalfSize[2]);

        auto position = sample.getPosition() - voxelCenter;
        position[0] = abs(position[0]);
        position[1] = abs(position[1]);
        position[2] = abs(position[2]);

        auto r2 = sample.getRadius() * sample.getRadius();
        auto minDist = 0.0f;

        for (uint8_t i = 0; i < 3; ++i)
        {
            if (position[i] > voxelHalfSize[i])
            {
                auto dist = position[i] - voxelHalfSize[i];
                minDist += (dist * dist);
            }
        }
        return minDist <= r2;
    }

private:
    Ultraliser::PointVoxelizer _voxelizer;
    Ultraliser::Vector3f _oneOverBoundsDims;
    Ultraliser::Vector3f _minBound;
    uint32_t _resolution;
    Ultraliser::Vector3f _voxelSize;
};
}

namespace Ultraliser
{
void NeuronMorphologyVoxelizer::voxelize(SparseOctree &octree, const NeuronMorphology &morphology)
{
    auto sampleVoxelizer = MorphologySampleVoxelizer(octree);

    auto somaMesh = Ultraliser::Mesh(&morphology);
    Ultraliser::MeshVoxelizer::voxelize(octree, somaMesh);

    auto sections = morphology.getSections();
    for (size_t i = 0; i < sections.size(); ++i)
    {
        auto section = sections[i];
        auto samples = section->getSamples();

        if (samples.empty())
        {
            continue;
        }

        if (samples.size() == 1)
        {
            sampleVoxelizer.voxelize(*samples[0]);
            continue;
        }

        // Rasterize a polyline representing the section samples
        Ultraliser::MeshVoxelizer::voxelize(octree, Ultraliser::Mesh(samples));

        // Rasterize the first and last samples as spheres to fill any gaps
        sampleVoxelizer.voxelize(*samples.front());
        sampleVoxelizer.voxelize(*samples.back());
    }
}
}