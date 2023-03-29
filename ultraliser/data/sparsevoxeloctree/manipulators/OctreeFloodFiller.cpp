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

#include "OctreeFloodFiller.h"

#include "PointSampler.h"
#include "PointVoxelizer.h"

namespace
{
class MoldFactory
{
public:
    static Ultraliser::SparseOctree createMold(const Ultraliser::SparseOctree &octree)
    {
        auto sampler = Ultraliser::PointSampler(octree);

        auto mold = Ultraliser::SparseOctree(octree.getBounds(), octree.getMaxDepth());
        auto voxelizer = Ultraliser::PointVoxelizer(mold);

        _fillZSlices(sampler, voxelizer);
        _fillYSlices(sampler, voxelizer);
        _fillXSlices(sampler, voxelizer);

        mold.compact();

        return mold;
    }

private:
    static void _fillZSlices(const Ultraliser::PointSampler &sampler, Ultraliser::PointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(2, 1, 0);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillYSlices(const Ultraliser::PointSampler &sampler, Ultraliser::PointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(1, 0, 2);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillXSlices(const Ultraliser::PointSampler &sampler, Ultraliser::PointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(0, 1, 2);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillSlices(
        const Ultraliser::PointSampler &sampler,
        Ultraliser::PointVoxelizer &mold,
        const Ultraliser::Vec3ui_32 &pattern)
    {
        auto resolution = sampler.getResolution();
        for (uint32_t slice = 0; slice < resolution; ++slice)
        {
            _fillSlice(sampler, mold, pattern, slice);
        }
    }

    static void _fillSlice(
        const Ultraliser::PointSampler &sampler,
        Ultraliser::PointVoxelizer &mold,
        const Ultraliser::Vec3ui_32 &pattern,
        uint32_t slice)
    {
        auto resolution = sampler.getResolution();
        for (uint32_t column = 0; column < resolution; ++column)
        {
            _fillRow(sampler, mold, pattern, slice, column);
        }
    }

    static void _fillRow(
        const Ultraliser::PointSampler &sampler,
        Ultraliser::PointVoxelizer &mold,
        const Ultraliser::Vec3ui_32 &pattern,
        uint32_t slice,
        uint32_t column)
    {
        auto resolution = sampler.getResolution();
        for (uint32_t i = 0; i < resolution; ++i)
        {
            auto sideA = Ultraliser::Vec3ui_32();
            sideA[pattern[0]] = slice;
            sideA[pattern[1]] = column;
            sideA[pattern[2]] = i;
            if (!sampler.sample(sideA))
            {
                mold.voxelize(sideA);
                continue;
            }

            auto sideB = sideA;
            sideB[pattern[2]] = resolution - i - 1;
            if (!sampler.sample(sideB))
            {
                mold.voxelize(sideB);
                continue;
            }

            break;
        }
    }
};

class MoldCastFactory
{
public:
    static void cast(const Ultraliser::SparseOctree &mold, Ultraliser::SparseOctree &hollowCast)
    {
        auto sampler = Ultraliser::PointSampler(mold);
        auto voxelizer = Ultraliser::PointVoxelizer(hollowCast);

        auto resolution = sampler.getResolution();
        for (uint32_t x = 0; x < resolution; ++x)
        {
            for (uint32_t y = 0; y < resolution; ++y)
            {
                for (uint32_t z = 0; z < resolution; ++z)
                {
                    auto point = Ultraliser::Vec3ui_32(x, y, z);
                    auto filled = sampler.sample(point);
                    if (filled)
                    {
                        continue;
                    }

                    voxelizer.voxelize(point);
                }
            }
        }
    }
};
}

namespace Ultraliser
{
SparseOctree OctreeFloodFiller::createMold(const SparseOctree &octree)
{
    return MoldFactory::createMold(octree);
}

void OctreeFloodFiller::fill(SparseOctree &octree)
{
    auto mold = MoldFactory::createMold(octree);
    MoldCastFactory::cast(mold, octree);
}
}