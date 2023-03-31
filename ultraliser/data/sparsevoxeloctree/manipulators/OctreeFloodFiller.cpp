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

#include "OctreeCompacter.h"
#include "OctreePointSampler.h"
#include "OctreePointVoxelizer.h"

namespace
{
class MoldFactory
{
public:
    static Ultraliser::SparseOctree createMold(const Ultraliser::SparseOctree &octree)
    {
        auto sampler = Ultraliser::OctreePointSampler(octree);

        auto mold = Ultraliser::SparseOctree(octree.getBounds(), octree.getMaxDepth());
        auto voxelizer = Ultraliser::OctreePointVoxelizer(mold);

        _fillZSlices(sampler, voxelizer);
        _fillYSlices(sampler, voxelizer);
        _fillXSlices(sampler, voxelizer);

        Ultraliser::OctreeCompacter::compact(mold);

        return mold;
    }

private:
    static void _fillZSlices(const Ultraliser::OctreePointSampler &sampler, Ultraliser::OctreePointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(2, 1, 0);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillYSlices(const Ultraliser::OctreePointSampler &sampler, Ultraliser::OctreePointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(1, 0, 2);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillXSlices(const Ultraliser::OctreePointSampler &sampler, Ultraliser::OctreePointVoxelizer &mold)
    {
        auto pattern = Ultraliser::Vec3ui_32(0, 2, 1);
        _fillSlices(sampler, mold, pattern);
    }

    static void _fillSlices(
        const Ultraliser::OctreePointSampler &sampler,
        Ultraliser::OctreePointVoxelizer &mold,
        const Ultraliser::Vec3ui_32 &pattern)
    {
        auto resolution = sampler.getResolution();
        for (uint32_t slice = 0; slice < resolution; ++slice)
        {
            _fillSlice(sampler, mold, pattern, slice);
        }
    }

    static void _fillSlice(
        const Ultraliser::OctreePointSampler &sampler,
        Ultraliser::OctreePointVoxelizer &mold,
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
        const Ultraliser::OctreePointSampler &sampler,
        Ultraliser::OctreePointVoxelizer &mold,
        const Ultraliser::Vec3ui_32 &pattern,
        uint32_t slice,
        uint32_t column)
    {
        auto resolution = sampler.getResolution();
        auto advanceA = true;
        auto advanceB = true;
        for (uint32_t i = 0; i < resolution; ++i)
        {
            auto sideA = Ultraliser::Vec3ui_32();
            sideA[pattern[0]] = slice;
            sideA[pattern[1]] = column;
            sideA[pattern[2]] = i;
            if (advanceA && sampler.sample(sideA))
            {
                advanceA = false;
            }

            auto sideB = sideA;
            sideB[pattern[2]] = resolution - i - 1;
            if (advanceB && sampler.sample(sideB))
            {
                advanceB = false;
            }

            if (advanceA)
            {
                mold.voxelize(sideA);
            }
            if (advanceB)
            {
                mold.voxelize(sideB);
            }
            if (!advanceA && !advanceB)
            {
                break;
            }
        }
    }
};

class MoldCastFactory
{
public:
    static void cast(const Ultraliser::SparseOctree &mold, Ultraliser::SparseOctree &hollowCast)
    {
        auto sampler = Ultraliser::OctreePointSampler(mold);
        auto voxelizer = Ultraliser::OctreePointVoxelizer(hollowCast);

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

class ShellFillingExtractor
{
public:
    static Ultraliser::SparseOctree extract(
        const Ultraliser::SparseOctree &mold,
        const Ultraliser::SparseOctree &hollowCast)
    {
        auto moldSampler = Ultraliser::OctreePointSampler(mold);
        auto hollowSampler = Ultraliser::OctreePointSampler(hollowCast);

        auto octree = Ultraliser::SparseOctree(hollowCast.getBounds(), hollowCast.getMaxDepth());
        auto voxelizer = Ultraliser::OctreePointVoxelizer(octree);

        auto resolution = moldSampler.getResolution();
        for (uint32_t x = 0; x < resolution; ++x)
        {
            for (uint32_t y = 0; y < resolution; ++y)
            {
                for (uint32_t z = 0; z < resolution; ++z)
                {
                    auto point = Ultraliser::Vec3ui_32(x, y, z);
                    auto filled = moldSampler.sample(point) || hollowSampler.sample(point);
                    if (filled)
                    {
                        continue;
                    }

                    voxelizer.voxelize(point);
                }
            }
        }

        return octree;
    }
};
}

namespace Ultraliser
{
SparseOctree OctreeFloodFiller::generateFilling(const SparseOctree &octree)
{
    auto mold = MoldFactory::createMold(octree);
    return ShellFillingExtractor::extract(mold, octree);
}

void OctreeFloodFiller::fill(SparseOctree &octree)
{
    auto mold = MoldFactory::createMold(octree);
    MoldCastFactory::cast(mold, octree);
}
}