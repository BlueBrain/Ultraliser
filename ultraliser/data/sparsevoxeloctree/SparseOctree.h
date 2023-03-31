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

#pragma once

#include <cstdint>

#include "Bounds.h"
#include "SparseOctreeNode.h"

namespace Ultraliser
{
class SparseOctree
{
public:
    /*
    class VolumeView
    {
    public:
        VolumeView(SparseOctree &octree);

        int64_t getWidth() const;
        int64_t getHeight() const;
        int64_t getDepth() const;
        uint64_t getValueUI64(int64_t x, int64_t y, int64_t z) const;
        void getVoxelBoundingBox(int64_t x, int64_t y, int64_t z, Vector3f &min, Vector3f &max) const;

    private:
        SparseOctree &_octree;
        OctreePointSampler _sampler;
        Vector3f _voxelSize;
    };
    */
public:
    /**
     * @brief  Construct an Sparse octree to voxelize the given axis-aligned space bounds until a given maximum depth.
     *
     * @param minBound Minimum XYZ of the axis aligned bounds to voxelize.
     * @param maxBound Maximum XYZ of the axis aligned bounds to voxelize.
     * @param maxDepth Max children depth to reach.
     */
    SparseOctree(const Bounds &bounds, uint8_t maxDepth);

    /**
     * @brief Returns the root node of the sparse octree
     *
     * @return const SparseOctreeNode&
     */
    SparseOctreeNode &getRoot();
    const SparseOctreeNode &getRoot() const;

    /**
     * @brief Returns the sparse octree max depth.
     *
     * @return uint8_t
     */
    uint8_t getMaxDepth() const;

    /**
     * @brief Returns the sparse octree spatial bounds
     *
     * @return const Bounds<float>&
     */
    const Bounds &getBounds() const;

protected:
    Bounds _bounds;
    uint8_t _maxDepth;
    SparseOctreeNode _root;
};
}