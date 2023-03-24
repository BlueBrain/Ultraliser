/* Copyright (c) 2020, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of SimCrusher
 * <LINK>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#pragma once

#include <cstdint>

#include <math/Vector3f.h>

#include "Bounds.h"
#include "SparseOctreeNode.h"

namespace Ultraliser
{
class SparseOctree
{
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

    /**
     * @brief Voxelizes a point into the sparse octree if its inside the octree bounds.
     *
     * @param point The point to voxelize
     */
    void voxelizePoint(const Vector3f &point);

    /**
     * @brief Voxelizes a triangle into the sparse octree, if it is fully inside the octree
     * bounds.
     *
     * @param a Triangle vertex 1
     * @param b Triangle vertex 2
     * @param c Triangle vertex 3
     */
    void voxelizeTriangle(const Vector3f &a, const Vector3f &b, const Vector3f &c);

    /**
     * @brief Samples the octree at a given position, and returns true if it falls within a leaf node.
     *
     * @param p 3D point to sample
     * @return true if the point is within any octree leaf node
     */
    bool sample(const Vector3f &p) const;

    /**
     * @brief Compacts the octree by collapsing parents with 8 children into leaf nodes recursively.
     */
    void compact();

protected:
    Bounds _bounds;
    uint8_t _maxDepth;
    SparseOctreeNode _root;
};
}