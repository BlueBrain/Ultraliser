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

#ifndef SIMCRUSHER_SPARSEVOXELOCTREE_H
#define SIMCRUSHER_SPARSEVOXELOCTREE_H

#include <cstdint>

#include <ultraliser/math/Vector3f.h>

#include "../util/Point3D.h"
#include "../util/Bounds.h"
#include "SparseOctreeNode.h"

namespace sc
{
class SparseOctree
{  
  public:
    /*
      * Constructs an empty octree with zero-dimension
      * and zero depth
      */ 
    SparseOctree();
    
    /*
      * Constructs a sparse octree which will allow to
      * voxelize the space within the given bounds at 
      * a maximun given resolution
      */
    SparseOctree(const Point3DF& minBound,
                  const Point3DF& maxBound,
                  const uint8_t maxDepth);

    /*
      * Initializes the octree to the given bounds 
      * and max depth. Any previous content is
      * discarded
      */ 
    void init(const Point3DF& minBound,
              const Point3DF& maxBound,
              const uint8_t maxDepth);

    /*
      * Returns a const reference to the root node
      */
    SparseOctreeNode& getRoot();

    /*
      * Returns a const reference to the root node
      */
    const SparseOctreeNode& getRoot() const;

    /*
      * Returns the octree maximun allowed depth
      */
    uint8_t getMaxDepth() const;

    /*
      * Returns the octree's bounds
      */
    const Bounds3DF& getBounds() const;

    /*
      * Returns the octree's 3D size
      */
    const Point3DF& get3DSize() const;

    /*
      * Attempts to voxelize the given value in the 3D space
      * represented by the octree if the given points falls within it
      */
    //void voxelizePoint(const Point3DF& p, const uint8_t value);

    /*
      * Attempts to sample the value at the given point if it
      * falls within the bounds of the octree
      */ 
    //uint8_t sample(const Point3DF& p) const;

    /*
      * Triggers the compact of the octree from the leaf nodes in a 
      * bottom-to-top manner
      */ 
    void compact();

    /*
      * Prints octree stats on standard output
      * (memory size, number of nodes, number of leaf nodes)
      */
    void printStats() const;

  protected:
    Bounds3DF _volumeBounds;
    Point3DF _volumeSize;

    uint8_t _maxDepth;

    SparseOctreeNode _root;
};
}

#endif