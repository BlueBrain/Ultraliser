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

#include "SparseOctree.h"
#include "../util/Util.h"

#include <cmath>
#include <iostream>
#include <queue>
#include <stack>
#include <stdexcept>

namespace sc
{
  SparseOctree::SparseOctree()
   : _volumeBounds(Bounds3DF())
   , _volumeSize(Point3DF(0.f))
   , _maxDepth(0u)
  {
  }

  SparseOctree::SparseOctree(const Point3DF& minBound,
                             const Point3DF& maxBound,
                             const uint8_t maxDepth)
   : _volumeBounds(minBound, maxBound)
   , _volumeSize(maxBound - minBound)
   , _maxDepth(maxDepth)
  {
  }

  void SparseOctree::init(const Point3DF& minBound,
                          const Point3DF& maxBound,
                          const uint8_t maxDepth)
  {
    _volumeBounds = Bounds3DF(minBound, maxBound);
    _volumeSize = maxBound - minBound;
    _maxDepth = maxDepth;

    _root = SparseOctreeNode();
  }

  SparseOctreeNode& SparseOctree::getRoot()
  {
    return _root;
  }

  const SparseOctreeNode& SparseOctree::getRoot() const
  {
    return _root;
  }

  uint8_t SparseOctree::getMaxDepth() const
  {
    return _maxDepth;
  }

  const Bounds3DF& SparseOctree::getBounds() const
  {
    return _volumeBounds;
  }

  const Point3DF& SparseOctree::get3DSize() const
  {
    return _volumeSize;
  }

  void SparseOctree::compact()
  {
    std::stack<SparseOctreeNode*> nonLeafNodes;

    std::queue<SparseOctreeNode*> searchQueue;
    searchQueue.push(&_root);

    // gather nodes with children
    while(!searchQueue.empty())
    {
      SparseOctreeNode* node = searchQueue.front();
      searchQueue.pop();

      // Only process non-leaf nodes
      if(node->getNumChildren() > 0)
      {
        for(auto& child : node->_children)
          searchQueue.push(&child);

        nonLeafNodes.push(node);
      }
    }

    // compact from bottom to top
    while(!nonLeafNodes.empty())
    {
      SparseOctreeNode* node = nonLeafNodes.top();
      nonLeafNodes.pop();
      node->compact();
    }
  }

  void SparseOctree::printStats() const
  {
    // Common constants
    const size_t nodeSize = sizeof(SparseOctreeNode);

    size_t totalSize = 0;
    size_t numNodes = 0;
    size_t numLeafs = 0;

    std::queue<const SparseOctreeNode*> searchQueue;
    searchQueue.push(&_root);
    while(!searchQueue.empty())
    {
      const SparseOctreeNode* node = searchQueue.front();
      searchQueue.pop();

      totalSize += nodeSize;
      numNodes++;
      numLeafs = node->getNumChildren() == 0? numLeafs + 1 : numLeafs;

      for(const auto& child : node->_children)
        searchQueue.push(&child);
    }

    std::cout << "Tree stats:\n"
              << "\tMin bound: " << _volumeBounds.min.x << ", " 
              << _volumeBounds.min.y << ", " << _volumeBounds.min.z << "\n"
              << "\tMax bound: " << _volumeBounds.max.x << ", " 
              << _volumeBounds.max.y << ", " << _volumeBounds.max.z << "\n"
              << "\tTotal number of nodes: \t" << numNodes << "\n"
              << "\tTotal number of leaf nodes: \t" << numLeafs << "\n"
              << "\tTotal memory size (bytes): \t" << totalSize << std::endl;
  }
}
