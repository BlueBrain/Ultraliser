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

#include "SparseOctree.h"

#include "common/SparseOctreeNodeSlot.h"

#include <queue>
#include <stack>

namespace Ultraliser
{
SparseOctree::SparseOctree(const Bounds &bounds, uint8_t maxDepth)
    : _bounds(bounds)
    , _maxDepth(maxDepth)
    , _root(SparseOctreeNodeSlot::root)
{
}

SparseOctreeNode &SparseOctree::getRoot()
{
    return _root;
}

const SparseOctreeNode &SparseOctree::getRoot() const
{
    return _root;
}

uint8_t SparseOctree::getMaxDepth() const
{
    return _maxDepth;
}

const Bounds &SparseOctree::getBounds() const
{
    return _bounds;
}

void SparseOctree::compact()
{
    std::stack<SparseOctreeNode *> nonLeafNodes;

    std::queue<SparseOctreeNode *> searchQueue;
    searchQueue.push(&_root);

    // gather nodes with children
    while (!searchQueue.empty())
    {
        SparseOctreeNode *node = searchQueue.front();
        searchQueue.pop();
        for (size_t i = 0; i < node->getNumChildren(); ++i)
        {
            searchQueue.push(&node->getChild(i));
        }

        nonLeafNodes.push(node);
    }

    // compact from bottom to top
    while (!nonLeafNodes.empty())
    {
        SparseOctreeNode *node = nonLeafNodes.top();
        nonLeafNodes.pop();
        node->compact();
    }
}
}
