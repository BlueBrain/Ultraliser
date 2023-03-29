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

#include "PointSampler.h"

#include <data/sparsevoxeloctree/common/OffsetTable.h>
#include <data/sparsevoxeloctree/common/PointToGridPosition.h>
#include <data/sparsevoxeloctree/common/SparseOctreeGridTraverser.h>

namespace Ultraliser
{
PointSampler::PointSampler(const SparseOctree &octree)
    : _root(octree.getRoot())
    , _boundMin(octree.getBounds().getMin())
    , _oneOverBoundsDims(1.f / octree.getBounds().getDimensions())
    , _maxDepth(octree.getMaxDepth())
    , _gridSize(1 << _maxDepth)
{
}

uint32_t PointSampler::getResolution() const
{
    return _gridSize;
}

bool PointSampler::sample(const Vector3f &point) const
{
    auto gridPos = PointToGridPosition::convert(point, _boundMin, _oneOverBoundsDims, _gridSize);
    return sample(gridPos);
}

bool PointSampler::sample(const Vec3ui_32 &gridPosition) const
{
    if (gridPosition.x() > _gridSize || gridPosition.y() > _gridSize || gridPosition.z() > _gridSize)
    {
        return false;
    }

    auto traverser = SparseOctreeGridTraverser(gridPosition, _gridSize);

    auto node = &_root;

    while (true)
    {
        auto mask = traverser.next();

        auto targetChild = mask.x() + mask.y() + mask.z();
        auto slot = static_cast<uint8_t>(1 << targetChild);
        auto childrenMask = node->getChildrenMask();

        if (!(childrenMask & slot))
        {
            return false;
        }

        if (node->getNumChildren() == 0)
        {
            break;
        }

        auto tableIdx = (static_cast<size_t>(childrenMask) << 3) + targetChild;
        auto childOffset = OffsetTable::of(tableIdx);
        node = &(node->getChild(childOffset));
    }

    return true;
}
}