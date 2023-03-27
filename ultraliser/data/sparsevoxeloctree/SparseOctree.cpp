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

namespace
{
/*
class PointVoxelizationAlgorithm
{
public:
    PointVoxelizationAlgorithm(
        Ultraliser::SparseOctreeNode &root,
        const Ultraliser::Bounds &rootBounds,
        uint8_t maxDepth)
        : _root(root)
        , _rootBounds(rootBounds)
        , _maxDepth(maxDepth)
    {
    }

    void tryVoxelize(const Ultraliser::Vector3f &p)
    {
        if (!_rootBounds.intersects(p))
        {
            return;
        }

        uint32_t gridSize = 1 << _maxDepth;

        auto oneOverDimensions = 1.f / _rootBounds.getDimensions();
        auto gridPos = PositionConverter::worldToGrid((p - _rootBounds.getMin()) * oneOverDimensions, gridSize);
        auto gridMin = Point3UI(0);
        auto gridMax = Point3UI(gridSize);

        auto node = &_root;

        for (uint32_t i = 0; i < _maxDepth; ++i)
        {
            auto gridDims = gridMax - gridMin;
            auto nodeCenter = gridMin + Point3UI(gridDims.x() >> 1, gridDims.y() >> 1, gridDims.z() >> 1);

            uint8_t xShift = gridPos.x() < nodeCenter.x() ? 0 : 1;
            uint8_t yShift = gridPos.y() < nodeCenter.y() ? 0 : 2;
            uint8_t zShift = gridPos.z() < nodeCenter.z() ? 0 : 4;

            size_t targetChild = xShift + yShift + zShift;

            auto tableIdx = (static_cast<size_t>(node->getChildrenMask()) << 3) + targetChild;
            auto childOffset = OffsetTable::values[tableIdx];

            // The point falls where there are no children
            if (childOffset == OffsetTable::noEntry && (i + 1) <= _maxDepth)
            {
                // uint8_t targetChildMask = (1 << (xShift + yShift + zShift));
                uint8_t targetChildMask = 1;
                targetChildMask <<= xShift;
                targetChildMask <<= yShift;
                targetChildMask <<= zShift;

                node->addChildNode(targetChildMask);
                tableIdx = (static_cast<size_t>(node->getChildrenMask()) << 3) + targetChild;
                childOffset = OffsetTable::values[tableIdx];
                assert(childOffset != OffsetTable::noEntry);
            }

            node = &(node->getChild(childOffset));

            gridMin[0] = xShift == 0 ? gridMin.x() : nodeCenter.x();
            gridMax[0] = xShift == 0 ? nodeCenter.x() : gridMax.x();
            gridMin[1] = yShift == 0 ? gridMin.y() : nodeCenter.y();
            gridMax[1] = yShift == 0 ? nodeCenter.y() : gridMax.y();
            gridMin[2] = zShift == 0 ? gridMin.z() : nodeCenter.z();
            gridMax[2] = zShift == 0 ? nodeCenter.z() : gridMax.z();
        }
    }

private:
    Ultraliser::SparseOctreeNode &_root;
    const Ultraliser::Bounds &_rootBounds;
    uint8_t _maxDepth;
};

struct TriangleVertices
{
    const Ultraliser::Vector3f &a;
    const Ultraliser::Vector3f &b;
    const Ultraliser::Vector3f &c;
};

class TriangleBoxIntersection
{
public:
    static bool test(const Ultraliser::Bounds &bounds, const TriangleVertices &triangle)
    {
        auto center = bounds.getCenter();
        auto halfSize = bounds.getDimensions() * 0.5f;

        auto &a = triangle.a;
        auto &b = triangle.b;
        auto &c = triangle.c;

        double centerArray[3] = {center.x(), center.y(), center.z()};
        double halfSizeArray[3] = {halfSize.x(), halfSize.y(), halfSize.z()};
        double vertices[3][3] = {{a.x(), a.y(), a.z()}, {b.x(), b.y(), b.z()}, {c.x(), c.y(), c.z()}};

        return Ultraliser::checkTriangleBoxIntersection(centerArray, halfSizeArray, vertices);
    }
};

class TriangleVoxelizationAlgorithm
{
private:
    struct TestNodeData
    {
        Ultraliser::SparseOctreeNode &parent;
        Ultraliser::Bounds bounds;
        uint8_t slot;
        uint8_t depth;
    };

public:
    TriangleVoxelizationAlgorithm(
        Ultraliser::SparseOctreeNode &root,
        const Ultraliser::Bounds &rootBounds,
        uint8_t maxDepth)
        : _root(root)
        , _rootBounds(rootBounds)
        , _maxDepth(maxDepth)
    {
    }

    void tryVoxelize(const TriangleVertices &triangle)
    {
        if (!TriangleBoxIntersection::test(_rootBounds, triangle))
        {
            return;
        }

        _tryVoxelize(_root, _rootBounds, triangle, 0);
        while (!_queue.empty())
        {
            auto &next = _queue.front();
            if (next.depth + 1 == _maxDepth)
            {
                _queue.pop();
                continue;
            }

            auto &node = *next.parent.findChild(next.slot);
            auto &bounds = next.bounds;
            _tryVoxelize(node, bounds, triangle, next.depth);
            _queue.pop();
        }
    }

private:
    void _tryVoxelize(
        Ultraliser::SparseOctreeNode &node,
        const Ultraliser::Bounds &nodeBounds,
        const TriangleVertices &triangle,
        uint8_t nodeDepth)
    {
        auto nextDepth = static_cast<uint8_t>(nodeDepth + 1);
        for (size_t i = 0; i < 8; ++i)
        {
            auto slot = static_cast<uint8_t>(1 << i);
            auto childBounds = SparseNodeBounds::fromParentBounds(nodeBounds, slot);
            if (!TriangleBoxIntersection::test(childBounds, triangle))
            {
                continue;
            }
            node.addChildNode(slot);
            _queue.push({node, childBounds, slot, nextDepth});
        }
    }

private:
    Ultraliser::SparseOctreeNode &_root;
    const Ultraliser::Bounds &_rootBounds;
    uint8_t _maxDepth;
    std::queue<TestNodeData> _queue;
};


class PointSampler
{
public:
    PointSampler(const Ultraliser::SparseOctreeNode &root, const Ultraliser::Bounds &rootBounds, uint8_t maxDepth)
        : _root(root)
        , _rootBounds(rootBounds)
        , _maxDepth(maxDepth)
    {
    }

    bool sample(const Ultraliser::Vector3f &p)
    {
        if (!_rootBounds.intersects(p))
        {
            return false;
        }

        uint32_t gridSize = 1 << _maxDepth;

        auto oneOverDimensions = 1.f / _rootBounds.getDimensions();
        auto gridPos = PositionConverter::worldToGrid((p - _rootBounds.getMin()) * oneOverDimensions, gridSize);
        auto gridMin = Point3UI(0);
        auto gridMax = Point3UI(gridSize);

        auto node = &_root;

        for (uint32_t i = 0; i < _maxDepth; ++i)
        {
            if (node->isLeaf())
            {
                return true;
            }

            auto gridDims = gridMax - gridMin;
            auto nodeCenter = gridMin + Point3UI(gridDims.x() >> 1, gridDims.y() >> 1, gridDims.z() >> 1);

            uint8_t xShift = gridPos.x() < nodeCenter.x() ? 0 : 1;
            uint8_t yShift = gridPos.y() < nodeCenter.y() ? 0 : 2;
            uint8_t zShift = gridPos.z() < nodeCenter.z() ? 0 : 4;

            size_t targetChild = xShift + yShift + zShift;

            auto tableIdx = (static_cast<size_t>(node->getChildrenMask()) << 3) + targetChild;
            auto childOffset = OffsetTable::values[tableIdx];

            // The point falls where there are no children
            if (childOffset == OffsetTable::noEntry)
            {
                return false;
            }

            node = &(node->getChild(childOffset));

            gridMin[0] = xShift == 0 ? gridMin.x() : nodeCenter.x();
            gridMax[0] = xShift == 0 ? nodeCenter.x() : gridMax.x();
            gridMin[1] = yShift == 0 ? gridMin.y() : nodeCenter.y();
            gridMax[1] = yShift == 0 ? nodeCenter.y() : gridMax.y();
            gridMin[2] = zShift == 0 ? gridMin.z() : nodeCenter.z();
            gridMax[2] = zShift == 0 ? nodeCenter.z() : gridMax.z();
        }

        return true;
    }

private:
    const Ultraliser::SparseOctreeNode &_root;
    const Ultraliser::Bounds &_rootBounds;
    uint8_t _maxDepth;
};
*/
}

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

        if (node->isLeaf())
        {
            continue;
        }

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
