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

#include <geometry/Intersection.h>
#include <math/Vector.h>

#include <array>
#include <cmath>
#include <iostream>
#include <queue>
#include <stack>
#include <stdexcept>

namespace
{

class OffsetTableComputer
{
public:
    static std::array<uint8_t, 2048> compute(const uint8_t noEntryValue)
    {
        auto result = std::array<uint8_t, 2048>();
        result.fill(noEntryValue);

        for (size_t i = 0; i < 256; ++i)
        {
            auto idx = i << 3;
            auto offset = 0ul;

            for (size_t j = 0; j < 8; ++j)
            {
                auto mask = (1ul << j);
                if ((i & mask))
                {
                    result[idx + j] = offset;
                    offset++;
                }
            }
        }

        return result;
    }
};

class OffsetTable
{
public:
    inline static const uint8_t noEntry = 8;
    inline static const std::array<uint8_t, 2048> values = OffsetTableComputer::compute(noEntry);
};

using Point3UI = Ultraliser::Vec<3, uint32_t>;

class PositionConverter
{
public:
    static Point3UI worldToGrid(const Ultraliser::Vector3f &normalizedPosition, uint32_t gridSize)
    {
        auto gridPosF = Ultraliser::Vector3f(UI2F(gridSize)) * normalizedPosition;
        return Point3UI(F2UI32(gridPosF.x()), F2UI32(gridPosF.y()), F2UI32(gridPosF.z()));
    }
};

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
        assert(_bounds.intersects(p));

        uint32_t gridSize = 1 << _maxDepth;

        auto oneOverDimensions = 1.f / _rootBounds.getDimensions();
        auto gridPos = PositionConverter::worldToGrid((p - _rootBounds.getMin()) * oneOverDimensions, gridSize);
        auto gridMin = Point3UI(0);
        auto gridMax = Point3UI(gridSize);
        auto gridDims = gridMax - gridMin;

        auto node = &_root;

        for (uint32_t i = 0; i < _maxDepth; ++i)
        {
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
                uint8_t targetChildMask = (1 << (xShift + yShift + zShift));
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

class SparseNodeBounds
{
public:
    static Ultraliser::Bounds fromParentBounds(const Ultraliser::Bounds &parentBounds, uint8_t slotMask)
    {
        auto left = slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::backTopLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontTopLeft;
        auto bottom = slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::backBottomRight
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomLeft
            || slotMask & Ultraliser::SparseOctreeNodeSlot::frontBottomRight;
        auto back = slotMask < Ultraliser::SparseOctreeNodeSlot::frontBottomLeft;

        auto xMult = left ? 0.f : 1.f;
        auto yMult = bottom ? 0.f : 1.f;
        auto zMult = back ? 0.f : 1.f;
        auto mult = Ultraliser::Vector3f(xMult, yMult, zMult);

        auto min = parentBounds.getMin();
        auto center = parentBounds.getCenter();
        auto halfSize = center - min;

        auto childMin = min + halfSize * mult;
        auto childMax = center + halfSize * mult;
        return Ultraliser::Bounds(childMin, childMax);
    }
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
        uint8_t index;
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

        _tryVoxelize(_root, _rootBounds, triangle);
        for (uint8_t depth = 1; depth < _maxDepth && !_queue.empty(); ++depth)
        {
            auto &next = _queue.front();
            auto &node = next.parent.getChild(next.index);
            auto &bounds = next.bounds;
            _tryVoxelize(node, bounds, triangle);
        }
    }

private:
    void _tryVoxelize(
        Ultraliser::SparseOctreeNode &node,
        const Ultraliser::Bounds &nodeBounds,
        const TriangleVertices &triangle)
    {
        for (size_t i = 0; i < 8; ++i)
        {
            auto slot = static_cast<uint8_t>(1 << i);
            auto childBounds = SparseNodeBounds::fromParentBounds(nodeBounds, slot);
            if (!TriangleBoxIntersection::test(childBounds, triangle))
            {
                continue;
            }
            auto index = static_cast<uint8_t>(node.getNumChildren());
            node.addChildNode(slot);
            _queue.push({node, childBounds, index});
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
        auto gridDims = gridMax - gridMin;

        auto node = &_root;

        for (uint32_t i = 0; i < _maxDepth; ++i)
        {
            if (node->isLeaf())
            {
                return true;
            }

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
}

namespace Ultraliser
{
SparseOctree::SparseOctree(const Bounds &bounds, uint8_t maxDepth)
    : _bounds(bounds)
    , _maxDepth(maxDepth)
    , _root(SparseOctreeNodeSlot::root)
{
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

void SparseOctree::voxelizePoint(const Vector3f &p)
{
    auto voxelizer = PointVoxelizationAlgorithm(_root, _bounds, _maxDepth);
    voxelizer.tryVoxelize(p);
}

void SparseOctree::voxelizeTriangle(const Vector3f &a, const Vector3f &b, const Vector3f &c)
{
    auto triangle = TriangleVertices{a, b, c};
    auto voxelizer = TriangleVoxelizationAlgorithm(_root, _bounds, _maxDepth);
    voxelizer.tryVoxelize(triangle);
}

bool SparseOctree::sample(const Vector3f &p) const
{
    auto sampler = PointSampler(_root, _bounds, _maxDepth);
    return sampler.sample(p);
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
