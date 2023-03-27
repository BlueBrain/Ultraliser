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

#include "TriangleVoxelizer.h"

#include <geometry/Intersection.h>

#include <data/sparsevoxeloctree/common/SparseOctreeNodeBounds.h>

#include <queue>

namespace
{
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
            auto childBounds = Ultraliser::SparseOctreeNodeBounds::fromParentBounds(nodeBounds, slot);
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
}

namespace Ultraliser
{
TriangleVoxelizer::TriangleVoxelizer(SparseOctree &octree)
    : _root(octree.getRoot())
    , _rootBounds(octree.getBounds())
    , _maxDepth(octree.getMaxDepth())
{
}

void TriangleVoxelizer::voxelize(const Vector3f &a, const Vector3f &b, const Vector3f &c)
{
    auto voxelizer = TriangleVoxelizationAlgorithm(_root, _rootBounds, _maxDepth);
    voxelizer.tryVoxelize({a, b, c});
}
}