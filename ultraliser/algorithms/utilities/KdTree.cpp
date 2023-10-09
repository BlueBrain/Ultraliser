/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Adrien Fleury <adrien.fleury@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include "KdTree.h"
#include <common/Common.h>

namespace
{
    using namespace Ultraliser;

    using Node = KdTree::Node;

    constexpr auto none = KdTree::none;

    struct SearchResult
    {
        size_t index = none;
        float distance = 0.0f;
    };

    size_t nextAxis(size_t axis)
    {
        return (axis + 1) % 3;
    }

    size_t buildTree(std::vector< Node > &nodes, size_t first, size_t last, size_t axis)
    {
        if (last <= first)
        {
            return none;
        }

        auto pivot = first + (last - first) / 2;
        auto begin = nodes.begin();

        std::nth_element(
            begin + first,
            begin + pivot,
            begin + last,
            [=](auto &left, auto &right)
            { return left.point[axis] < right.point[axis]; });

        axis = nextAxis(axis);

        nodes[pivot].left = buildTree(nodes, first, pivot, axis);
        nodes[pivot].right = buildTree(nodes, pivot + 1, last, axis);

        return pivot;
    }

    size_t buildTree(std::vector<Node> &nodes)
    {
        auto first = 0;
        auto last = nodes.size();
        auto axis = 0;
        return buildTree(nodes, first, last, axis);
    }

    SearchResult search(
        const std::vector<Node> &nodes,
        const Vector3f &point,
        const SearchResult &best,
        size_t index,
        size_t axis)
    {
        auto result = best;

        if (index == none)
        {
            return result;
        }

        auto &node = nodes[index];
        auto &position = node.point;

        auto delta = position - point;
        auto distance = delta.abs();

        if (best.index == none || distance < best.distance)
        {
            result.index = index;
            result.distance = distance;
        }

        if (distance == 0)
        {
            return result;
        }

        auto dx = delta[axis];
        auto left = dx > 0;

        axis = nextAxis(axis);

        auto next = left ? node.left : node.right;
        result = search(nodes, point, result, next, axis);

        if (result.distance <= dx * dx)
        {
            return result;
        }

        next = left ? node.right : node.left;
        return search(nodes, point, result, next, axis);
    }
}

namespace Ultraliser
{
    KdTree KdTree::from(const std::vector<Vector3f> &points)
    {
        auto nodes = std::vector<Node>(points.begin(), points.end());
        auto root = buildTree(nodes);
        return KdTree(root, std::move(nodes));
    }

    KdTree::KdTree(const size_t &root, std::vector<Node> nodes)
        : _root(root),
          _nodes(std::move(nodes))
    {
        /// EMPTY
    }

    KdTree::NearestPoint KdTree::findNearestPoint(const Vector3f &point) const
    {
        if (_nodes.empty() || _root == none)
        {
            LOG_ERROR("Invalid Tree");
        }

        auto result = search(_nodes, point, {}, _root, 0);

        if (result.index == none)
        {
            LOG_ERROR("Unexpected Search Failure");
        }

        return { _nodes[result.index].point, result.distance };
    }
}
