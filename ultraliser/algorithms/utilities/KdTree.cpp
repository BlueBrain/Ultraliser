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

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace
{
    using namespace Ultraliser;

    using Node = KdTree::Node;

    constexpr auto none = KdTree::none;

    struct SearchResult
    {
        std::size_t index = none;
        float distance = 0.0f;
    };

    struct SearchContext
    {
        std::size_t index;
        std::size_t axis = 0;
        SearchResult result;

        explicit SearchContext(std::size_t root)
            : index(root)
        {
        }
    };

    std::size_t nextAxis(std::size_t axis)
    {
        return (axis + 1) % 3;
    }

    std::size_t buildTree(std::vector<Node> &nodes, std::size_t first, std::size_t last, std::size_t axis)
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

    std::size_t buildTree(std::vector<Node> &nodes)
    {
        auto first = 0;
        auto last = nodes.size();
        auto axis = 0;
        return buildTree(nodes, first, last, axis);
    }

    float computeDistance(const Vector3f &a, const Vector3f &b)
    {
        auto distance = 0.0f;
        for (auto i = std::size_t(0); i < 3; ++i)
        {
            auto d = a[i] - b[i];
            distance += d * d;
        }
        return std::sqrt(distance);
    }

    void search(SearchContext &context, const std::vector<Node> &nodes, const Vector3f &point)
    {
        auto &index = context.index;
        if (index == none)
        {
            return;
        }

        auto &node = nodes[index];
        auto &position = node.point;
        auto distance = computeDistance(position, point);

        auto &result = context.result;
        if (result.index == none || distance < result.distance)
        {
            result.index = index;
            result.distance = distance;
        }

        if (distance == 0)
        {
            return;
        }

        auto &axis = context.axis;
        distance = position[axis] - point[axis];

        axis = nextAxis(axis);

        auto left = distance > 0;
        index = left ? node.left : node.right;
        search(context, nodes, point);

        if (distance * distance >= result.distance)
        {
            return;
        }

        index = left ? node.right : node.left;
        search(context, nodes, point);
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

    KdTree::KdTree(std::size_t root, std::vector<Node> nodes)
        : _root(root),
          _nodes(std::move(nodes))
    {
    }

    KdTree::NearestPoint KdTree::findNearestPoint(const Vector3f &point) const
    {
        if (_nodes.empty() || _root == none)
        {
            throw std::runtime_error("Invalid tree");
        }

        auto context = SearchContext(_root);
        search(context, _nodes, point);
        auto &result = context.result;

        if (result.index == none)
        {
            throw std::runtime_error("Unexpected search failure");
        }

        return {_nodes[result.index].point, result.distance};
    }
}
