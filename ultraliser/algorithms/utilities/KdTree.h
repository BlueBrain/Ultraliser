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

#pragma once

#include <cstddef>
#include <vector>
#include <limits>

#include <math/Vector3f.h>

namespace Ultraliser
{
    constexpr auto none = std::numeric_limits<std::size_t>::max();

    struct Node
    {
        Vector3f point = {0, 0, 0};
        std::size_t left = none;
        std::size_t right = none;

        Node() = default;

        explicit Node(const Vector3f &value)
            : point(value)
        {
        }
    };

    struct KdTree
    {
        std::size_t root = none;
        std::vector<Node> nodes;
    };

    KdTree createTree(const std::vector<Vector3f> &points);

    const Vector3f &findNearestPoint(const KdTree &tree, const Vector3f &point);
}
