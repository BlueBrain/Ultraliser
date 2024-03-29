/***************************************************************************************************
 * Copyright (c) 2015 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * All rights reserved. Do not distribute without permission.
 * Author(s):
 *      Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of Brayns <https://github.com/BlueBrain/Brayns>
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301 USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/
#include "HeaderUtils.h"

namespace libNRRD
{
Vector3ui HeaderUtils::get3DSize(const NRRDHeader &header)
{
    const auto &sizes = header.sizes;
    assert(sizes.size() >= 3);

    Vector3ui result;
    const auto start = sizes.size() - 3;
    for (size_t i = start; i < sizes.size(); ++i)
    {
        result[i - start] = sizes[i];
    }
    return result;
}

Vector3f HeaderUtils::get3DDimensions(const NRRDHeader &header)
{
    Vector3f result(1.f);

    const auto &spaceDirections = header.spaceDirections;
    if (!spaceDirections)
    {
        return result;
    }

    const auto &directions = *spaceDirections;
    assert(directions.size() == 3);

    for (size_t i = 0; i < directions.size(); ++i)
    {
        const auto &direction = directions[i];
        Vector3f vector;
        for (size_t j = 0; j < 3; ++j)
        {
            vector[j] = direction[j];
        }

        result[i] = glm::length(vector);
    }
    return result;
}

size_t HeaderUtils::getVoxelDimension(const NRRDHeader &header)
{
    if (header.dimensions <= 3)
    {
        return 1;
    }

    return header.sizes.front();
}
}
