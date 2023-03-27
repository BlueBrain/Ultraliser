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

#include "Bounds.h"

#include <algorithm>
#include <cassert>
#include <cstddef>

namespace Ultraliser
{
Bounds::Bounds(const Vector3f &min, const Vector3f &max)
    : _min(min)
    , _max(max)
{
    assert(_min < _max);
}

void Bounds::expand(const Vector3f &point)
{
    for (size_t i = 0; i < 3; ++i)
    {
        _min[i] = std::min(_min[i], point[i]);
        _max[i] = std::max(_max[i], point[i]);
    }
}

void Bounds::expand(const Bounds &other)
{
    for (size_t i = 0; i < 3; ++i)
    {
        _min[i] = std::min(_min[i], other._min[i]);
        _max[i] = std::max(_max[i], other._max[i]);
    }
}

bool Bounds::intersects(const Vector3f &point) const
{
    return _min.x() <= point.x() && point.x() <= _max.x() && _min.y() <= point.y() && point.y() <= _max.y()
        && _min.z() <= point.z() && point.z() <= _max.z();
}

bool Bounds::intersects(const Bounds &other) const
{
    // slab comparsion
    auto xStart = _min.x();
    auto xEnd = _max.x();
    auto otherXStart = other._min.x();
    auto otherXEnd = other._max.x();

    if (xStart > otherXEnd || xEnd < otherXStart)
    {
        return false;
    }

    auto yStart = _min.y();
    auto yEnd = _max.y();
    auto otherYStart = other._min.y();
    auto otherYEnd = other._max.y();

    if (yStart > otherYEnd || yEnd < otherYStart)
    {
        return false;
    }

    auto zStart = _min.z();
    auto zEnd = _max.z();
    auto otherZStart = other._min.z();
    auto otherZEnd = other._max.z();

    if (zStart > otherZEnd || zEnd < otherZStart)
    {
        return false;
    }

    return true;
}

const Vector3f &Bounds::getMin() const noexcept
{
    return _min;
}

const Vector3f &Bounds::getMax() const noexcept
{
    return _max;
}

Vector3f Bounds::getCenter() const
{
    return (_min + _max) * 0.5f;
}

Vector3f Bounds::getDimensions() const
{
    return _max - _min;
}
}