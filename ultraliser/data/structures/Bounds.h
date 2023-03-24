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

#pragma once

#include <math/Vector3f.h>

#include <limits>

namespace Ultraliser
{
class Bounds
{
public:
    Bounds() = default;
    Bounds(const Vector3f &min, const Vector3f &max);

    void expand(const Vector3f &point);
    void expand(const Bounds &other);

    bool intersects(const Vector3f &point) const;
    bool intersects(const Bounds &other) const;

    const Vector3f &getMin() const noexcept;
    const Vector3f &getMax() const noexcept;
    Vector3f getCenter() const;
    Vector3f getDimensions() const;

private:
    Vector3f _min = Vector3f(std::numeric_limits<float>::max());
    Vector3f _max = Vector3f(std::numeric_limits<float>::lowest());
};
}