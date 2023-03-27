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

#include "SparseOctreeNodeBounds.h"

#include "SparseOctreeNodeSlot.h"

namespace Ultraliser
{
Bounds SparseOctreeNodeBounds::fromParentBounds(const Bounds &parentBounds, uint8_t slotMask)
{
    auto left = slotMask & SparseOctreeNodeSlot::backBottomLeft || slotMask & SparseOctreeNodeSlot::backTopLeft
        || slotMask & SparseOctreeNodeSlot::frontBottomLeft || slotMask & SparseOctreeNodeSlot::frontTopLeft;

    auto bottom = slotMask & SparseOctreeNodeSlot::backBottomLeft || slotMask & SparseOctreeNodeSlot::backBottomRight
        || slotMask & SparseOctreeNodeSlot::frontBottomLeft || slotMask & SparseOctreeNodeSlot::frontBottomRight;

    auto back = slotMask < SparseOctreeNodeSlot::frontBottomLeft;

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
}