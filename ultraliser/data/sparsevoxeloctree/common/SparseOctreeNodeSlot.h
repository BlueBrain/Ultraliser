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

#include <cstdint>

namespace Ultraliser
{
struct SparseOctreeNodeSlot
{
    static constexpr uint8_t root = 0;
    static constexpr uint8_t backBottomLeft = 1;
    static constexpr uint8_t backBottomRight = 1 << 1;
    static constexpr uint8_t backTopLeft = 1 << 2;
    static constexpr uint8_t backTopRight = 1 << 3;
    static constexpr uint8_t frontBottomLeft = 1 << 4;
    static constexpr uint8_t frontBottomRight = 1 << 5;
    static constexpr uint8_t frontTopLeft = 1 << 6;
    static constexpr uint8_t frontTopRight = 1 << 7;
};
}