/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
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

#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The VoxelXYZ struct
 */
template <class T>
struct VoxelXYZ
{
public:

    /**
     * @brief VoxelXYZ
     */
    VoxelXYZ(const T& i, const T& j, const T& k) { x = i; y = j; z = k; }

public:

    /**
     * @brief x
     */
    T x;

    /**
     * @brief y
     */
    T y;

    /**
     * @brief z
     */
    T z;
};

/**
 * @brief VoxelXYZUI8
 */
typedef VoxelXYZ< uint8_t > VoxelXYZUI8;

/**
 * @brief VoxelXYZUI16
 */
typedef VoxelXYZ< uint16_t > VoxelXYZUI16;

/**
 * @brief VoxelXYZUI32
 */
typedef VoxelXYZ< uint32_t > VoxelXYZUI32;

/**
 * @brief VoxelXYZUI64
 */
typedef VoxelXYZ< uint64_t > VoxelXYZUI64;

/**
 * @brief VoxelXYZF32
 */
typedef VoxelXYZ< float > VoxelXYZF32;

/**
 * @brief VoxelXYZF64
 */
typedef VoxelXYZ< double > VoxelXYZF64;

/**
 * @brief VoxelsXYZUI8
 */
typedef std::vector< VoxelXYZUI8 > VoxelsXYZUI8;

/**
 * @brief VoxelsXYZUI16
 */
typedef std::vector< VoxelXYZUI16 > VoxelsXYZUI16;

/**
 * @brief VoxelXYZUI32
 */
typedef std::vector< VoxelXYZUI32 > VoxelsXYZUI32;

/**
 * @brief VoxelsXYZUI64
 */
typedef std::vector< VoxelXYZUI64 > VoxelsXYZUI64;

/**
 * @brief VoxelXYZF32
 */
typedef std::vector< VoxelXYZF32 > VoxelsXYZF32;

/**
 * @brief VoxelsXYZF64
 */
typedef std::vector< VoxelXYZF64 > VoxelsXYZF64;

}
