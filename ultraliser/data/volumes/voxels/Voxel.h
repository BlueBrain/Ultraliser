/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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
 * @brief The Voxel struct
 */
template <class T>
struct Voxel
{
public:

    /**
     * @brief Voxel
     * @param voxelValue
     */
    Voxel(const T voxelValue = 0)
    {
        value = voxelValue;
    }

public:

    /**
     * @brief value
     */
    T value;
};

/**
 * @brief VoxelsUI8
 */
typedef std::vector< Voxel< uint8_t > > VoxelsUI8;

/**
 * @brief VoxelsUI16
 */
typedef std::vector< Voxel< uint16_t > > VoxelsUI16;

/**
 * @brief VoxelsUI32
 */
typedef std::vector< Voxel< uint32_t > > VoxelsUI32;

/**
 * @brief VoxelsUI64
 */
typedef std::vector< Voxel< uint64_t > > VoxelsUI64;

/**
 * @brief VoxelsF32
 */
typedef std::vector< Voxel< float > > VoxelsF32;

/**
 * @brief VoxelsF64
 */
typedef std::vector< Voxel< double > > VoxelsF64;

}
