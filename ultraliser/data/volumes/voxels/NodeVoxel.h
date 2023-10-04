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
 * @brief The NodeVoxel struct
 */
template< class T >
struct NodeVoxel
{

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

    /**
     * @brief index
     */
    size_t index;
};

/**
 * @brief NodeVoxelUI8
 */
typedef NodeVoxel< uint8_t > NodeVoxelUI8;

/**
 * @brief NodeVoxelUI16
 */
typedef NodeVoxel< uint16_t > NodeVoxelUI16;

/**
 * @brief NodeVoxelUI32
 */
typedef NodeVoxel< uint32_t > NodeVoxelUI32;

/**
 * @brief NodeVoxelUI64
 */
typedef NodeVoxel< uint64_t > NodeVoxelUI64;

/**
 * @brief NodeVoxelsUI8
 */
typedef std::vector< NodeVoxelUI8 > NodeVoxelsUI8;

/**
 * @brief NodeVoxelsUI16
 */
typedef std::vector< NodeVoxelUI16 > NodeVoxelsUI16;

/**
 * @brief NodeVoxelsUI32
 */
typedef std::vector< NodeVoxelUI32 > NodeVoxelsUI32;

/**
 * @brief NodeVoxelsUI64
 */
typedef std::vector< NodeVoxelUI64 > NodeVoxelsUI64;

}
