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
 * @brief The ThinningVoxel struct
 */
template< class T >
struct ThinningVoxel
{
public:

    /**
     * @brief ThinningVoxel
     */
    ThinningVoxel();
    ~ThinningVoxel() { }

    /**
     * @brief ThinningVoxel
     * @param i
     * @param j
     * @param k
     * @param deletable
     */
    ThinningVoxel(const T i, const T j, const T k,
                  const bool& deletable = false,
                  const bool& border = false,
                  const bool& active = true)
    {
        this->x = i; this->y = j; this->z = k;
        this->deletable = deletable;
        this->border = border;
        this->active = active;
    }

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
     * @brief deletable
     */
    bool deletable = false;

    /**
     * @brief border
     */
    bool border = false;

    /**
     * @brief active
     */
    bool active = true;
};

/**
 * @brief ThinningVoxelUI8
 */
typedef ThinningVoxel< uint8_t > ThinningVoxelUI8;

/**
 * @brief ThinningVoxelUI16
 */
typedef ThinningVoxel< uint16_t > ThinningVoxelUI16;

/**
 * @brief ThinningVoxelUI32
 */
typedef ThinningVoxel< uint32_t > ThinningVoxelUI32;

/**
 * @brief ThinningVoxelUI64
 */
typedef ThinningVoxel< uint64_t > ThinningVoxelUI64;

/**
 * @brief ThinningVoxelsUI8
 */
typedef std::vector< ThinningVoxelUI8 > ThinningVoxelsUI8;

/**
 * @brief ThinningVoxelsUI16
 */
typedef std::vector< ThinningVoxelUI16 > ThinningVoxelsUI16;

/**
 * @brief ThinningVoxelsUI32
 */
typedef std::vector< ThinningVoxelUI32 > ThinningVoxelsUI32;

/**
 * @brief ThinningVoxelsUI64
 */
typedef std::vector< ThinningVoxelUI64 > ThinningVoxelsUI64;

/**
 * @brief ThinningVoxelsUI8
 */
typedef std::vector< ThinningVoxelUI8* > ThinningVoxelUI8List;

/**
 * @brief ThinningVoxelsUI16
 */
typedef std::vector< ThinningVoxelUI16* > ThinningVoxelsUI16List;

/**
 * @brief ThinningVoxelsUI32
 */
typedef std::vector< ThinningVoxelUI32* > ThinningVoxelsUI32List;

/**
 * @brief ThinningVoxelsUI64
 */
typedef std::vector< ThinningVoxelUI64* > ThinningVoxelsUI64List;

}
