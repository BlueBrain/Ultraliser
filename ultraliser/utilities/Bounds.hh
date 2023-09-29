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

#include <common/Headers.hh>

namespace Ultraliser
{

/**
 * @brief The Bounds3D struct
 */
template < typename T >
struct Bounds3D
{
public:

    /**
     * @brief Bounds3D
     * @param x1
     * @param x2
     * @param y1
     * @param y2
     * @param z1
     * @param z2
     */
    Bounds3D(const T& x1, const T& x2, const T& y1, const T& y2, const T& z1, const T& z2)
    {
        if (x1 < x2) { this->x1 = x1; this->x2 = x2; }
        else { this->x1 = x2; this->x2 = x1; }

        if (y1 < y2) { this->y1 = y1; this->y2 = y2; }
        else { this->y1 = y2; this->y2 = y1; }

        if (z1 < z2) { this->z1 = z1; this->z2 = z2; }
        else { this->z1 = z2; this->z2 = z1; }
    }

    /**
     * @brief getWidth
     * @return
     */
    T getWidth() const { return this->x2 - this->x1; }

    /**
     * @brief getHeight
     * @return
     */
    T getHeight() const { return this->y2 - this->y1; }

    /**
     * @brief getDepth
     * @return
     */
    T getDepth() const { return this->z2 - this->z1; }

    /**
     * @brief getDiagonal
     * @return
     */
    T getDiagonal() const
    {
        const auto& width = getWidth();
        const auto& height = getHeight();
        const auto& depth = getDepth();
        return std::sqrt((width * width) + (height * height) + (depth * depth));
    }

public:

    /**
     * @brief x1
     */
    T x1;

    /**
     * @brief x2
     */
    T x2;

    /**
     * @brief y1
     */
    T y1;

    /**
     * @brief y2
     */
    T y2;

    /**
     * @brief z1
     */
    T z1;

    /**
     * @brief z2
     */
    T z2;
};

typedef struct Bounds3D< int8_t > Bounds3D_i8;
typedef struct Bounds3D< int16_t > Bounds3D_i16;
typedef struct Bounds3D< int32_t > Bounds3D_i32;
typedef struct Bounds3D< int64_t > Bounds3D_i64;

typedef struct Bounds3D< uint8_t > Bounds3D_ui8;
typedef struct Bounds3D< uint16_t > Bounds3D_ui16;
typedef struct Bounds3D< uint32_t > Bounds3D_ui32;
typedef struct Bounds3D< uint64_t > Bounds3D_ui64;

typedef struct Bounds3D< size_t > Bounds3D_size_t;

typedef struct Bounds3D< float > Bounds3D_f;
typedef struct Bounds3D< double > Bounds3D_d;

}
