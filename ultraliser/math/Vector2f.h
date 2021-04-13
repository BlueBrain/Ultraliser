/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_MATH_VECTOR_2F_H
#define ULTRALISER_MATH_VECTOR_2F_H

#include <cmath>

namespace Ultraliser
{

class Vector3f;

/**
 * @brief The Vector2f class
 */
class Vector2f
{
public:
    
    /**
     * @brief ZERO
     */
    static const Vector2f ZERO;

    /**
     * @brief UP
     */
    static const Vector2f UP;

    /**
     * @brief RIGHT
     */
    static const Vector2f RIGHT;

    /**
     * @brief LEFT
     */
    static const Vector2f LEFT;

    /**
     * @brief DOWN
     */
    static const Vector2f DOWN;

public:

    /**
     * @brief Vector2f
     * @param f
     */
    Vector2f(float f = 0.f);

    /**
     * @brief Vector2f
     * @param x
     * @param y
     */
    Vector2f(float x, float y);

    /**
     * @brief Vector2f
     * @param rv
     */
    Vector2f(const Vector2f& rv);


public:

    /**
     * @brief x
     * @return
     */
    float& x();

    /**
     * @brief y
     * @return
     */
    float& y();

    /**
     * @brief x
     * @return
     */
    float x() const;

    /**
     * @brief y
     * @return
     */
    float y() const;

    /**
     * @brief xy
     * @return
     */
    Vector2f xy() const;

    /**
     * @brief yx
     * @return
     */
    Vector2f yx() const;

    /**
     * @brief xx
     * @return
     */
    Vector2f xx() const;

    /**
     * @brief yy
     * @return
     */
    Vector2f yy() const;

    /**
     * @brief normal
     * @return (-y, x)
     */
    Vector2f normal() const;

    /**
     * @brief abs
     * @return
     */
    float abs() const;

    /**
     * @brief absSquared
     * @return
     */
    float absSquared() const;

    /**
     * @brief normalize
     */
    void normalize();

    /**
     * @brief normalized
     * @return
     */
    Vector2f normalized() const;

    /**
     * @brief negate
     */
    void negate();

    /**
     * @brief operator const float*
     * Automatic type conversion for OpenGL.
     */
    operator const float* () const;

    /**
     * @brief operator float*
     * Automatic type conversion for OpenGL.
     */
    operator float* ();

    /**
     * @brief print
     */
    void print() const;

public:

    /**
     * @brief operator =
     * @param rv
     * @return
     */
    Vector2f& operator = (const Vector2f& rv);

    /**
     * @brief operator []
     * @param i
     * @return
     */
    const float& operator [] (int i) const;

    /**
     * @brief operator []
     * @param i
     * @return
     */
    float& operator [] (int i);

    /**
     * @brief operator +=
     * @param v
     * @return
     */
    Vector2f& operator += (const Vector2f& v);

    /**
     * @brief operator -=
     * @param v
     * @return
     */
    Vector2f& operator -= (const Vector2f& v);

    /**
     * @brief operator *=
     * @param f
     * @return
     */
    Vector2f& operator *= (float f);

public:

    /**
     * @brief dot
     * @param v0
     * @param v1
     * @return
     */
    static float dot(const Vector2f& v0, const Vector2f& v1);

    /**
     * @brief cross
     * @param v0
     * @param v1
     * @return
     */
    static Vector3f cross(const Vector2f& v0, const Vector2f& v1);

    /**
     * @brief lerp
     * @param v0
     * @param v1
     * @param alpha
     * @return v0 * (1 - alpha) * v1 * alpha
     */
    static Vector2f lerp(const Vector2f& v0, const Vector2f& v1, float alpha);

private:

    /**
     * @brief _elements
     */
    float _elements[2];
};

/**
 * @brief operator +
 * @param v0
 * @param v1
 * @return
 */
Vector2f operator + (const Vector2f& v0, const Vector2f& v1);

/**
 * @brief operator -
 * @param v0
 * @param v1
 * @return
 */
Vector2f operator - (const Vector2f& v0, const Vector2f& v1);

/**
 * @brief operator *
 * @param v0
 * @param v1
 * @return
 */
Vector2f operator * (const Vector2f& v0, const Vector2f& v1);

/**
 * @brief operator /
 * @param v0
 * @param v1
 * @return
 */
Vector2f operator / (const Vector2f& v0, const Vector2f& v1);

/**
 * @brief operator -
 * Unary negation
 * @param v
 * @return
 */
Vector2f operator - (const Vector2f& v);

/**
 * @brief operator *
 * @param f
 * @param v
 * @return
 */
Vector2f operator * (float f, const Vector2f& v);

/**
 * @brief operator *
 * @param v
 * @param f
 * @return
 */
Vector2f operator * (const Vector2f& v, float f);

/**
 * @brief operator /
 * @param v
 * @param f
 * @return
 */
Vector2f operator / (const Vector2f& v, float f);

/**
 * @brief operator ==
 * @param v0
 * @param v1
 * @return
 */
bool operator == (const Vector2f& v0, const Vector2f& v1);

/**
 * @brief operator !=
 * @param v0
 * @param v1
 * @return
 */
bool operator != (const Vector2f& v0, const Vector2f& v1);

}

#endif // ULTRALISER_MATH_VECTOR_2F_H
