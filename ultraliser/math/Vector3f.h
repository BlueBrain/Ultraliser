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

namespace Ultraliser
{

class Vector2f;

/**
 * @brief The Vector3f class
 */
class Vector3f
{
public:

    /**
     * @brief Vector3f
     * @param f
     */
    Vector3f(float f = 0.f);

    /**
     * @brief Vector3f
     * @param x
     * @param y
     * @param z
     */
    Vector3f(float x, float y, float z);

    /**
     * @brief Vector3f
     * @param xy
     * @param z
     */
    Vector3f(const Vector2f& xy, float z);

    /**
     * @brief Vector3f
     * @param x
     * @param yz
     */
    Vector3f(float x, const Vector2f& yz);

    /**
     * @brief Vector3f
     * @param rv
     */
    Vector3f(const Vector3f& rv);

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
     * @brief z
     * @return
     */
    float& z();

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
     * @brief z
     * @return
     */
    float z() const;

    /**
     * @brief xy
     * @return
     */
    Vector2f xy() const;

    /**
     * @brief xz
     * @return
     */
    Vector2f xz() const;

    /**
     * @brief yz
     * @return
     */
    Vector2f yz() const;

    /**
     * @brief xyz
     * @return
     */
    Vector3f xyz() const;

    /**
     * @brief yzx
     * @return
     */
    Vector3f yzx() const;

    /**
     * @brief zxy
     * @return
     */
    Vector3f zxy() const;

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
    Vector3f normalized() const;

    /**
     * @brief homogenized
     * @return
     */
    Vector2f homogenized() const;

    /**
     * @brief negate
     */
    void negate();

    /**
     * @brief getLargestDimension
     * Gets the largest dimension in the vector.
     *
     * @return
     * The largest dimension in the vector
     */
    float getLargestDimension() const;

    /**
     * @brief operator const float *
     * Automatic type conversion for OpenGL.
     */
    operator const float* () const;

    /**
     * @brief operator float *
     * Automatic type conversion for OpenGL.
     */
    operator float* ();

    /**
     * @brief print
     */
    void print() const;

    /**
     * @brief isNan
     */
    bool isNan() const;

    /**
     * @brief isZero
     * @return
     */
    bool isZero() const;
    
    /**
     * @brief orthogonal
     * Gets an orthogonal vector.
     * @return
     */
    Vector3f orthogonal() const;

public:

    /**
     * @brief operator =
     * @param rv
     * @return
     */
    Vector3f& operator = (const Vector3f& rv);


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
    Vector3f& operator += (const Vector3f& v);

    /**
     * @brief operator -=
     * @param v
     * @return
     */
    Vector3f& operator -= (const Vector3f& v);

    /**
     * @brief operator *=
     * @param f
     * @return
     */
    Vector3f& operator *= (float f);

    /**
     * @brief operator /=
     * @param f
     * @return
     */
    Vector3f& operator /= (float f);

    const float* data() const
    {
        return _elements;
    }

public:

    /**
     * @brief dot
     * @param v0
     * @param v1
     * @return
     */
    static float dot(const Vector3f& v0, const Vector3f& v1);

    /**
     * @brief cross
     * @param v0
     * @param v1
     * @return
     */
    static Vector3f cross(const Vector3f& v0, const Vector3f& v1);
    
    /**
     * @brief lerp
     * Computes the linear interpolation between v0 and v1 by alpha \in [0,1].
     * @param v0
     * @param v1
     * @param alpha
     * @return v0 * (1 - alpha) * v1 * alpha
     */
    static Vector3f lerp(const Vector3f& v0, const Vector3f& v1, float alpha);

    /**
     * @brief cubicInterpolate
     * Computes the cubic catmull-rom interpolation between p0, p1, p2, p3
     * by t \in [0,1].
     * Guarantees that at t = 0, the result is p0 and at p1, the result is p2.
     * @param p0
     * @param p1
     * @param p2
     * @param p3
     * @param t
     * @return
     */
    static Vector3f cubicInterpolate(const Vector3f& p0,
                                      const Vector3f& p1,
                                      const Vector3f& p2,
                                      const Vector3f& p3,
                                      float t);

public:

    /**
     * @brief ZERO
     */
    static const Vector3f ZERO;

    /**
     * @brief ONE
     */
    static const Vector3f ONE;

    /**
     * @brief UP
     */
    static const Vector3f UP;

    /**
     * @brief RIGHT
     */
    static const Vector3f RIGHT;

    /**
     * @brief FORWARD
     */
    static const Vector3f FORWARD;

private:

    /**
     * @brief _elements
     */
    float _elements[ 3 ];

};

/**
 * @brief operator +
 * @param v0
 * @param v1
 * @return
 */
Vector3f operator + (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator -
 * @param v0
 * @param v1
 * @return
 */
Vector3f operator - (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator *
 * @param v0
 * @param v1
 * @return
 */
Vector3f operator * (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator /
 * @param v0
 * @param v1
 * @return
 */
Vector3f operator / (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator -
 * @param v
 * @return
 */
Vector3f operator - (const Vector3f& v);

/**
 * @brief operator *
 * @param f
 * @param v
 * @return
 */
Vector3f operator * (float f, const Vector3f& v);

/**
 * @brief operator *
 * @param v
 * @param f
 * @return
 */
Vector3f operator * (const Vector3f& v, float f);

/**
 * @brief operator /
 * @param v
 * @param f
 * @return
 */
Vector3f operator / (const Vector3f& v, float f);

/**
 * @brief operator ==
 * @param v0
 * @param v1
 * @return
 */
bool operator == (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator >
 * @param v0
 * @param v1
 * @return
 */
bool operator > (const Vector3f& v0, const Vector3f& v1);

/**
 * @brief operator <
 * @param v0
 * @param v1
 * @return
 */
bool operator < (const Vector3f& v0, const Vector3f& v1);


/**
 * @brief operator !=
 * @param v0
 * @param v1
 * @return
 */
bool operator != (const Vector3f& v0, const Vector3f& v1);

}
