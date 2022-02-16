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

#ifndef ULTRALISER_MATH_VECTOR_4F_H
#define ULTRALISER_MATH_VECTOR_4F_H

namespace Ultraliser
{

class Vector2f;
class Vector3f;

/**
 * @brief The Vector4f class
 */
class Vector4f
{
public:

    /**
     * @brief Vector4f
     * @param f
     */
    Vector4f(float f = 0.f);

    /**
     * @brief Vector4f
     * @param fx
     * @param fy
     * @param fz
     * @param fw
     */
    Vector4f(float fx, float fy, float fz, float fw);

    /**
     * @brief Vector4f
     * @param buffer
     */
    Vector4f(float buffer[ 4 ]);

    /**
     * @brief Vector4f
     * @param xy
     * @param z
     * @param w
     */
    Vector4f(const Vector2f& xy, float z, float w);

    /**
     * @brief Vector4f
     * @param x
     * @param yz
     * @param w
     */
    Vector4f(float x, const Vector2f& yz, float w);

    /**
     * @brief Vector4f
     * @param x
     * @param y
     * @param zw
     */
    Vector4f(float x, float y, const Vector2f& zw);

    /**
     * @brief Vector4f
     * @param xy
     * @param zw
     */
    Vector4f(const Vector2f& xy, const Vector2f& zw);

    /**
     * @brief Vector4f
     * @param xyz
     * @param w
     */
    Vector4f(const Vector3f& xyz, float w);

    /**
     * @brief Vector4f
     * @param x
     * @param yzw
     */
    Vector4f(float x, const Vector3f& yzw);

    /**
     * @brief Vector4f
     * @param rv
     */
    Vector4f(const Vector4f& rv);

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
     * @brief w
     * @return
     */
    float& w();

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
     * @brief w
     * @return
     */
    float w() const;

    /**
     * @brief xy
     * @return
     */
    Vector2f xy() const;

    /**
     * @brief yz
     * @return
     */
    Vector2f yz() const;

    /**
     * @brief zw
     * @return
     */
    Vector2f zw() const;

    /**
     * @brief wx
     * @return
     */
    Vector2f wx() const;

    /**
     * @brief xyz
     * @return
     */
    Vector3f xyz() const;

    /**
     * @brief yzw
     * @return
     */
    Vector3f yzw() const;

    /**
     * @brief zwx
     * @return
     */
    Vector3f zwx() const;

    /**
     * @brief wxy
     * @return
     */
    Vector3f wxy() const;

    /**
     * @brief xyw
     * @return
     */
    Vector3f xyw() const;

    /**
     * @brief yzx
     * @return
     */
    Vector3f yzx() const;

    /**
     * @brief zwy
     * @return
     */
    Vector3f zwy() const;

    /**
     * @brief wxz
     * @return
     */
    Vector3f wxz() const;

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
    Vector4f normalized() const;

    /**
     * @brief homogenize
     * If v.z != 0, v = v / v.w
     */
    void homogenize();

    /**
     * @brief homogenized
     * @return
     */
    Vector4f homogenized() const;

    /**
     * @brief negate
     */
    void negate();

public:

    /**
     * @brief operator =
     * @param rv
     * @return
     */
    Vector4f& operator = (const Vector4f& rv);

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

public:

    /**
     * @brief dot
     * @param v0
     * @param v1
     * @return
     */
    static float dot(const Vector4f& v0, const Vector4f& v1);

    /**
     * @brief lerp
     * @param v0
     * @param v1
     * @param alpha
     * @return
     */
    static Vector4f lerp(const Vector4f& v0, const Vector4f& v1, float alpha);

private:

    /**
     * @brief _elements
     */
    float _elements[4];

};

/**
 * @brief operator +
 * @param v0
 * @param v1
 * @return
 */
Vector4f operator + (const Vector4f& v0, const Vector4f& v1);

/**
 * @brief operator -
 * @param v0
 * @param v1
 * @return
 */
Vector4f operator - (const Vector4f& v0, const Vector4f& v1);

/**
 * @brief operator *
 * @param v0
 * @param v1
 * @return
 */
Vector4f operator * (const Vector4f& v0, const Vector4f& v1);

/**
 * @brief operator /
 * @param v0
 * @param v1
 * @return
 */
Vector4f operator / (const Vector4f& v0, const Vector4f& v1);

/**
 * @brief operator -
 * @param v
 * @return
 */
Vector4f operator - (const Vector4f& v);

/**
 * @brief operator *
 * @param f
 * @param v
 * @return
 */
Vector4f operator * (float f, const Vector4f& v);

/**
 * @brief operator *
 * @param v
 * @param f
 * @return
 */
Vector4f operator * (const Vector4f& v, float f);

/**
 * @brief operator /
 * @param v
 * @param f
 * @return
 */
Vector4f operator / (const Vector4f& v, float f);

/**
 * @brief operator ==
 * @param v0
 * @param v1
 * @return
 */
bool operator == (const Vector4f& v0, const Vector4f& v1);

/**
 * @brief operator !=
 * @param v0
 * @param v1
 * @return
 */
bool operator != (const Vector4f& v0, const Vector4f& v1);

}

#endif // ULTRALISER_MATH_VECTOR_4F_H
