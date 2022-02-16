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

#ifndef ULTRALISER_MATH_QUAT4F_H
#define ULTRALISER_MATH_QUAT4F_H

#include <math/Matrix3f.h>

namespace Ultraliser
{

class Vector3f;
class Vector4f;

/**
 * @brief The Quat4f class
 */
class Quat4f
{
public:

    /**
     * @brief Quat4f
     */
    Quat4f();

    /**
     * @brief Quat4f
     * Q = w + x * i + y * j + z * k
     * @param w
     * @param x
     * @param y
     * @param z
     */
    Quat4f(float w, float x, float y, float z);

    /**
     * @brief Quat4f
     * @param rq
     */
    Quat4f(const Quat4f& rq);

    /**
     * @brief operator =
     * @param rq
     * @return
     */
    Quat4f& operator = (const Quat4f& rq);

    /**
     * @brief Quat4f
     * Returns a quaternion with 0 real part.
     * @param v
     */
    Quat4f(const Vector3f& v);

    /**
     * @brief Quat4f
     * Copies the components of a Vector4f directly into this quaternion.
     * @param v
     */
    Quat4f(const Vector4f& v);

    /**
     * @brief operator []
     * @param i
     * @return the ith element.
     */
    const float& operator [] (int i) const;

    /**
     * @brief operator []
     * @param i
     * @return
     */
    float& operator [] (int i);

    /**
     * @brief w
     * @return
     */
    float w() const;

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
     * @brief xyz
     * @return
     */
    Vector3f xyz() const;

    /**
     * @brief wxyz
     * @return
     */
    Vector4f wxyz() const;

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
    Quat4f normalized() const;

    /**
     * @brief conjugate
     */
    void conjugate();

    /**
     * @brief conjugated
     * @return
     */
    Quat4f conjugated() const;

    /**
     * @brief invert
     */
    void invert();

    /**
     * @brief inverse
     * @return
     */
    Quat4f inverse() const;

    /**
     * @brief log
     * @return
     */
    Quat4f log() const;

    /**
     * @brief exp
     * @return
     */
    Quat4f exp() const;

    /**
     * @brief rotate
     * Rotate vector v by quaternion
     * @return
     */
    Vector3f rotate(const Vector3f& v) const;

    /**
     * @brief getAxisAngle
     * @param radiansOut
     * @return unit vector for rotation and radians about the unit vector.
     */
    Vector3f getAxisAngle(float* radiansOut);

    /**
     * @brief setAxisAngle
     * Sets this quaternion to be a rotation of fRadians about v = < fx, fy, fz >,
     * v need not necessarily be unit length
     * @param radians
     * @param axis
     */
    void setAxisAngle(float radians, const Vector3f& axis);

    /**
     * @brief print
     */
    void print();

    /**
     * @brief dot
     * Quaternion dot product (a la vector)
     * @param q0
     * @param q1
     * @return
     */
    static float dot(const Quat4f& q0, const Quat4f& q1);

    /**
     * @brief lerp
     * Linear (stupid) interpolation
     * @param q0
     * @param q1
     * @param alpha
     * @return
     */
    static Quat4f lerp(const Quat4f& q0, const Quat4f& q1, float alpha);

    /**
     * @brief slerp
     * Spherical linear interpolation
     * @param a
     * @param b
     * @param t
     * @param allowFlip
     * @return
     */
    static Quat4f slerp(const Quat4f& a, const Quat4f& b, float t, bool allowFlip = true);

    /**
     * @brief squad
     * Spherical quadratic interoplation between a and b at point t  given
     * quaternion tangents tanA and tanB (can be computed using squadTangent).
     * @param a
     * @param tanA
     * @param tanB
     * @param b
     * @param t
     * @return
     */
    static Quat4f squad(const Quat4f& a,
                         const Quat4f& tanA,
                         const Quat4f& tanB,
                         const Quat4f& b,
                         float t);

    /**
     * @brief cubicInterpolate
     * @param q0
     * @param q1
     * @param q2
     * @param q3
     * @param t
     * @return
     */
    static Quat4f cubicInterpolate(const Quat4f& q0,
                                    const Quat4f& q1,
                                    const Quat4f& q2,
                                    const Quat4f& q3,
                                    float t);

    /**
     * @brief logDifference
     * Log-difference between a and b, used for squadTangent
     * @param a
     * @param b
     * @return log(a^-1 b)
     */
    static Quat4f logDifference(const Quat4f& a, const Quat4f& b);

    /**
     * @brief squadTangent
     * Computes a tangent at center, defined by the before and after quaternions
     * Useful for squad()
     * @param before
     * @param center
     * @param after
     * @return
     */
    static Quat4f squadTangent(const Quat4f& before,
                                const Quat4f& center,
                                const Quat4f& after);

    /**
     * @brief fromRotationMatrix
     * @param m
     * @return
     */
    static Quat4f fromRotationMatrix(const Matrix3f& m);

    /**
     * @brief fromRotatedBasis
     * @param x
     * @param y
     * @param z
     * @return
     */
    static Quat4f fromRotatedBasis(const Vector3f& x,
                                    const Vector3f& y,
                                    const Vector3f& z);

    /**
     * @brief fromTwoVectors
     * @param v0
     * @param v1
     * @return
     */
    static Quat4f fromTwoVectors(const Vector3f& v0,
                                 const Vector3f& v1);

    /**
     * @brief randomRotation
     * @param u0
     * @param u1
     * @param u2
     * @return a unit quaternion that's a uniformly distributed rotation
     * given u[i] is a uniformly distributed random number in [0,1].
     * Taken from Graphics Gems II.
     */
    static Quat4f randomRotation(float u0, float u1, float u2);

public:

    /**
     * @brief ZERO
     */
    static const Quat4f ZERO;

    /**
     * @brief IDENTITY
     */
    static const Quat4f IDENTITY;

private:

    /**
     * @brief _elements
     */
    float _elements[ 4 ];
};

/**
 * @brief operator +
 * @param q0
 * @param q1
 * @return
 */
Quat4f operator + (const Quat4f& q0, const Quat4f& q1);

/**
 * @brief operator -
 * @param q0
 * @param q1
 * @return
 */
Quat4f operator - (const Quat4f& q0, const Quat4f& q1);

/**
 * @brief operator *
 * @param q0
 * @param q1
 * @return
 */
Quat4f operator * (const Quat4f& q0, const Quat4f& q1);

/**
 * @brief operator *
 * @param f
 * @param q
 * @return
 */
Quat4f operator * (float f, const Quat4f& q);

/**
 * @brief operator *
 * @param q
 * @param f
 * @return
 */
Quat4f operator * (const Quat4f& q, float f);

}
#endif // ULTRALISER_MATH_QUAT4F_H
