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

#ifndef ULTRALISER_MATH_MATRIX4F_H
#define ULTRALISER_MATH_MATRIX4F_H

#include <common/Common.h>

namespace Ultraliser
{

class Matrix2f;
class Matrix3f;
class Quat4f;
class Vector3f;
class Vector4f;

/**
 * @brief The Matrix4f class
 * 4x4 Matrix, stored in column major order (OpenGL style).
 */
class Matrix4f
{
public:

    /**
     * @brief Matrix4f
     * Fill a 4x4 matrix with "fill".  Default to 0.
     * @param fill
     */
    Matrix4f(float fill = 0.f);

    /**
     * @brief Matrix4f
     * @param m00
     * @param m01
     * @param m02
     * @param m03
     * @param m10
     * @param m11
     * @param m12
     * @param m13
     * @param m20
     * @param m21
     * @param m22
     * @param m23
     * @param m30
     * @param m31
     * @param m32
     * @param m33
     */
    Matrix4f(float m00, float m01, float m02, float m03,
             float m10, float m11, float m12, float m13,
             float m20, float m21, float m22, float m23,
             float m30, float m31, float m32, float m33);

    /**
     * @brief Matrix4f
     * setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2 v3]
     * otherwise, sets the rows.
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     * @param setColumns
     */
    Matrix4f(const Vector4f& v0,
             const Vector4f& v1,
             const Vector4f& v2,
             const Vector4f& v3,
             bool setColumns = true);

    /**
     * @brief Matrix4f
     * @param rm
     */
    Matrix4f(const Matrix4f& rm);

    /**
     * @brief operator =
     * @param rm
     * @return
     */
    Matrix4f& operator = (const Matrix4f& rm);

    /**
     * @brief operator /=
     * @param d
     * @return
     */
    Matrix4f& operator/=(float d);

    /**
     * @brief operator ()
     * @param i
     * @param j
     * @return
     */
    const float& operator () (int i, int j) const;

    /**
     * @brief operator ()
     * @param i
     * @param j
     * @return
     */
    float& operator () (int i, int j);

    /**
     * @brief getRow
     * @param i
     * @return
     */
    Vector4f getRow(int i) const;

    /**
     * @brief setRow
     * @param i
     * @param v
     */
    void setRow(int i, const Vector4f& v);

    /**
     * @brief getCol
     * Get column j (mod 4).
     * @param j
     * @return
     */
    Vector4f getCol(int j) const;

    /**
     * @brief setCol
     * @param j
     * @param v
     */
    void setCol(int j, const Vector4f& v);

    /**
     * @brief getSubmatrix2x2
     * Gets the 2x2 submatrix of this matrix to m starting with upper
     * left corner at (i0, j0).
     * @param i0
     * @param j0
     * @return
     */
    Matrix2f getSubmatrix2x2(int i0, int j0) const;

    /**
     * @brief getSubmatrix3x3
     * Gets the 3x3 submatrix of this matrix to m starting with upper left
     * corner at (i0, j0).
     * @param i0
     * @param j0
     * @return
     */
    Matrix3f getSubmatrix3x3(int i0, int j0) const;

    /**
     * @brief setSubmatrix2x2
     * Sets a 2x2 submatrix of this matrix to m starting with upper
     * left corner at (i0, j0).
     * @param i0
     * @param j0
     * @param m
     */
    void setSubmatrix2x2(int i0, int j0, const Matrix2f& m);

    /**
     * @brief setSubmatrix3x3
     * Sets a 3x3 submatrix of this matrix to m starting with upper left
     * corner at (i0, j0).
     * @param i0
     * @param j0
     * @param m
     */
    void setSubmatrix3x3(int i0, int j0, const Matrix3f& m);

    /**
     * @brief determinant
     * @return
     */
    float determinant() const;

    /**
     * @brief inverse
     * @param pbIsSingular
     * @param epsilon
     * @return
     */
    Matrix4f inverse(bool* pbIsSingular = nullptr, float epsilon = 0.f) const;

    /**
     * @brief transpose
     */
    void transpose();

    /**
     * @brief transposed
     * @return
     */
    Matrix4f transposed() const;

public:

    /**
     * @brief operator float *
     *
     */
    operator float* (); //

    /**
     * @brief operator const float *
     * Automatic type conversion for OpenGL.
     */
    operator const float* () const;

    /**
     * @brief print
     */
    void print();

    /**
     * @brief ones
     * @return
     */
    static Matrix4f ones();

    /**
     * @brief identity
     * @return
     */
    static Matrix4f identity();

    /**
     * @brief translation
     * @param x
     * @param y
     * @param z
     * @return
     */
    static Matrix4f translation(float x, float y, float z);

    /**
     * @brief translation
     * @param rTranslation
     * @return
     */
    static Matrix4f translation(const Vector3f& rTranslation);

    /**
     * @brief rotateX
     * @param radians
     * @return
     */
    static Matrix4f rotateX(float radians);

    /**
     * @brief rotateY
     * @param radians
     * @return
     */
    static Matrix4f rotateY(float radians);

    /**
     * @brief rotateZ
     * @param radians
     * @return
     */
    static Matrix4f rotateZ(float radians);

    /**
     * @brief rotation
     * @param rDirection
     * @param radians
     * @return
     */
    static Matrix4f rotation(const Vector3f& rDirection, float radians);

    /**
     * @brief scaling
     * @param sx
     * @param sy
     * @param sz
     * @return
     */
    static Matrix4f scaling(float sx, float sy, float sz);

    /**
     * @brief uniformScaling
     * @param s
     * @return
     */
    static Matrix4f uniformScaling(float s);

    /**
     * @brief lookAt
     * @param eye
     * @param center
     * @param up
     * @return
     */
    static Matrix4f lookAt(const Vector3f& eye,
                           const Vector3f& center,
                           const Vector3f& up);

    /**
     * @brief orthographicProjection
     * @param width
     * @param height
     * @param zNear
     * @param zFar
     * @param directX
     * @return
     */
    static Matrix4f orthographicProjection(float width,
                                           float height,
                                           float zNear,
                                           float zFar,
                                           bool directX);

    /**
     * @brief orthographicProjection
     * @param left
     * @param right
     * @param bottom
     * @param top
     * @param zNear
     * @param zFar
     * @param directX
     * @return
     */
    static Matrix4f orthographicProjection(float left,
                                           float right,
                                           float bottom,
                                           float top,
                                           float zNear,
                                           float zFar,
                                           bool directX);

    /**
     * @brief perspectiveProjection
     * @param fLeft
     * @param fRight
     * @param fBottom
     * @param fTop
     * @param fZNear
     * @param fZFar
     * @param directX
     * @return
     */
    static Matrix4f perspectiveProjection(float fLeft,
                                          float fRight,
                                          float fBottom,
                                          float fTop,
                                          float fZNear,
                                          float fZFar,
                                          bool directX);

    /**
     * @brief perspectiveProjection
     * @param fovYRadians
     * @param aspect
     * @param zNear
     * @param zFar
     * @param directX
     * @return
     */
    static Matrix4f perspectiveProjection(float fovYRadians,
                                          float aspect,
                                          float zNear,
                                          float zFar,
                                          bool directX);

    /**
     * @brief infinitePerspectiveProjection
     * @param fLeft
     * @param fRight
     * @param fBottom
     * @param fTop
     * @param fZNear
     * @param directX
     * @return
     */
    static Matrix4f infinitePerspectiveProjection(float fLeft,
                                                  float fRight,
                                                  float fBottom,
                                                  float fTop,
                                                  float fZNear,
                                                  bool directX);

    /**
     * @brief rotation
     * @param q
     * @return the rotation matrix represented by a quaternion uses a
     * normalized version of q.
     */
    static Matrix4f rotation(const Quat4f& q);

    /**
     * @brief randomRotation
     * @param u0
     * @param u1
     * @param u2
     * @return an orthogonal matrix that's a uniformly distributed rotation
     * given u[i] is a uniformly distributed random number in [0,1].
     */
    static Matrix4f randomRotation(float u0, float u1, float u2);

private:

    /**
     * @brief _elements
     */
    float _elements[ 16 ];

};

/**
 * @brief operator *
 * Matrix-Vector multiplication, 4x4 * 4x1 ==> 4x1.
 * @param m
 * @param v
 * @return
 */
Vector4f operator * (const Matrix4f& m, const Vector4f& v);

/**
 * @brief operator *
 * Matrix-Matrix multiplication.
 * @param x
 * @param y
 * @return
 */
Matrix4f operator * (const Matrix4f& x, const Matrix4f& y);

}

#endif // ULTRALISER_MATH_MATRIX4F_H
