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

#ifndef ULTRALISER_MATH_MATRIX3F_H
#define ULTRALISER_MATH_MATRIX3F_H

#include <cstdio>

namespace Ultraliser
{

class Matrix2f;
class Quat4f;
class Vector3f;

/**
 * @brief The Matrix3f class
 * 3x3 Matrix, stored in column major order (OpenGL style).
 */
class Matrix3f
{
public:

    /**
     * @brief Matrix3f
     * Fill a 3x3 matrix with "fill", default to 0.
     * @param fill
     */
    Matrix3f(float fill = 0.f);

    /**
     * @brief Matrix3f
     * @param m00
     * @param m01
     * @param m02
     * @param m10
     * @param m11
     * @param m12
     * @param m20
     * @param m21
     * @param m22
     */
    Matrix3f(float m00, float m01, float m02,
              float m10, float m11, float m12,
              float m20, float m21, float m22);

    /**
     * @brief Matrix3f
     * setColumns = true ==> sets the columns of the matrix to be [v0 v1 v2]
     * otherwise, sets the rows
     * @param v0
     * @param v1
     * @param v2
     * @param setColumns
     */
    Matrix3f(const Vector3f& v0,
              const Vector3f& v1,
              const Vector3f& v2,
              bool setColumns = true);

    /**
     * @brief Matrix3f
     * @param rm
     */
    Matrix3f(const Matrix3f& rm);

public:

    /**
     * @brief getRow
     * @param i
     * @return
     */
    Vector3f getRow(int i) const;

    /**
     * @brief setRow
     * @param i
     * @param v
     */
    void setRow(int i, const Vector3f& v);

    /**
     * @brief getCol
     * @param j
     * @return
     */
    Vector3f getCol(int j) const;

    /**
     * @brief setCol
     * @param j
     * @param v
     */
    void setCol(int j, const Vector3f& v);

    /**
     * @brief getSubmatrix2x2
     * Gets the 2x2 submatrix of this matrix to m starting with upper left
     * corner at (i0, j0)
     * @param i0
     * @param j0
     * @return
     */
    Matrix2f getSubmatrix2x2(int i0, int j0) const;

    /**
     * @brief setSubmatrix2x2
     * Sets a 2x2 submatrix of this matrix to m starting with upper left corner
     * at (i0, j0).
     * @param i0
     * @param j0
     * @param m
     */
    void setSubmatrix2x2(int i0, int j0, const Matrix2f& m);

    /**
     * @brief determinant
     * @return
     */
    float determinant() const;

    /**
     * @brief inverse
     * @param pbIsSingular
     * @param epsilon
     * @todo invert in place as well
     * @return
     */
    Matrix3f inverse(bool* pbIsSingular = nullptr, float epsilon = 0.f) const;

    /**
     * @brief transpose
     */
    void transpose();

    /**
     * @brief transposed
     * @return
     */
    Matrix3f transposed() const;

public:

    /**
     * @brief operator =
     * @param rm
     * @return
     */
    Matrix3f& operator = (const Matrix3f& rm);

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
     * @brief operator float *
     * Automatic type conversion for OpenGL.
     */
    operator float* ();

    /**
     * @brief print
     */
    void print();

public:

    /**
     * @brief determinant3x3
     * @param m00
     * @param m01
     * @param m02
     * @param m10
     * @param m11
     * @param m12
     * @param m20
     * @param m21
     * @param m22
     * @return
     */
    static float determinant3x3(float m00, float m01, float m02,
                                 float m10, float m11, float m12,
                                 float m20, float m21, float m22);


    /**
     * @brief ones
     * @return
     */
    static Matrix3f ones();

    /**
     * @brief identity
     * @return
     */
    static Matrix3f identity();

    /**
     * @brief rotateX
     * @param radians
     * @return
     */
    static Matrix3f rotateX(float radians);

    /**
     * @brief rotateY
     * @param radians
     * @return
     */
    static Matrix3f rotateY(float radians);

    /**
     * @brief rotateZ
     * @param radians
     * @return
     */
    static Matrix3f rotateZ(float radians);

    /**
     * @brief scaling
     * @param sx
     * @param sy
     * @param sz
     * @return
     */
    static Matrix3f scaling(float sx, float sy, float sz);

    /**
     * @brief uniformScaling
     * @param s
     * @return
     */
    static Matrix3f uniformScaling(float s);

    /**
     * @brief rotation
     * @param rDirection
     * @param radians
     * @return
     */
    static Matrix3f rotation(const Vector3f& rDirection, float radians);

    /**
     * @brief rotation
     * @param rq
     * @return the rotation matrix represented by a unit quaternion if q is
     * not normalized, it it normalized first.
     */
    static Matrix3f rotation(const Quat4f& rq);

private:

    /**
     * @brief _elements
     */
    float _elements[9];

};

/**
 * @brief operator *
 * Matrix-Vector multiplication, 3x3 * 3x1 ==> 3x1
 * @param m
 * @param v
 * @return
 */
Vector3f operator * (const Matrix3f& m, const Vector3f& v);

/**
 * @brief operator *
 * Matrix-Matrix multiplication
 * @param x
 * @param y
 * @return
 */
Matrix3f operator * (const Matrix3f& x, const Matrix3f& y);

}

#endif // ULTRALISER_MATH_MATRIX3F_H
