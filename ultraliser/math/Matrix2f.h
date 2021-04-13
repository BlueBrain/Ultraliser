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

#ifndef ULTRALISER_MATH_MATRIX2F_H
#define ULTRALISER_MATH_MATRIX2F_H

#include <cstdio>

namespace Ultraliser
{

class Vector2f;

/**
 * @brief The Matrix2f class
 * 2x2 Matrix, stored in column major order (OpenGL style).
 */
class Matrix2f
{
public:

    /**
     * @brief Matrix2f
     * Fill a 2x2 matrix with "fill", default to 0.
     * @param fill
     */
    Matrix2f(float fill = 0.f);

    /**
     * @brief Matrix2f
     * @param m00
     * @param m01
     * @param m10
     * @param m11
     */
    Matrix2f(float m00, float m01,
              float m10, float m11);

    /**
     * @brief Matrix2f
     * SetColumns = true ==> sets the columns of the matrix to be [v0 v1]
     * otherwise, sets the rows
     * @param v0
     * @param v1
     * @param setColumns
     */
    Matrix2f(const Vector2f& v0, const Vector2f& v1, bool setColumns = true);

    /**
     * @brief Matrix2f
     * @param rm
     */
    Matrix2f(const Matrix2f& rm);

public:

    /**
     * @brief getRow
     * @param i
     * @return
     */
    Vector2f getRow(int i) const;

    /**
     * @brief setRow
     * @param i
     * @param v
     */
    void setRow(int i, const Vector2f& v);

    /**
     * @brief getCol
     * @param j
     * @return
     */
    Vector2f getCol(int j) const;

    /**
     * @brief setCol
     * @param j
     * @param v
     */
    void setCol(int j, const Vector2f& v);

    /**
     * @brief determinant
     * @return
     */
    float determinant();

    /**
     * @brief inverse
     * @param pbIsSingular
     * @param epsilon
     * @return
     */
    Matrix2f inverse(bool* pbIsSingular = nullptr, float epsilon = 0.f);

    /**
     * @brief transpose
     */
    void transpose();

    /**
     * @brief transposed
     * @return
     */
    Matrix2f transposed() const;

    /**
     * @brief operator float *
     * Automatic type conversion for OpenGL.
     */
    operator float* ();

    /**
     * @brief print
     */
    void print();

    /**
     * @brief determinant2x2
     * @param m00
     * @param m01
     * @param m10
     * @param m11
     * @return
     */
    static float determinant2x2(float m00, float m01,
                                 float m10, float m11);

    /**
     * @brief ones
     * @return
     */
    static Matrix2f ones();

    /**
     * @brief identity
     * @return
     */
    static Matrix2f identity();

    /**
     * @brief rotation
     * @param degrees
     * @return
     */
    static Matrix2f rotation(float degrees);

public:

    /**
     * @brief operator =
     * @param rm
     * @return
     */
    Matrix2f& operator = (const Matrix2f& rm);

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

private:

    /**
     * @brief _elements
     */
    float _elements[ 4 ];
};

/**
 * @brief operator *
 * @param f
 * @param m
 * @return
 */
Matrix2f operator * (float f, const Matrix2f& m);

/**
 * @brief operator *
 * @param m
 * @param f
 * @return
 */
Matrix2f operator * (const Matrix2f& m, float f);

/**
 * @brief operator *
 * Matrix-Vector multiplication, 2x2 * 2x1 ==> 2x1
 * @param m
 * @param v
 * @return
 */
Vector2f operator * (const Matrix2f& m, const Vector2f& v);

/**
 * @brief operator *
 * Matrix-Matrix multiplication.
 * @param x
 * @param y
 * @return
 */
Matrix2f operator * (const Matrix2f& x, const Matrix2f& y);

}

#endif // ULTRALISER_MATH_MATRIX2F_H
