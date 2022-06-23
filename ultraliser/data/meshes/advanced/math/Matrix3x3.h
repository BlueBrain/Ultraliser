/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marco Attene < IMATI-GE / CNR >
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
 *
 * The content of this file is based on MeshFix. The code has been modified under the terms of
 * the GNU General Public License as published by the Free Software Foundation either version 3 of
 * the License, or (at your option) any later version.
 * MeshFix has a dual license for free and commercial use. For further information, please refer
 * to the original repository at < https://github.com/MarcoAttene/MeshFix-V2.1>.
 **************************************************************************************************/

#pragma once

namespace Ultraliser
{

/**
 * @brief The Matrix3x3 class
 * Generic 3x3 matrix.
 * Elements are stored in a row-dominant order, thus for example, M[4] is the first element of the
 * second row.
 */
class Matrix3x3
{
public:

    /**
     * @brief elements
     * Actual values of the matrix.
     */
    double elements[9];

    /**
     * @brief Matrix3x3
     * Contructs a null matrix.
     */
    Matrix3x3()
    {
        elements[0] = elements[1] = elements[2] = 0.0;
        elements[3] = elements[4] = elements[5] = 0.0;
        elements[6] = elements[7] = elements[8] = 0.0;
    }

    /**
     * @brief Matrix3x3
     * Constructs a fully initialized matrix.
     * @param a11
     * @param a12
     * @param a13
     * @param a21
     * @param a22
     * @param a23
     * @param a31
     * @param a32
     * @param a33
     */
    Matrix3x3(const double &a11, const double &a12, const double &a13,
              const double &a21, const double &a22, const double &a23,
              const double &a31, const double &a32, const double &a33);

    /**
     * @brief Matrix3x3
     * Constructs a 3x3 matrix as the product of Transpose(v1, v2, v3) and (w1, w2, w3).
     * @param v1
     * @param v2
     * @param v3
     * @param w1
     * @param w2
     * @param w3
     */
    Matrix3x3(const double &v1, const double &v2, const double &v3,
              const double &w1, const double &w2, const double &w3);

    /**
     * @brief Matrix3x3
     * Constructs a 3x3 matrix as the product of Transpose(a, b, c) and (a, b, c).
     * @param a
     * @param b
     * @param c
     */
    Matrix3x3(const double &a, const double &b, const double &c);

    /**
     * @brief isSymmetric
     * Returns TRUE if the matrix is symmetric.
     * @return
     */
    bool isSymmetric() const
    {
        return (elements[2] == elements[4] &&
                elements[3] == elements[7] &&
                elements[6] == elements[8]);
    }

    /**
     * @brief operator =
     * Initializes all elements to the given number.
     * @param d
     */
    void operator=(const double &number)
    {
        elements[0] = elements[1] = elements[2] = number;
        elements[3] = elements[4] = elements[5] = number;
        elements[6] = elements[7] = elements[8] = number;
    }

    /**
     * @brief operator +=
     * Sum another matrix.
     */
    void operator+=(const Matrix3x3 &m);

    /**
     * @brief operator -=
     * Subtract another matrix.
     */
    void operator-=(const Matrix3x3 &m);

    /**
     * @brief operator *=
     * Multiply by a scalar.
     */
    void operator*=(const double &number);

    /**
     * @brief operator /=
     * Divide by a scalar.
     * @param d
     */
    void operator/=(const double &number) { operator *= (1.0 / number); }

    /**
     * @brief operator +
     * Returns the sum of this and another matrix.
     * @return
     */
    Matrix3x3 operator+(const Matrix3x3 &m) const;

    /**
     * @brief operator *
     * Returns the product of this matrix with a scalar.
     * @return
     */
    Matrix3x3 operator*(const double &number) const;

    /**
     * @brief operator *
     * Returns the product of this and another matrix (rows by columns).
     * @return
     */
    Matrix3x3 operator*(const Matrix3x3 &m) const;

    /**
     * @brief operator ~
     * Returns the transpose of this matrix.
     * @return
     */
    Matrix3x3 operator~() const;

    /**
     * @brief transpose
     * Return the matrix transpose.
     * @return
     */
    Matrix3x3 transpose() const;

    /**
     * @brief lrMultiply
     * Returns Transpose(a, b, c) * M * (a, b, c)
     * Returns the (scalar) result of multiplying the matrix on the left and
     * on the right by the vector (a, b, c).
     * @param a
     * @param b
     * @param c
     * @return
     */
    double lrMultiply(const double &a, const double &b, const double &c) const;

    /**#ifndef ULTRALISER_DATA_STRUCTURES_MATRIX_3X3_H
#define ULTRALISER_DATA_STRUCTURES_MATRIX_3X3_H

     * @brief lrMultiply
     * Returns the (scalar) result of v * M * w.
     * @param v1
     * @param v2
     * @param v3
     * @param w1
     * @param w2
     * @param w3
     * @return
     */
    double lrMultiply(const double &v1, const double &v2, const double &v3,
                      const double &w1, const double &w2, const double &w3) const;
};

}

