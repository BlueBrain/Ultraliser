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

#ifndef ULTRALISER_DATA_STRUCTURES_MATRIX_4X4_H
#define ULTRALISER_DATA_STRUCTURES_MATRIX_4X4_H

#include <stdio.h>
#include <float.h>
#include <data/structures/List.h>

namespace Ultraliser
{

/**
 * @brief The Matrix4x4 class
 * Generic 4x4 matrix.
 */
class Matrix4x4
{
public:

    /**
     * @brief elements
     * Actual matrix coefficients
     */
    double elements[4][4];

    /**
     * @brief Matrix4x4
     * Constructs a diagonal matrix with @number values on the diagonal.
     *
     * @param number
     * Initialization value, default 0.0.
     */
    Matrix4x4(const double& number=0.f);


    /**
     * @brief Matrix4x4
     * Constructs a fully initialized matrix (parameters are in row dominant order .
     * Example: M[0][0], M[0][1], ...).
     *
     * @param a11
     * @param a12
     * @param a13
     * @param a14
     * @param a21
     * @param a22
     * @param a23
     * @param a24
     * @param a31
     * @param a32
     * @param a33
     * @param a34
     * @param a41
     * @param a42
     * @param a43
     * @param a44
     */
    Matrix4x4(const double &a11, const double &a12, const double &a13, const double &a14,
              const double &a21, const double &a22, const double &a23, const double &a24,
              const double &a31, const double &a32, const double &a33, const double &a34,
              const double &a41, const double &a42, const double &a43, const double &a44);

    /**
     * @brief setRotation
     * otation matrix from a quaternion.
     *
     * @param rx
     * @param ry
     * @param rz
     * @param rw
     */
    void setRotation(const double &rx, const double &ry, const double &rz, const double &rw);

    /**
     * @brief setTranslation
     * Translation matrix from a vector.
     *
     * @param x
     * The X-coordinate of the vector.
     * @param y
     * The Y-coordinate of the vector.
     * @param z
     * The Z-coordinate of the vector.
     */
    void setTranslation(const double &x, const double &y, const double &z);

    /**
     * @brief operator *
     * Matrix multiplication operator.
     *
     * @param m
     * Input matrix
     * @return
     * Returns the product of this and another matrix (rows by columns).
     */
    Matrix4x4 operator*(const Matrix4x4& m) const;

    /**
     * @brief transform
     * Transform the input vector by left-multiplication with the matrix.
     *
     * @param x
     * The X-coordinate of the vector.
     * @param y
     * The Y-coordinate of the vector.
     * @param z
     * The Z-coordinate of the vector.
     */
    void transform(double *x, double *y, double *z);
};

}

#endif // ULTRALISER_DATA_STRUCTURES_MATRIX_4X4_H

