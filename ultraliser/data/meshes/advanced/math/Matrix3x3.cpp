/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#include <math.h>
#include "Matrix3x3.h"

namespace Ultraliser
{

Matrix3x3::Matrix3x3(const double &a11, const double &a12, const double &a13,
                     const double &a21, const double &a22, const double &a23,
                     const double &a31, const double &a32, const double &a33)
{
    elements[0] = a11; elements[1] = a12; elements[2] = a13;
    elements[3] = a21; elements[4] = a22; elements[5] = a23;
    elements[6] = a31; elements[7] = a32; elements[8] = a33;
}

Matrix3x3::Matrix3x3(const double &v1, const double &v2, const double &v3,
                     const double &w1, const double &w2, const double &w3)
{
    elements[0] = v1*w1; elements[1] = v1*w2; elements[2] = v1*w3;
    elements[3] = v2*w1; elements[4] = v2*w2; elements[5] = v2*w3;
    elements[6] = v3*w1; elements[7] = v3*w2; elements[8] = v3*w3;
}

Matrix3x3::Matrix3x3(const double &x, const double &y, const double &z)
{
    elements[0] = x * x;  elements[1] = x * y;  elements[2] = x * z;
    elements[3] = elements[1]; elements[4] = y * y;  elements[5] = y * z;
    elements[6] = elements[2]; elements[7] = elements[5]; elements[8] = z * z;
}

void Matrix3x3::operator+=(const Matrix3x3 &m)
{
    elements[0] += m.elements[0]; elements[1] += m.elements[1]; elements[2] += m.elements[2];
    elements[3] += m.elements[3]; elements[4] += m.elements[4]; elements[5] += m.elements[5];
    elements[6] += m.elements[6]; elements[7] += m.elements[7]; elements[8] += m.elements[8];
}

void Matrix3x3::operator-=(const Matrix3x3 &m)
{
    elements[0] -= m.elements[0]; elements[1] -= m.elements[1]; elements[2] -= m.elements[2];
    elements[3] -= m.elements[3]; elements[4] -= m.elements[4]; elements[5] -= m.elements[5];
    elements[6] -= m.elements[6]; elements[7] -= m.elements[7]; elements[8] -= m.elements[8];
}

void Matrix3x3::operator*=(const double &number)
{
    elements[0] *= number; elements[1] *= number; elements[2] *= number;
    elements[3] *= number; elements[4] *= number; elements[5] *= number;
    elements[6] *= number; elements[7] *= number; elements[8] *= number;
}

Matrix3x3 Matrix3x3::operator+(const Matrix3x3 &m) const
{
    return Matrix3x3(
            elements[0] + m.elements[0], elements[1] + m.elements[1], elements[2] + m.elements[2],
            elements[3] + m.elements[3], elements[4] + m.elements[4], elements[5] + m.elements[5],
            elements[6] + m.elements[6], elements[7] + m.elements[7], elements[8] + m.elements[8]);
}

Matrix3x3 Matrix3x3::operator * (const double &number) const
{
    return Matrix3x3(
                elements[0] * number, elements[1] * number, elements[2] * number,
                elements[3] * number, elements[4] * number, elements[5] * number,
                elements[6] * number, elements[7] * number, elements[8] * number);
}

Matrix3x3 Matrix3x3::operator * (const Matrix3x3& q) const
{
    return Matrix3x3(
                elements[0] * q.elements[0] +
                elements[1] * q.elements[3] +
                elements[2] * q.elements[6],

                elements[0] * q.elements[1] +
                elements[1] * q.elements[4] +
                elements[2] * q.elements[7],

                elements[0] * q.elements[2] +
                elements[1] * q.elements[5] +
                elements[2] * q.elements[8],

                elements[3] * q.elements[0] +
                elements[4] * q.elements[3] +
                elements[5] * q.elements[6],

                elements[3] * q.elements[1] +
                elements[4] * q.elements[4] +
                elements[5] * q.elements[7],

                elements[3] * q.elements[2] +
                elements[4] * q.elements[5] +
                elements[5] * q.elements[8],

                elements[6] * q.elements[0] +
                elements[7] * q.elements[3] +
                elements[8] * q.elements[6],

                elements[6] * q.elements[1] +
                elements[7] * q.elements[4] +
                elements[8] * q.elements[7],

                elements[6] * q.elements[2] +
                elements[7] * q.elements[5] +
                elements[8] * q.elements[8]);
}

Matrix3x3 Matrix3x3::operator~() const
{
    return Matrix3x3(elements[0], elements[3], elements[6],
                     elements[1], elements[4], elements[7],
                     elements[2], elements[5], elements[8]);
}

double Matrix3x3::lrMultiply(const double &x, const double &y, const double &z) const
{
    return (x * (x * elements[0] + y * elements[3] + z * elements[6]) +
            y * (x * elements[1] + y * elements[4] + z * elements[7]) +
            z * (x * elements[2] + y * elements[5] + z * elements[8]));
}

double Matrix3x3::lrMultiply(const double &v1, const double &v2, const double &v3,
                             const double &w1, const double &w2, const double &w3) const
{
    return (w1 * (v1 * elements[0] + v2 * elements[3] + v3 * elements[6]) +
            w2 * (v1 * elements[1] + v2 * elements[4] + v3 * elements[7]) +
            w3 * (v1 * elements[2] + v2 * elements[5] + v3 * elements[8]));
}

Matrix3x3 Matrix3x3::transpose() const
{
    return Matrix3x3(elements[0], elements[3], elements[6],
                     elements[1], elements[4], elements[7],
                     elements[2], elements[5], elements[8]);
}

}
