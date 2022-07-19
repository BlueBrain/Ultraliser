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

#include <math.h>
#include "Matrix.h"

namespace Ultraliser
{

Matrix4x4::Matrix4x4(const double &number)
{
    elements[0][0] = elements[1][1] = elements[2][2] = elements[3][3] = number;
    elements[1][0] = elements[0][1] = elements[1][2] = elements[1][3] = 0;
    elements[2][0] = elements[2][1] = elements[0][2] = elements[2][3] = 0;
    elements[3][0] = elements[3][1] = elements[3][2] = elements[0][3] = 0;
}

Matrix4x4::Matrix4x4(const double &a11, const double &a12, const double &a13, const double &a14,
                     const double &a21, const double &a22, const double &a23, const double &a24,
                     const double &a31, const double &a32, const double &a33, const double &a34,
                     const double &a41, const double &a42, const double &a43, const double &a44)
{
    elements[0][0] = a11; elements[0][1] = a12; elements[0][2] = a13; elements[0][3] = a14;
    elements[1][0] = a21; elements[1][1] = a22; elements[1][2] = a23; elements[1][3] = a24;
    elements[2][0] = a31; elements[2][1] = a32; elements[2][2] = a33; elements[2][3] = a34;
    elements[3][0] = a41; elements[3][1] = a42; elements[3][2] = a43; elements[3][3] = a44;
}

void Matrix4x4::setRotation(const double &rx, const double &ry, const double &rz, const double &rw)
{
    elements[0][0] = rw * rw + rx * rx - ry * ry - rz * rz;
    elements[0][1] = 2 * rx * ry + 2 * rw * rz;
    elements[0][2] = 2 * rx * rz - 2 * rw * ry;
    elements[0][3] = 0.f;

    elements[1][0] = 2 * rx * ry-2 * rw * rz;
    elements[1][1] = rw * rw - rx * rx + ry * ry - rz * rz;
    elements[1][2] = 2 * ry * rz + 2 * rw * rx;
    elements[1][3] = 0.f;

    elements[2][0] = 2 * rx * rz + 2 * rw * ry;
    elements[2][1] = 2 * ry * rz - 2 * rw * rx;
    elements[2][2] = rw * rw - rx * rx - ry * ry + rz * rz;
    elements[2][3] = 0.f;

    elements[3][0] = 0.f;
    elements[3][1] = 0.f;
    elements[3][2] = 0.f;
    elements[3][3] = rw * rw + rx * rx + ry * ry + rz * rz;
}

void Matrix4x4::setTranslation(const double &x, const double &y, const double &z)
{
    elements[0][0] = elements[1][1] = elements[2][2] = elements[3][3] = 1;
    elements[1][0] = elements[0][1] = elements[1][2] = 0;
    elements[2][0] = elements[2][1] = elements[0][2] = 0;
    elements[3][0] = elements[3][1] = elements[3][2] = 0;
    elements[0][3] = x;
    elements[1][3] = y;
    elements[2][3] = z;
}

Matrix4x4 Matrix4x4::operator*(const Matrix4x4& q) const
{
    Matrix4x4 m;

    for (size_t i = 0; i < 4; ++i)
    {
        for (size_t j = 0; j < 4; ++j)
        {
            m.elements[i][j] = elements[i][0] * q.elements[0][j] +
                               elements[i][1] * q.elements[1][j] +
                               elements[i][2] * q.elements[2][j] +
                               elements[i][3] * q.elements[3][j];
        }
    }

    return m;
}

void Matrix4x4::transform(double *x, double *y, double *z)
{
    double w, a = *x, b = *y, c = *z;

    *x = elements[0][0] * a + elements[0][1] * b + elements[0][2] * c + elements[0][3];
    *y = elements[1][0] * a + elements[1][1] * b + elements[1][2] * c + elements[1][3];
    *z = elements[2][0] * a + elements[2][1] * b + elements[2][2] * c + elements[2][3];
    w  = elements[3][0] * a + elements[3][1] * b + elements[3][2] * c + elements[3][3];

    (*x) /= w; (*y) /= w; (*z) /= w;
}

}
