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

#include "Matrix3f.h"
#include "Matrix2f.h"
#include "Quat4f.h"
#include "Vector3f.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

Matrix3f::Matrix3f(float fill)
{
    for (int i = 0; i < 9; ++i)
    {
        _elements[i] = fill;
    }
}


Matrix3f::Matrix3f(float m00, float m01, float m02,
                   float m10, float m11, float m12,
                   float m20, float m21, float m22)
{
    _elements[0] = m00;
    _elements[1] = m10;
    _elements[2] = m20;

    _elements[3] = m01;
    _elements[4] = m11;
    _elements[5] = m21;

    _elements[6] = m02;
    _elements[7] = m12;
    _elements[8] = m22;
}

Matrix3f::Matrix3f(const Vector3f& v0,
                   const Vector3f& v1,
                   const Vector3f& v2,
                   bool setColumns)
{
    if (setColumns)
    {
        setCol(0, v0);
        setCol(1, v1);
        setCol(2, v2);
    }
    else
    {
        setRow(0, v0);
        setRow(1, v1);
        setRow(2, v2);
    }
}

Matrix3f::Matrix3f(const Matrix3f& rm)
{
    memcpy(_elements, rm._elements, sizeof(_elements));
}

Matrix3f& Matrix3f::operator = (const Matrix3f& rm)
{
    if (this != &rm)
    {
        memcpy(_elements, rm._elements, sizeof(_elements));
    }
    return *this;
}

const float& Matrix3f::operator () (int i, int j) const
{
    return _elements[j * 3 + i];
}

float& Matrix3f::operator () (int i, int j)
{
    return _elements[j * 3 + i];
}

Vector3f Matrix3f::getRow(int i) const
{
    return Vector3f(_elements[i], _elements[i + 3], _elements[i + 6]);
}

void Matrix3f::setRow(int i, const Vector3f& v)
{
    _elements[i] = v.x();
    _elements[i + 3] = v.y();
    _elements[i + 6] = v.z();
}

Vector3f Matrix3f::getCol(int j) const
{
    int colStart = 3 * j;

    return Vector3f(_elements[colStart],
                    _elements[colStart + 1],
                    _elements[colStart + 2]);
}

void Matrix3f::setCol(int j, const Vector3f& v)
{
    int colStart = 3 * j;

    _elements[colStart] = v.x();
    _elements[colStart + 1] = v.y();
    _elements[colStart + 2] = v.z();
}

Matrix2f Matrix3f::getSubmatrix2x2(int i0, int j0) const
{
    Matrix2f out;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            out(i, j) = (*this)(i + i0, j + j0);
        }
    }

    return out;
}

void Matrix3f::setSubmatrix2x2(int i0, int j0, const Matrix2f& m)
{
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            (*this)(i + i0, j + j0) = m(i, j);
        }
    }
}

float Matrix3f::determinant() const
{
    return Matrix3f::determinant3x3
            (_elements[0], _elements[3], _elements[6],
            _elements[1], _elements[4], _elements[7],
            _elements[2], _elements[5], _elements[8]);
}

Matrix3f Matrix3f::inverse(bool* pbIsSingular, float epsilon) const
{
    float m00 = _elements[0];
    float m10 = _elements[1];
    float m20 = _elements[2];

    float m01 = _elements[3];
    float m11 = _elements[4];
    float m21 = _elements[5];

    float m02 = _elements[6];
    float m12 = _elements[7];
    float m22 = _elements[8];

    float cofactor00 =  Matrix2f::determinant2x2(m11, m12, m21, m22);
    float cofactor01 = -Matrix2f::determinant2x2(m10, m12, m20, m22);
    float cofactor02 =  Matrix2f::determinant2x2(m10, m11, m20, m21);

    float cofactor10 = -Matrix2f::determinant2x2(m01, m02, m21, m22);
    float cofactor11 =  Matrix2f::determinant2x2(m00, m02, m20, m22);
    float cofactor12 = -Matrix2f::determinant2x2(m00, m01, m20, m21);

    float cofactor20 =  Matrix2f::determinant2x2(m01, m02, m11, m12);
    float cofactor21 = -Matrix2f::determinant2x2(m00, m02, m10, m12);
    float cofactor22 =  Matrix2f::determinant2x2(m00, m01, m10, m11);

    float determinant = m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02;

    bool isSingular = (std::fabs(determinant) < epsilon);
    if (isSingular)
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = true;
        }
        return Matrix3f();
    }
    else
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = false;
        }

        float reciprocalDeterminant = 1.0f / determinant;

        return Matrix3f(cofactor00 * reciprocalDeterminant,
                        cofactor10 * reciprocalDeterminant,
                        cofactor20 * reciprocalDeterminant,
                        cofactor01 * reciprocalDeterminant,
                        cofactor11 * reciprocalDeterminant,
                        cofactor21 * reciprocalDeterminant,
                        cofactor02 * reciprocalDeterminant,
                        cofactor12 * reciprocalDeterminant,
                        cofactor22 * reciprocalDeterminant);
    }
}

void Matrix3f::transpose()
{
    float temp;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = i + 1; j < 3; ++j)
        {
            temp = (*this)(i, j);
            (*this)(i, j) = (*this)(j, i);
            (*this)(j, i) = temp;
        }
    }
}

Matrix3f Matrix3f::transposed() const
{
    Matrix3f out;
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            out(j, i) = (*this)(i, j);
        }
    }

    return out;
}

Matrix3f::operator float* ()
{
    return _elements;
}

void Matrix3f::print()
{
    printf("[%.4f %.4f %.4f]\n[%.4f %.4f %.4f]\n[%.4f %.4f %.4f]\n",
            F2D(_elements[0]), F2D(_elements[3]), F2D(_elements[6]),
            F2D(_elements[1]), F2D(_elements[4]), F2D(_elements[7]),
            F2D(_elements[2]), F2D(_elements[5]), F2D(_elements[8]));
}

float Matrix3f::determinant3x3(float m00, float m01, float m02,
                                float m10, float m11, float m12,
                                float m20, float m21, float m22)
{
    return(m00 * (m11 * m22 - m12 * m21)
         - m01 * (m10 * m22 - m12 * m20)
         + m02 * (m10 * m21 - m11 * m20));
}

Matrix3f Matrix3f::ones()
{
    Matrix3f m;
    for (int i = 0; i < 9; ++i)
    {
        m._elements[i] = 1;
    }

    return m;
}

Matrix3f Matrix3f::identity()
{
    Matrix3f m;

    m(0, 0) = 1;
    m(1, 1) = 1;
    m(2, 2) = 1;

    return m;
}

Matrix3f Matrix3f::rotateX(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(1, 0, 0, 0, c, -s, 0, s, c);
}

Matrix3f Matrix3f::rotateY(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(c, 0, s, 0, 1, 0, -s, 0, c);
}

Matrix3f Matrix3f::rotateZ(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix3f(c, -s, 0, s, c, 0, 0, 0, 1);
}

Matrix3f Matrix3f::scaling(float sx, float sy, float sz)
{
    return Matrix3f(sx, 0, 0, 0, sy, 0, 0, 0, sz);
}

Matrix3f Matrix3f::uniformScaling(float s)
{
    return Matrix3f(s, 0, 0, 0, s, 0, 0, 0, s);
}

Matrix3f Matrix3f::rotation(const Vector3f& rDirection, float radians)
{
    Vector3f normalizedDirection = rDirection.normalized();

    float cosTheta = cos(radians);
    float sinTheta = sin(radians);

    float x = normalizedDirection.x();
    float y = normalizedDirection.y();
    float z = normalizedDirection.z();

    return Matrix3f(x * x * (1.0f - cosTheta) + cosTheta,
                     y * x * (1.0f - cosTheta) - z * sinTheta,
                     z * x * (1.0f - cosTheta) + y * sinTheta,
                     x * y * (1.0f - cosTheta) + z * sinTheta,
                     y * y * (1.0f - cosTheta) + cosTheta,
                     z * y * (1.0f - cosTheta) - x * sinTheta,
                     x * z * (1.0f - cosTheta) - y * sinTheta,
                     y * z * (1.0f - cosTheta) + x * sinTheta,
                     z * z * (1.0f - cosTheta) + cosTheta);
}

Matrix3f Matrix3f::rotation(const Quat4f& rq)
{
    Quat4f q = rq.normalized();

    float xx = q.x() * q.x();
    float yy = q.y() * q.y();
    float zz = q.z() * q.z();

    float xy = q.x() * q.y();
    float zw = q.z() * q.w();

    float xz = q.x() * q.z();
    float yw = q.y() * q.w();

    float yz = q.y() * q.z();
    float xw = q.x() * q.w();

    return Matrix3f(1.0f - 2.0f * (yy + zz),
                     2.0f * (xy - zw),
                     2.0f * (xz + yw),
                     2.0f * (xy + zw),
                     1.0f - 2.0f * (xx + zz),
                     2.0f * (yz - xw),
                     2.0f * (xz - yw),
                     2.0f * (yz + xw),
                     1.0f - 2.0f * (xx + yy));
}

Vector3f operator * (const Matrix3f& m, const Vector3f& v)
{
    Vector3f output(0, 0, 0);

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            output[i] += m(i, j) * v[j];
        }
    }

    return output;
}

Matrix3f operator * (const Matrix3f& x, const Matrix3f& y)
{
    Matrix3f product; // Zeroes

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            for (int k = 0; k < 3; ++k)
            {
                product(i, k) += x(i, j) * y(j, k);
            }
        }
    }

    return product;
}

}
