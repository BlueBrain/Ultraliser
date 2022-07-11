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

#include "Vector2f.h"
#include "Matrix2f.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

Matrix2f::Matrix2f(float fill)
{
    for (int i = 0; i < 4; ++i)
    {
        _elements[i] = fill;
    }
}

Matrix2f::Matrix2f(float m00, float m01,
                   float m10, float m11)
{
    _elements[0] = m00;
    _elements[1] = m10;

    _elements[2] = m01;
    _elements[3] = m11;
}

Matrix2f::Matrix2f(const Vector2f& v0, const Vector2f& v1, bool setColumns)
{
    if (setColumns)
    {
        setCol(0, v0);
        setCol(1, v1);
    }
    else
    {
        setRow(0, v0);
        setRow(1, v1);
    }
}

Matrix2f::Matrix2f(const Matrix2f& rm)
{
    memcpy(_elements, rm._elements, sizeof(_elements));
}

Matrix2f& Matrix2f::operator = (const Matrix2f& rm)
{
    if (this != &rm)
    {
        memcpy(_elements, rm._elements, sizeof(_elements) );
    }
    return *this;
}

const float& Matrix2f::operator () (int i, int j) const
{
    return _elements[j * 2 + i];
}

float& Matrix2f::operator () (int i, int j)
{
    return _elements[j * 2 + i];
}

Vector2f Matrix2f::getRow(int i) const
{
    return Vector2f(_elements[i], _elements[i + 2]);
}

void Matrix2f::setRow(int i, const Vector2f& v)
{
    _elements[i] = v.x();
    _elements[i + 2] = v.y();
}

Vector2f Matrix2f::getCol(int j) const
{
    int colStart = 2 * j;

    return Vector2f(_elements[colStart], _elements[colStart + 1]);
}

void Matrix2f::setCol(int j, const Vector2f& v)
{
    int colStart = 2 * j;

    _elements[colStart] = v.x();
    _elements[colStart + 1] = v.y();
}

float Matrix2f::determinant()
{
    return Matrix2f::determinant2x2(_elements[0], _elements[2],
            _elements[1], _elements[3]);
}

Matrix2f Matrix2f::inverse(bool* pbIsSingular, float epsilon)
{
    float determinant = _elements[0] * _elements[3] -
                        _elements[2] * _elements[1];

    bool isSingular = (std::fabs(determinant) < epsilon);
    if (isSingular)
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = true;
        }
        return Matrix2f();
    }
    else
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = false;
        }

        float reciprocalDeterminant = 1.0f / determinant;

        return Matrix2f (_elements[3] * reciprocalDeterminant,
                -_elements[2] * reciprocalDeterminant,
                -_elements[1] * reciprocalDeterminant,
                _elements[0] * reciprocalDeterminant);
    }
}

void Matrix2f::transpose()
{
    float m01 = (*this)(0, 1);
    float m10 = (*this)(1, 0);

    (*this)(0, 1) = m10;
    (*this)(1, 0) = m01;
}

Matrix2f Matrix2f::transposed() const
{
    return Matrix2f((*this)(0, 0), (*this)(1, 0),
                    (*this)(0, 1), (*this)(1, 1));

}

Matrix2f::operator float* ()
{
    return _elements;
}

void Matrix2f::print()
{
    printf("[%.4f %.4f]\n[%.4f %.4f]\n", F2D(_elements[0]), F2D(_elements[2]),
                                         F2D(_elements[1]), F2D(_elements[3]));
}

float Matrix2f::determinant2x2(float m00, float m01,
                               float m10, float m11)
{
    return(m00 * m11 - m01 * m10);
}

Matrix2f Matrix2f::ones()
{
    Matrix2f m;
    for (int i = 0; i < 4; ++i)
    {
        m._elements[i] = 1;
    }

    return m;
}

Matrix2f Matrix2f::identity()
{
    Matrix2f m;

    m(0, 0) = 1;
    m(1, 1) = 1;

    return m;
}

Matrix2f Matrix2f::rotation(float degrees)
{
    float c = cos(degrees);
    float s = sin(degrees);

    return Matrix2f(c, -s, s, c);
}

Matrix2f operator * (float f, const Matrix2f& m)
{
    Matrix2f output;

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            output(i, j) = f * m(i, j);
        }
    }

    return output;
}

Matrix2f operator * (const Matrix2f& m, float f)
{
    return f * m;
}

Vector2f operator * (const Matrix2f& m, const Vector2f& v)
{
    Vector2f output(0, 0);

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            output[i] += m(i, j) * v[j];
        }
    }

    return output;
}

Matrix2f operator * (const Matrix2f& x, const Matrix2f& y)
{
    Matrix2f product; // Zeroes

    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            for (int k = 0; k < 2; ++k)
            {
                product(i, k) += x(i, j) * y(j, k);
            }
        }
    }

    return product;
}

}
