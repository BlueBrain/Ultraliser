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
#include "Vector3f.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

const Vector2f Vector2f::ZERO       = Vector2f(0.f, 0.f);
const Vector2f Vector2f::UP         = Vector2f(0.f, 1.f);
const Vector2f Vector2f::RIGHT      = Vector2f(1.f, 0.f);
const Vector2f Vector2f::LEFT       = Vector2f(-1.f, 0.f);
const Vector2f Vector2f::DOWN       = Vector2f(0.f, -1.f);



Vector2f::Vector2f(float f)
{
    _elements[0] = f;
    _elements[1] = f;
}

Vector2f::Vector2f(float x, float y)
{
    _elements[0] = x;
    _elements[1] = y;
}

Vector2f::Vector2f(const Vector2f& rv)
{
    _elements[0] = rv[0];
    _elements[1] = rv[1];
}

Vector2f& Vector2f::operator = (const Vector2f& rv)
{
    if (this != &rv)
    {
        _elements[0] = rv[0];
        _elements[1] = rv[1];
    }
    return *this;
}

const float& Vector2f::operator [] (int i) const
{
    return _elements[i];
}

float& Vector2f::operator [] (int i)
{
    return _elements[i];
}

float& Vector2f::x()
{
    return _elements[0];
}

float& Vector2f::y()
{
    return _elements[1];
}

float Vector2f::x() const
{
    return _elements[0];
}

float Vector2f::y() const
{
    return _elements[1];
}

Vector2f Vector2f::xy() const
{
    return *this;
}

Vector2f Vector2f::yx() const
{
    return Vector2f(_elements[1], _elements[0]);
}

Vector2f Vector2f::xx() const
{
    return Vector2f(_elements[0], _elements[0]);
}

Vector2f Vector2f::yy() const
{
    return Vector2f(_elements[1], _elements[1]);
}

Vector2f Vector2f::normal() const
{
    return Vector2f(-_elements[1], _elements[0]);
}

float Vector2f::abs() const
{
    return sqrt(absSquared());
}

float Vector2f::absSquared() const
{
    return _elements[0] * _elements[0] + _elements[1] * _elements[1];
}

void Vector2f::normalize()
{
    float norm = abs();
    _elements[0] /= norm;
    _elements[1] /= norm;
}

Vector2f Vector2f::normalized() const
{
    float norm = abs();
    return Vector2f(_elements[0] / norm, _elements[1] / norm);
}

void Vector2f::negate()
{
    _elements[0] = -_elements[0];
    _elements[1] = -_elements[1];
}

Vector2f::operator const float* () const
{
    return _elements;
}

Vector2f::operator float* ()
{
    return _elements;
}

void Vector2f::print() const
{
    printf("< %.4f, %.4f >\n", F2D(_elements[0]), F2D(_elements[1]));
}

Vector2f& Vector2f::operator += (const Vector2f& v)
{
    _elements[0] += v._elements[0];
    _elements[1] += v._elements[1];
    return *this;
}

Vector2f& Vector2f::operator -= (const Vector2f& v)
{
    _elements[0] -= v._elements[0];
    _elements[1] -= v._elements[1];
    return *this;
}

Vector2f& Vector2f::operator *= (float f)
{
    _elements[0] *= f;
    _elements[1] *= f;
    return *this;
}

float Vector2f::dot(const Vector2f& v0, const Vector2f& v1)
{
    return v0[0] * v1[0] + v0[1] * v1[1];
}

Vector3f Vector2f::cross(const Vector2f& v0, const Vector2f& v1)
{
    return Vector3f(0, 0, v0.x() * v1.y() - v0.y() * v1.x());
}

Vector2f Vector2f::lerp(const Vector2f& v0, const Vector2f& v1, float alpha)
{
    return alpha * (v1 - v0) + v0;
}

Vector2f operator + (const Vector2f& v0, const Vector2f& v1)
{
    return Vector2f(v0.x() + v1.x(), v0.y() + v1.y());
}

Vector2f operator - (const Vector2f& v0, const Vector2f& v1)
{
    return Vector2f(v0.x() - v1.x(), v0.y() - v1.y());
}

Vector2f operator * (const Vector2f& v0, const Vector2f& v1)
{
    return Vector2f(v0.x() * v1.x(), v0.y() * v1.y());
}

Vector2f operator / (const Vector2f& v0, const Vector2f& v1)
{
    return Vector2f(v0.x() * v1.x(), v0.y() * v1.y());
}

Vector2f operator - (const Vector2f& v)
{
    return Vector2f(-v.x(), -v.y());
}

Vector2f operator * (float f, const Vector2f& v)
{
    return Vector2f(f * v.x(), f * v.y());
}

Vector2f operator * (const Vector2f& v, float f)
{
    return Vector2f(f * v.x(), f * v.y());
}

Vector2f operator / (const Vector2f& v, float f)
{
    return Vector2f(v.x() / f, v.y() / f);
}

bool operator == (const Vector2f& v0, const Vector2f& v1)
{
    return(isEqual(v0.x(), v1.x()) && isEqual(v0.y(), v1.y()));
}

bool operator != (const Vector2f& v0, const Vector2f& v1)
{
    return !(v0 == v1);
}

}
