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

#include "Vector2f.h"
#include "Vector3f.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

const Vector3f Vector3f::ZERO       = Vector3f(0.f, 0.f, 0.f);
const Vector3f Vector3f::UP         = Vector3f(0.f, 1.f, 0.f);
const Vector3f Vector3f::RIGHT      = Vector3f(1.f, 0.f, 0.f);
const Vector3f Vector3f::FORWARD    = Vector3f(0.f, 0.f, -1.f);

Vector3f::Vector3f(float f)
{
    _elements[0] = f;
    _elements[1] = f;
    _elements[2] = f;
}

Vector3f::Vector3f(float x, float y, float z)
{
    _elements[0] = x;
    _elements[1] = y;
    _elements[2] = z;
}

Vector3f::Vector3f(const Vector2f& xy, float z)
{
    _elements[0] = xy.x();
    _elements[1] = xy.y();
    _elements[2] = z;
}

Vector3f::Vector3f(float x, const Vector2f& yz)
{
    _elements[0] = x;
    _elements[1] = yz.x();
    _elements[2] = yz.y();
}

Vector3f::Vector3f(const Vector3f& rv)
{
    _elements[0] = rv[0];
    _elements[1] = rv[1];
    _elements[2] = rv[2];
}

Vector3f& Vector3f::operator = (const Vector3f& rv)
{
    if (this != &rv)
    {
        _elements[0] = rv[0];
        _elements[1] = rv[1];
        _elements[2] = rv[2];
    }
    return *this;
}

const float& Vector3f::operator [] (int i) const
{
    return _elements[i];
}

float& Vector3f::operator [] (int i)
{
    return _elements[i];
}

float& Vector3f::x()
{
    return _elements[0];
}

float& Vector3f::y()
{
    return _elements[1];
}

float& Vector3f::z()
{
    return _elements[2];
}

float Vector3f::x() const
{
    return _elements[0];
}

float Vector3f::y() const
{
    return _elements[1];
}

float Vector3f::z() const
{
    return _elements[2];
}

Vector2f Vector3f::xy() const
{
    return Vector2f(_elements[0], _elements[1]);
}

Vector2f Vector3f::xz() const
{
    return Vector2f(_elements[0], _elements[2]);
}

Vector2f Vector3f::yz() const
{
    return Vector2f(_elements[1], _elements[2]);
}

Vector3f Vector3f::xyz() const
{
    return Vector3f(_elements[0], _elements[1], _elements[2]);
}

Vector3f Vector3f::yzx() const
{
    return Vector3f(_elements[1], _elements[2], _elements[0]);
}

Vector3f Vector3f::zxy() const
{
    return Vector3f(_elements[2], _elements[0], _elements[1]);
}

float Vector3f::distance(const Vector3f& v) const
{
    return sqrt((v._elements[0] - _elements[0]) * (v._elements[0] - _elements[0]) +
                (v._elements[1] - _elements[1]) * (v._elements[1] - _elements[1]) +
                (v._elements[2] - _elements[2]) * (v._elements[2] - _elements[2]));
}

float Vector3f::abs() const
{
    return sqrt(_elements[0] * _elements[0] +
                _elements[1] * _elements[1] +
                _elements[2] * _elements[2]);
}

float Vector3f::absSquared() const
{
    return(_elements[0] * _elements[0] +
           _elements[1] * _elements[1] +
           _elements[2] * _elements[2]);
}

void Vector3f::normalize()
{
    float norm = abs();
    _elements[0] /= norm;
    _elements[1] /= norm;
    _elements[2] /= norm;
}

Vector3f Vector3f::normalized() const
{
    float norm = abs();
    return Vector3f (_elements[0] / norm,
                     _elements[1] / norm,
                     _elements[2] / norm);
}

Vector2f Vector3f::homogenized() const
{
    return Vector2f (_elements[0] / _elements[2], _elements[1] / _elements[2]);
}

void Vector3f::negate()
{
    _elements[0] = -_elements[0];
    _elements[1] = -_elements[1];
    _elements[2] = -_elements[2];
}

float Vector3f::getLargestDimension() const
{
    // Initially, start with X
    float largestDimension = _elements[0];

    // Go with Y
    if (_elements[1] > largestDimension)
        largestDimension = _elements[1];

    // Go with Z
    if (_elements[2] > largestDimension)
        largestDimension = _elements[2];

    // Return the largest dimension
    return largestDimension;
}

Vector3f::operator const float* () const
{
    return _elements;
}

Vector3f::operator float* ()
{
    return _elements;
}

void Vector3f::print() const
{
    printf("< %.4f, %.4f, %.4f >\n", F2D(_elements[0]),
                                     F2D(_elements[1]),
                                     F2D(_elements[2]));
}

bool Vector3f::isNan() const
{
    return std::isnan(_elements[0]) ||
           std::isnan(_elements[1]) ||
           std::isnan(_elements[2]);
}

Vector3f Vector3f::orthogonal() const
{
    // Normalize
    Vector3f normalized = this->normalized();

    // Absolutes
    const float x = std::abs(normalized.x());
    const float y = std::abs(normalized.y());
    const float z = std::abs(normalized.z());

    // Based on the components of the vector, get the other vector
    Vector3f other;
    if (x < y)
    {
        if (x < z)
        {
            other = RIGHT;
        }
        else
        {
            other = FORWARD;
        }
    }
    else
    {
        if (y < z)
        {
            other = UP;
        }
        else
        {
            other = FORWARD;
        }
    }
    return cross(normalized, other);
}

Vector3f& Vector3f::operator += (const Vector3f& v)
{
    _elements[0] += v._elements[0];
    _elements[1] += v._elements[1];
    _elements[2] += v._elements[2];
    return *this;
}

Vector3f& Vector3f::operator -= (const Vector3f& v)
{
    _elements[0] -= v._elements[0];
    _elements[1] -= v._elements[1];
    _elements[2] -= v._elements[2];
    return *this;
}

Vector3f& Vector3f::operator *= (float f)
{
    _elements[0] *= f;
    _elements[1] *= f;
    _elements[2] *= f;
    return *this;
}

Vector3f& Vector3f::operator /= (float f)
{
    _elements[0] /= f;
    _elements[1] /= f;
    _elements[2] /= f;
    return *this;
}

float Vector3f::dot(const Vector3f& v0, const Vector3f& v1)
{
    return v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2];
}

Vector3f Vector3f::cross(const Vector3f& v0, const Vector3f& v1)
{
    return Vector3f(v0.y() * v1.z() - v0.z() * v1.y(),
                    v0.z() * v1.x() - v0.x() * v1.z(),
                    v0.x() * v1.y() - v0.y() * v1.x());
}

Vector3f Vector3f::lerp(const Vector3f& v0, const Vector3f& v1, float alpha)
{
    return alpha * (v1 - v0) + v0;
}

Vector3f Vector3f::cubicInterpolate(const Vector3f& p0,
                                     const Vector3f& p1,
                                     const Vector3f& p2,
                                     const Vector3f& p3, float t)
{
    // Geometric construction:
    //            t
    //   (t+1)/2     t/2
    // t+1        t	        t-1

    // Bottom level
    Vector3f p0p1 = Vector3f::lerp(p0, p1, t + 1);
    Vector3f p1p2 = Vector3f::lerp(p1, p2, t);
    Vector3f p2p3 = Vector3f::lerp(p2, p3, t - 1);

    // Middle level
    Vector3f p0p1_p1p2 = Vector3f::lerp(p0p1, p1p2, 0.5f * (t + 1));
    Vector3f p1p2_p2p3 = Vector3f::lerp(p1p2, p2p3, 0.5f * t);

    // Top level
    return Vector3f::lerp(p0p1_p1p2, p1p2_p2p3, t);
}

Vector3f operator + (const Vector3f& v0, const Vector3f& v1)
{
    return Vector3f(v0[0] + v1[0], v0[1] + v1[1], v0[2] + v1[2]);
}

Vector3f operator - (const Vector3f& v0, const Vector3f& v1)
{
    return Vector3f(v0[0] - v1[0], v0[1] - v1[1], v0[2] - v1[2]);
}

Vector3f operator * (const Vector3f& v0, const Vector3f& v1)
{
    return Vector3f(v0[0] * v1[0], v0[1] * v1[1], v0[2] * v1[2]);
}

Vector3f operator / (const Vector3f& v0, const Vector3f& v1)
{
    return Vector3f(v0[0] / v1[0], v0[1] / v1[1], v0[2] / v1[2]);
}

Vector3f operator - (const Vector3f& v)
{
    return Vector3f(-v[0], -v[1], -v[2]);
}

Vector3f operator * (float f, const Vector3f& v)
{
    return Vector3f(v[0] * f, v[1] * f, v[2] * f);
}

Vector3f operator * (const Vector3f& v, float f)
{
    return Vector3f(v[0] * f, v[1] * f, v[2] * f);
}

Vector3f operator / (const Vector3f& v, float f)
{
    return Vector3f(v[0] / f, v[1] / f, v[2] / f);
}

bool operator == (const Vector3f& v0, const Vector3f& v1)
{
    return(isEqual(v0.x(), v1.x()) &&
           isEqual(v0.y(), v1.y()) &&
           isEqual(v0.z(), v1.z()));
}

bool operator > (const Vector3f& v0, const Vector3f& v1)
{
    return(v0.x() > v1.x() && v0.y() > v1.y() && v0.z() > v1.z());
}

bool operator < (const Vector3f& v0, const Vector3f& v1)
{
    return(v0.x() < v1.x() && v0.y() < v1.y() && v0.z() < v1.z());
}

bool operator != (const Vector3f& v0, const Vector3f& v1)
{
    return !(v0 == v1);
}

}
