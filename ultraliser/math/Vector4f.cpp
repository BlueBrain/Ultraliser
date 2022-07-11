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

#include "Vector4f.h"
#include "Vector2f.h"
#include "Vector3f.h"
#include <utilities/TypeConversion.h>
#include <math/Functions.h>

namespace Ultraliser
{

Vector4f::Vector4f(float f)
{
    _elements[0] = f;
    _elements[1] = f;
    _elements[2] = f;
    _elements[3] = f;
}

Vector4f::Vector4f(float fx, float fy, float fz, float fw)
{
    _elements[0] = fx;
    _elements[1] = fy;
    _elements[2] = fz;
    _elements[3] = fw;
}

Vector4f::Vector4f(float buffer[4])
{
    _elements[0] = buffer[0];
    _elements[1] = buffer[1];
    _elements[2] = buffer[2];
    _elements[3] = buffer[3];
}

Vector4f::Vector4f(const Vector2f& xy, float z, float w)
{
    _elements[0] = xy.x();
    _elements[1] = xy.y();
    _elements[2] = z;
    _elements[3] = w;
}

Vector4f::Vector4f(float x, const Vector2f& yz, float w)
{
    _elements[0] = x;
    _elements[1] = yz.x();
    _elements[2] = yz.y();
    _elements[3] = w;
}

Vector4f::Vector4f(float x, float y, const Vector2f& zw)
{
    _elements[0] = x;
    _elements[1] = y;
    _elements[2] = zw.x();
    _elements[3] = zw.y();
}

Vector4f::Vector4f(const Vector2f& xy, const Vector2f& zw)
{
    _elements[0] = xy.x();
    _elements[1] = xy.y();
    _elements[2] = zw.x();
    _elements[3] = zw.y();
}

Vector4f::Vector4f(const Vector3f& xyz, float w)
{
    _elements[0] = xyz.x();
    _elements[1] = xyz.y();
    _elements[2] = xyz.z();
    _elements[3] = w;
}

Vector4f::Vector4f(float x, const Vector3f& yzw)
{
    _elements[0] = x;
    _elements[1] = yzw.x();
    _elements[2] = yzw.y();
    _elements[3] = yzw.z();
}

Vector4f::Vector4f(const Vector4f& rv)
{
    _elements[0] = rv._elements[0];
    _elements[1] = rv._elements[1];
    _elements[2] = rv._elements[2];
    _elements[3] = rv._elements[3];
}

Vector4f& Vector4f::operator = (const Vector4f& rv)
{
    if (this != &rv)
    {
        _elements[0] = rv._elements[0];
        _elements[1] = rv._elements[1];
        _elements[2] = rv._elements[2];
        _elements[3] = rv._elements[3];
    }
    return *this;
}

const float& Vector4f::operator [] (int i) const
{
    return _elements[i];
}

float& Vector4f::operator [] (int i)
{
    return _elements[i];
}

float& Vector4f::x()
{
    return _elements[0];
}

float& Vector4f::y()
{
    return _elements[1];
}

float& Vector4f::z()
{
    return _elements[2];
}

float& Vector4f::w()
{
    return _elements[3];
}

float Vector4f::x() const
{
    return _elements[0];
}

float Vector4f::y() const
{
    return _elements[1];
}

float Vector4f::z() const
{
    return _elements[2];
}

float Vector4f::w() const
{
    return _elements[3];
}

Vector2f Vector4f::xy() const
{
    return Vector2f(_elements[0], _elements[1]);
}

Vector2f Vector4f::yz() const
{
    return Vector2f(_elements[1], _elements[2]);
}

Vector2f Vector4f::zw() const
{
    return Vector2f(_elements[2], _elements[3]);
}

Vector2f Vector4f::wx() const
{
    return Vector2f(_elements[3], _elements[0]);
}

Vector3f Vector4f::xyz() const
{
    return Vector3f(_elements[0], _elements[1], _elements[2]);
}

Vector3f Vector4f::yzw() const
{
    return Vector3f(_elements[1], _elements[2], _elements[3]);
}

Vector3f Vector4f::zwx() const
{
    return Vector3f(_elements[2], _elements[3], _elements[0]);
}

Vector3f Vector4f::wxy() const
{
    return Vector3f(_elements[3], _elements[0], _elements[1]);
}

Vector3f Vector4f::xyw() const
{
    return Vector3f(_elements[0], _elements[1], _elements[3]);
}

Vector3f Vector4f::yzx() const
{
    return Vector3f(_elements[1], _elements[2], _elements[0]);
}

Vector3f Vector4f::zwy() const
{
    return Vector3f(_elements[2], _elements[3], _elements[1]);
}

Vector3f Vector4f::wxz() const
{
    return Vector3f(_elements[3], _elements[0], _elements[2]);
}

float Vector4f::abs() const
{
    const double value = sqrt(F2D(_elements[0] * _elements[0]) +
                              F2D(_elements[1] * _elements[1]) +
                              F2D(_elements[2] * _elements[2]) +
                              F2D(_elements[3] * _elements[3]));
    return D2F(value);
}

float Vector4f::absSquared() const
{
    const float value = _elements[0] * _elements[0] +
                        _elements[1] * _elements[1] +
                        _elements[2] * _elements[2] +
                        _elements[3] * _elements[3];
    return value;
}

void Vector4f::normalize()
{
    float norm = sqrt(_elements[0] * _elements[0] +
                      _elements[1] * _elements[1] +
                      _elements[2] * _elements[2] +
                      _elements[3] * _elements[3]);

    _elements[0] = _elements[0] / norm;
    _elements[1] = _elements[1] / norm;
    _elements[2] = _elements[2] / norm;
    _elements[3] = _elements[3] / norm;
}

Vector4f Vector4f::normalized() const
{
    float length = abs();
    return Vector4f (_elements[0] / length, _elements[1] / length,
            _elements[2] / length, _elements[3] / length);
}

void Vector4f::homogenize()
{
    if (!isZero(_elements[3]))
    {
        _elements[0] /= _elements[3];
        _elements[1] /= _elements[3];
        _elements[2] /= _elements[3];
        _elements[3] = 1;
    }
}

Vector4f Vector4f::homogenized() const
{
    if (!isZero(_elements[3]))
    {
        return Vector4f (_elements[0] / _elements[3],
                         _elements[1] / _elements[3],
                         _elements[2] / _elements[3], 1);
    }
    else
    {
        return Vector4f (_elements[0], _elements[1], _elements[2], _elements[3]);
    }
}

void Vector4f::negate()
{
    _elements[0] = -_elements[0];
    _elements[1] = -_elements[1];
    _elements[2] = -_elements[2];
    _elements[3] = -_elements[3];
}

Vector4f::operator const float* () const
{
    return _elements;
}

Vector4f::operator float* ()
{
    return _elements;
}

void Vector4f::print() const
{
    printf("< %.4f, %.4f, %.4f, %.4f >\n", F2D(_elements[0]),
                                           F2D(_elements[1]),
                                           F2D(_elements[2]),
                                           F2D(_elements[3]));
}

float Vector4f::dot(const Vector4f& v0, const Vector4f& v1)
{
    return v0.x() * v1.x() +
           v0.y() * v1.y() +
           v0.z() * v1.z() +
           v0.w() * v1.w();
}

Vector4f Vector4f::lerp(const Vector4f& v0, const Vector4f& v1, float alpha)
{
    return alpha * (v1 - v0) + v0;
}

Vector4f operator + (const Vector4f& v0, const Vector4f& v1)
{
    return Vector4f(v0.x() + v1.x(),
                     v0.y() + v1.y(),
                     v0.z() + v1.z(),
                     v0.w() + v1.w());
}

Vector4f operator - (const Vector4f& v0, const Vector4f& v1)
{
    return Vector4f(v0.x() - v1.x(),
                     v0.y() - v1.y(),
                     v0.z() - v1.z(),
                     v0.w() - v1.w());
}

Vector4f operator * (const Vector4f& v0, const Vector4f& v1)
{
    return Vector4f(v0.x() * v1.x(),
                     v0.y() * v1.y(),
                     v0.z() * v1.z(),
                     v0.w() * v1.w());
}

Vector4f operator / (const Vector4f& v0, const Vector4f& v1)
{
    return Vector4f(v0.x() / v1.x(),
                     v0.y() / v1.y(),
                     v0.z() / v1.z(),
                     v0.w() / v1.w());
}

Vector4f operator - (const Vector4f& v)
{
    return Vector4f(-v.x(), -v.y(), -v.z(), -v.w());
}

Vector4f operator * (float f, const Vector4f& v)
{
    return Vector4f(f * v.x(), f * v.y(), f * v.z(), f * v.w());
}

Vector4f operator * (const Vector4f& v, float f)
{
    return Vector4f(f * v.x(), f * v.y(), f * v.z(), f * v.w());
}

Vector4f operator / (const Vector4f& v, float f)
{
    return Vector4f(v[0] / f, v[1] / f, v[2] / f, v[3] / f);
}

bool operator == (const Vector4f& v0, const Vector4f& v1)
{
    return(isEqual(v0.x(), v1.x()) &&
           isEqual(v0.y(), v1.y()) &&
           isEqual(v0.z(), v1.z()) &&
           isEqual(v0.w(), v1.w()));
}

bool operator != (const Vector4f& v0, const Vector4f& v1)
{
    return !(v0 == v1);
}

}
