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

#include "Matrix4f.h"
#include "Matrix2f.h"
#include "Matrix3f.h"
#include "Quat4f.h"
#include "Vector3f.h"
#include "Vector4f.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

Matrix4f::Matrix4f(float fill)
{
    for (int i = 0; i < 16; ++i)
    {
        _elements[i] = fill;
    }
}

Matrix4f::Matrix4f(float m00, float m01, float m02, float m03,
                    float m10, float m11, float m12, float m13,
                    float m20, float m21, float m22, float m23,
                    float m30, float m31, float m32, float m33)
{
    _elements[0] = m00;
    _elements[1] = m10;
    _elements[2] = m20;
    _elements[3] = m30;

    _elements[4] = m01;
    _elements[5] = m11;
    _elements[6] = m21;
    _elements[7] = m31;

    _elements[8] = m02;
    _elements[9] = m12;
    _elements[10] = m22;
    _elements[11] = m32;

    _elements[12] = m03;
    _elements[13] = m13;
    _elements[14] = m23;
    _elements[15] = m33;
}

Matrix4f& Matrix4f::operator/=(float d)
{
    for (int ii=0;ii<16;ii++){
        _elements[ii]/=d;
    }
    return *this;
}

Matrix4f::Matrix4f(const Vector4f& v0,
                    const Vector4f& v1,
                    const Vector4f& v2,
                    const Vector4f& v3,
                    bool setColumns)
{
    if (setColumns)
    {
        setCol(0, v0);
        setCol(1, v1);
        setCol(2, v2);
        setCol(3, v3);
    }
    else
    {
        setRow(0, v0);
        setRow(1, v1);
        setRow(2, v2);
        setRow(3, v3);
    }
}

Matrix4f::Matrix4f(const Matrix4f& rm)
{
    memcpy(_elements, rm._elements, sizeof(_elements));
}

Matrix4f& Matrix4f::operator = (const Matrix4f& rm)
{
    if (this != &rm)
    {
        memcpy(_elements, rm._elements, sizeof(_elements));
    }
    return *this;
}

const float& Matrix4f::operator () (int i, int j) const
{
    return _elements[j * 4 + i];
}

float& Matrix4f::operator () (int i, int j)
{
    return _elements[j * 4 + i];
}

Vector4f Matrix4f::getRow(int i) const
{
    return Vector4f(_elements[i],
                    _elements[i + 4],
                    _elements[i + 8],
                    _elements[i + 12]);
}

void Matrix4f::setRow(int i, const Vector4f& v)
{
    _elements[i] = v.x();
    _elements[i + 4] = v.y();
    _elements[i + 8] = v.z();
    _elements[i + 12] = v.w();
}

Vector4f Matrix4f::getCol(int j) const
{
    int colStart = 4 * j;

    return Vector4f(_elements[colStart],
                    _elements[colStart + 1],
                    _elements[colStart + 2],
                    _elements[colStart + 3]);
}

void Matrix4f::setCol(int j, const Vector4f& v)
{
    int colStart = 4 * j;

    _elements[colStart] = v.x();
    _elements[colStart + 1] = v.y();
    _elements[colStart + 2] = v.z();
    _elements[colStart + 3] = v.w();
}

Matrix2f Matrix4f::getSubmatrix2x2(int i0, int j0) const
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

Matrix3f Matrix4f::getSubmatrix3x3(int i0, int j0) const
{
    Matrix3f out;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            out(i, j) = (*this)(i + i0, j + j0);
        }
    }

    return out;
}

void Matrix4f::setSubmatrix2x2(int i0, int j0, const Matrix2f& m)
{
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 2; ++j)
        {
            (*this)(i + i0, j + j0) = m(i, j);
        }
    }
}

void Matrix4f::setSubmatrix3x3(int i0, int j0, const Matrix3f& m)
{
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            (*this)(i + i0, j + j0) = m(i, j);
        }
    }
}

float Matrix4f::determinant() const
{
    float m00 = _elements[0];
    float m10 = _elements[1];
    float m20 = _elements[2];
    float m30 = _elements[3];

    float m01 = _elements[4];
    float m11 = _elements[5];
    float m21 = _elements[6];
    float m31 = _elements[7];

    float m02 = _elements[8];
    float m12 = _elements[9];
    float m22 = _elements[10];
    float m32 = _elements[11];

    float m03 = _elements[12];
    float m13 = _elements[13];
    float m23 = _elements[14];
    float m33 = _elements[15];

    float cofactor00 =  Matrix3f::determinant3x3(m11, m12, m13,
                                                  m21, m22, m23,
                                                  m31, m32, m33);
    float cofactor01 = -Matrix3f::determinant3x3(m12, m13, m10,
                                                  m22, m23, m20,
                                                  m32, m33, m30);
    float cofactor02 =  Matrix3f::determinant3x3(m13, m10, m11,
                                                  m23, m20, m21,
                                                  m33, m30, m31);
    float cofactor03 = -Matrix3f::determinant3x3(m10, m11, m12,
                                                  m20, m21, m22,
                                                  m30, m31, m32);

    return(m00 * cofactor00 + m01 * cofactor01 + m02 * cofactor02 + m03 * cofactor03);
}

Matrix4f Matrix4f::inverse(bool* pbIsSingular, float epsilon) const
{
    float m00 = _elements[0];
    float m10 = _elements[1];
    float m20 = _elements[2];
    float m30 = _elements[3];

    float m01 = _elements[4];
    float m11 = _elements[5];
    float m21 = _elements[6];
    float m31 = _elements[7];

    float m02 = _elements[8];
    float m12 = _elements[9];
    float m22 = _elements[10];
    float m32 = _elements[11];

    float m03 = _elements[12];
    float m13 = _elements[13];
    float m23 = _elements[14];
    float m33 = _elements[15];

    float cofactor00 =  Matrix3f::determinant3x3(m11, m12, m13,
                                                 m21, m22, m23,
                                                 m31, m32, m33);
    float cofactor01 = -Matrix3f::determinant3x3(m12, m13, m10,
                                                 m22, m23, m20,
                                                 m32, m33, m30);
    float cofactor02 =  Matrix3f::determinant3x3(m13, m10, m11,
                                                 m23, m20, m21,
                                                 m33, m30, m31);
    float cofactor03 = -Matrix3f::determinant3x3(m10, m11, m12,
                                                 m20, m21, m22,
                                                 m30, m31, m32);
    
    float cofactor10 = -Matrix3f::determinant3x3(m21, m22, m23,
                                                 m31, m32, m33,
                                                 m01, m02, m03);
    float cofactor11 =  Matrix3f::determinant3x3(m22, m23, m20,
                                                 m32, m33, m30,
                                                 m02, m03, m00);
    float cofactor12 = -Matrix3f::determinant3x3(m23, m20, m21,
                                                 m33, m30, m31,
                                                 m03, m00, m01);
    float cofactor13 =  Matrix3f::determinant3x3(m20, m21, m22,
                                                 m30, m31, m32,
                                                 m00, m01, m02);
    
    float cofactor20 =  Matrix3f::determinant3x3(m31, m32, m33,
                                                 m01, m02, m03,
                                                 m11, m12, m13);
    float cofactor21 = -Matrix3f::determinant3x3(m32, m33, m30,
                                                 m02, m03, m00,
                                                 m12, m13, m10);
    float cofactor22 =  Matrix3f::determinant3x3(m33, m30, m31,
                                                 m03, m00, m01,
                                                 m13, m10, m11);
    float cofactor23 = -Matrix3f::determinant3x3(m30, m31, m32,
                                                 m00, m01, m02,
                                                 m10, m11, m12);
    
    float cofactor30 = -Matrix3f::determinant3x3(m01, m02, m03,
                                                 m11, m12, m13,
                                                 m21, m22, m23);
    float cofactor31 =  Matrix3f::determinant3x3(m02, m03, m00,
                                                 m12, m13, m10,
                                                 m22, m23, m20);
    float cofactor32 = -Matrix3f::determinant3x3(m03, m00, m01,
                                                 m13, m10, m11,
                                                 m23, m20, m21);
    float cofactor33 =  Matrix3f::determinant3x3(m00, m01, m02,
                                                 m10, m11, m12,
                                                 m20, m21, m22);

    float determinant = m00 * cofactor00 + m01 * cofactor01 +
            m02 * cofactor02 + m03 * cofactor03;

    bool isSingular = (std::fabs(determinant) < epsilon);
    if (isSingular)
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = true;
        }
        return Matrix4f();
    }
    else
    {
        if (pbIsSingular != nullptr)
        {
            *pbIsSingular = false;
        }

        float reciprocalDeterminant = 1.0f / determinant;

        return Matrix4f(cofactor00 * reciprocalDeterminant,
                        cofactor10 * reciprocalDeterminant,
                        cofactor20 * reciprocalDeterminant,
                        cofactor30 * reciprocalDeterminant,

                        cofactor01 * reciprocalDeterminant,
                        cofactor11 * reciprocalDeterminant,
                        cofactor21 * reciprocalDeterminant,
                        cofactor31 * reciprocalDeterminant,

                        cofactor02 * reciprocalDeterminant,
                        cofactor12 * reciprocalDeterminant,
                        cofactor22 * reciprocalDeterminant,
                        cofactor32 * reciprocalDeterminant,

                        cofactor03 * reciprocalDeterminant,
                        cofactor13 * reciprocalDeterminant,
                        cofactor23 * reciprocalDeterminant,
                        cofactor33 * reciprocalDeterminant);
    }
}

void Matrix4f::transpose()
{
    float temp;

    for (int i = 0; i < 3; ++i)
    {
        for (int j = i + 1; j < 4; ++j)
        {
            temp = (*this)(i, j);
            (* this)(i, j) = (*this)(j, i);
            (*this)(j, i) = temp;
        }
    }
}

Matrix4f Matrix4f::transposed() const
{
    Matrix4f out;
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            out(j, i) = (*this)(i, j);
        }
    }

    return out;
}

Matrix4f::operator float* ()
{
    return _elements;
}

Matrix4f::operator const float* ()const
{
    return _elements;
}

void Matrix4f::print()
{
    printf("[%.4f %.4f %.4f %.4f]\n"
           "[%.4f %.4f %.4f %.4f]\n"
           "[%.4f %.4f %.4f %.4f]\n"
           "[%.4f %.4f %.4f %.4f]\n",
    F2D(_elements[0]), F2D(_elements[4]), F2D(_elements[8]), F2D(_elements[12]),
    F2D(_elements[1]), F2D(_elements[5]), F2D(_elements[9]), F2D(_elements[13]),
    F2D(_elements[2]), F2D(_elements[6]), F2D(_elements[10]), F2D(_elements[14]),
    F2D(_elements[3]), F2D(_elements[7]), F2D(_elements[11]), F2D(_elements[15]));
}

Matrix4f Matrix4f::ones()
{
    Matrix4f m;
    for (int i = 0; i < 16; ++i)
    {
        m._elements[i] = 1;
    }

    return m;
}

Matrix4f Matrix4f::identity()
{
    Matrix4f m;

    m(0, 0) = 1;
    m(1, 1) = 1;
    m(2, 2) = 1;
    m(3, 3) = 1;

    return m;
}

Matrix4f Matrix4f::translation(float x, float y, float z)
{
    return Matrix4f(1, 0, 0, x,
                    0, 1, 0, y,
                    0, 0, 1, z,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::translation(const Vector3f& rTranslation)
{
    return Matrix4f(1, 0, 0, rTranslation.x(),
                    0, 1, 0, rTranslation.y(),
                    0, 0, 1, rTranslation.z(),
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::rotateX(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix4f(1, 0, 0, 0,
                    0, c, -s, 0,
                    0, s, c, 0,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::rotateY(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix4f(c, 0, s, 0,
                    0, 1, 0, 0,
                    -s, 0, c, 0,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::rotateZ(float radians)
{
    float c = cos(radians);
    float s = sin(radians);

    return Matrix4f(c, -s, 0, 0,
                    s, c, 0, 0,
                    0, 0, 1, 0,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::rotation(const Vector3f& rDirection, float radians)
{
    Vector3f normalizedDirection = rDirection.normalized();

    float cosTheta = cos(radians);
    float sinTheta = sin(radians);

    float x = normalizedDirection.x();
    float y = normalizedDirection.y();
    float z = normalizedDirection.z();

    return Matrix4f(x * x * (1.0f - cosTheta) + cosTheta,
                    y * x * (1.0f - cosTheta) - z * sinTheta,
                    z * x * (1.0f - cosTheta) + y * sinTheta,
                    0.0f,

                    x * y * (1.0f - cosTheta) + z * sinTheta,
                    y * y * (1.0f - cosTheta) + cosTheta,
                    z * y * (1.0f - cosTheta) - x * sinTheta,
                    0.0f,

                    x * z * (1.0f - cosTheta) - y * sinTheta,
                    y * z * (1.0f - cosTheta) + x * sinTheta,
                    z * z * (1.0f - cosTheta) + cosTheta,
                    0.0f,

                    0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix4f Matrix4f::rotation(const Quat4f& q)
{
    Quat4f qq = q.normalized();

    float xx = qq.x() * qq.x();
    float yy = qq.y() * qq.y();
    float zz = qq.z() * qq.z();

    float xy = qq.x() * qq.y();
    float zw = qq.z() * qq.w();

    float xz = qq.x() * qq.z();
    float yw = qq.y() * qq.w();

    float yz = qq.y() * qq.z();
    float xw = qq.x() * qq.w();

    return Matrix4f(1.0f - 2.0f * (yy + zz),
                    2.0f * (xy - zw),
                    2.0f * (xz + yw),
                    0.0f,

                    2.0f * (xy + zw),
                    1.0f - 2.0f * (xx + zz),
                    2.0f * (yz - xw),
                    0.0f,

                    2.0f * (xz - yw),
                    2.0f * (yz + xw),
                    1.0f - 2.0f * (xx + yy),
                    0.0f,

                    0.0f, 0.0f, 0.0f, 1.0f);
}

Matrix4f Matrix4f::scaling(float sx, float sy, float sz)
{
    return Matrix4f(sx, 0, 0, 0,
                    0, sy, 0, 0,
                    0, 0, sz, 0,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::uniformScaling(float s)
{
    return Matrix4f(s, 0, 0, 0,
                    0, s, 0, 0,
                    0, 0, s, 0,
                    0, 0, 0, 1);
}

Matrix4f Matrix4f::randomRotation(float u0, float u1, float u2)
{
    return Matrix4f::rotation(Quat4f::randomRotation(u0, u1, u2));
}

Matrix4f Matrix4f::lookAt(const Vector3f& eye, const Vector3f& center, const Vector3f& up)
{
    // z is negative forward
    Vector3f z = (eye - center).normalized();
    Vector3f y = up;
    Vector3f x = Vector3f::cross(y, z);

    // the x, y, and z vectors define the orthonormal coordinate system
    // the affine part defines the overall translation
    Matrix4f view;

    view.setRow(0, Vector4f(x, -Vector3f::dot(x, eye)));
    view.setRow(1, Vector4f(y, -Vector3f::dot(y, eye)));
    view.setRow(2, Vector4f(z, -Vector3f::dot(z, eye)));
    view.setRow(3, Vector4f(0, 0, 0, 1));

    return view;
}

Matrix4f Matrix4f::orthographicProjection(float width, float height,
                                          float zNear, float zFar, bool directX)
{
    Matrix4f m;

    m(0, 0) = 2.0f / width;
    m(1, 1) = 2.0f / height;
    m(3, 3) = 1.0f;

    m(0, 3) = -1;
    m(1, 3) = -1;

    if (directX)
    {
        m(2, 2) = 1.0f / (zNear - zFar);
        m(2, 3) = zNear / (zNear - zFar);
    }
    else
    {
        m(2, 2) = 2.0f / (zNear - zFar);
        m(2, 3) = (zNear + zFar) / (zNear - zFar);
    }

    return m;
}

Matrix4f Matrix4f::orthographicProjection(float left, float right,
                                          float bottom, float top,
                                          float zNear, float zFar, bool directX)
{
    Matrix4f m;

    m(0, 0) = 2.0f / (right - left);
    m(1, 1) = 2.0f / (top - bottom);
    m(3, 3) = 1.0f;

    m(0, 3) = (left + right) / (left - right);
    m(1, 3) = (top + bottom) / (bottom - top);

    if (directX)
    {
        m(2, 2) = 1.0f / (zNear - zFar);
        m(2, 3) = zNear / (zNear - zFar);
    }
    else
    {
        m(2, 2) = 2.0f / (zNear - zFar);
        m(2, 3) = (zNear + zFar) / (zNear - zFar);
    }

    return m;
}

Matrix4f Matrix4f::perspectiveProjection(float fLeft, float fRight,
                                          float fBottom, float fTop,
                                          float fZNear, float fZFar,
                                          bool directX)
{
    Matrix4f projection; // Zero matrix

    projection(0, 0) = (2.0f * fZNear) / (fRight - fLeft);
    projection(1, 1) = (2.0f * fZNear) / (fTop - fBottom);
    projection(0, 2) = (fRight + fLeft) / (fRight - fLeft);
    projection(1, 2) = (fTop + fBottom) / (fTop - fBottom);
    projection(3, 2) = -1;

    if (directX)
    {
        projection(2, 2) = fZFar / (fZNear - fZFar);
        projection(2, 3) = (fZNear * fZFar) / (fZNear - fZFar);
    }
    else
    {
        projection(2, 2) = (fZNear + fZFar) / (fZNear - fZFar);
        projection(2, 3) = (2.0f * fZNear * fZFar) / (fZNear - fZFar);
    }

    return projection;
}

Matrix4f Matrix4f::perspectiveProjection(float fovYRadians, float aspect,
                                         float zNear, float zFar, bool directX)
{
    Matrix4f m; // Zero matrix

    float yScale = 1.f / tanf(0.5f * fovYRadians);
    float xScale = yScale / aspect;

    m(0, 0) = xScale;
    m(1, 1) = yScale;
    m(3, 2) = -1;

    if (directX)
    {
        m(2, 2) = zFar / (zNear - zFar);
        m(2, 3) = zNear * zFar / (zNear - zFar);
    }
    else
    {
        m(2, 2) = (zFar + zNear) / (zNear - zFar);
        m(2, 3) = 2.f * zFar * zNear / (zNear - zFar);
    }

    return m;
}

Matrix4f Matrix4f::infinitePerspectiveProjection(float fLeft, float fRight,
                                                 float fBottom, float fTop,
                                                 float fZNear, bool directX)
{
    Matrix4f projection;

    projection(0, 0) = (2.0f * fZNear) / (fRight - fLeft);
    projection(1, 1) = (2.0f * fZNear) / (fTop - fBottom);
    projection(0, 2) = (fRight + fLeft) / (fRight - fLeft);
    projection(1, 2) = (fTop + fBottom) / (fTop - fBottom);
    projection(3, 2) = -1;

    if (directX)
    {
        projection(2, 2) = -1.0f;
        projection(2, 3) = -fZNear;
    }
    else
    {
        projection(2, 2) = -1.0f;
        projection(2, 3) = -2.0f * fZNear;
    }

    return projection;
}

Vector4f operator * (const Matrix4f& m, const Vector4f& v)
{
    Vector4f output(0, 0, 0, 0);

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            output[i] += m(i, j) * v[j];
        }
    }

    return output;
}

Matrix4f operator * (const Matrix4f& x, const Matrix4f& y)
{
    Matrix4f product; // Zeroes

    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            for (int k = 0; k < 4; ++k)
            {
                product(i, k) += x(i, j) * y(j, k);
            }
        }
    }

    return product;
}

}
