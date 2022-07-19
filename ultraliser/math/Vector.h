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

#pragma once

#include <math/Functions.h>
#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The Vec struct
 * Defines a thin wrapper around fixed size C-style arrays, using template
 * parameters, which is useful for dealing with vectors of different dimensions.
 * For example, float[3] is equivalent to Vec< 3, float >.
 *
 * Entries in the vector are accessed with the overloaded [] operator, so for
 * example if x is a Vec< 3, float >, then the middle entry is x[1].
 *
 * @note
 * Arithmetic operators are appropriately overloaded, and functions are defined
 * for additional operations (such as dot-products, norms, cross-products, etc.)
 */
template< size_t N, class T >
struct Vec
{
public:

    /**
     * @brief v
     */
    T v[N];

public:

    /**
     * @brief Vec<N, T>
     */
    Vec< N, T >(void)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief Vec<N, T>
     * @param value
     */
    explicit Vec< N, T >(T value)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] = value;
        }
    }

    /**
     * @brief Vec<N, T>
     * @param source
     */
    template< class S>
    explicit Vec< N, T >(const S* source)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] = static_cast<T>(source[i]);
        }
    }

    /**
     * @brief Vec<N, T>
     * @param source
     */
    template < class S >
    explicit Vec< N, T >(const Vec< N,S>& source)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] = static_cast<T>(source[i]);
        }
    }

    /**
     * @brief Vec<N, T>
     * @param v0
     * @param v1
     */
    Vec< N, T >(T v0, T v1)
    {
        assert(N == 2);
        v[0] = v0;
        v[1] = v1;
    }

    /**
     * @brief Vec<N, T>
     * @param v0
     * @param v1
     * @param v2
     */
    Vec< N, T >(T v0, T v1, T v2)
    {
        assert(N == 3);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
    }

    /**
     * @brief Vec<N, T>
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     */
    Vec< N, T >(T v0, T v1, T v2, T v3)
    {
        assert(N == 4);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
    }

    /**
     * @brief Vec<N, T>
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     */
    Vec< N, T >(T v0, T v1, T v2, T v3, T v4)
    {
        assert(N == 5);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
    }

    /**
     * @brief Vec<N, T>
     * @param v0
     * @param v1
     * @param v2
     * @param v3
     * @param v4
     * @param v5
     */
    Vec< N, T >(T v0, T v1, T v2, T v3, T v4, T v5)
    {
        assert(N == 6);
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
        v[3] = v3;
        v[4] = v4;
        v[5] = v5;
    }

    /**
     * @brief operator []
     * @param index
     * @return
     */
    T& operator[](int32_t index)
    {
        assert(0 <= index && (uint32_t) index < N);
        return v[index];
    }

    /**
     * @brief operator []
     * @param index
     * @return
     */
    const T& operator[](int32_t index) const
    {
        assert(0 <= index && index < N);
        return v[index];
    }

    /**
     * @brief nonZero
     * @return
     */
    bool nonZero(void) const
    {
        for (size_t i = 0; i < N; ++i)
        {
            if(v[i])
            {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief x
     * @return
     */
    T x(void) const
    {
        return v[0];
    }

    /**
     * @brief y
     * @return
     */
    T y(void) const
    {
        return v[1];
    }

    /**
     * @brief z
     * @return
     */
    T z(void) const
    {
        return v[2];
    }

    /**
     * @brief w
     * @return
     */
    T w(void) const
    {
        return v[3];
    }

    /**
     * @brief operator +=
     * @param w
     * @return
     */
    Vec< N, T > operator+= (const Vec< N, T >& w)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] +=w [i];
        }

        return *this;
    }

    /**
     * @brief operator +
     * @param w
     * @return
     */
    Vec< N, T > operator+(const Vec< N, T >& w) const
    {
        Vec< N, T > sum(*this);
        sum += w;
        return sum;
    }

    /**
     * @brief operator -=
     * @param w
     * @return
     */
    Vec< N, T > operator-= (const Vec< N, T >& w)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] -= w[i];
        }
        return *this;
    }

    /**
     * @brief operator -
     * @return
     */
    Vec< N, T > operator-(void) const
    {
        Vec< N, T > negative;
        for (size_t i = 0; i < N; ++i)
        {
            negative.v[i] = -v[i];
        }
        return negative;
    }

    /**
     * @brief operator -
     * @param w
     * @return
     */
    Vec< N, T > operator-(const Vec< N, T >& w) const
    {
        Vec< N, T > diff(*this);
        diff -= w;
        return diff;
    }

    /**
     * @brief operator *=
     * @param a
     * @return
     */
    Vec< N, T > operator*= (T a)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] *= a;
        }
        return *this;
    }

    /**
     * @brief operator *
     * @param a
     * @return
     */
    Vec< N, T > operator*(T a) const
    {
        Vec< N, T > w(*this);
        w *= a;
        return w;
    }

    /**
     * @brief operator *=
     * @param w
     * @return
     */
    Vec< N, T > operator*= (const Vec< N, T >& w)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] *= w.v[i];
        }
        return *this;
    }

    /**
     * @brief operator *
     * @param w
     * @return
     */
    Vec< N, T > operator*(const Vec< N, T >& w) const
    {
        Vec< N, T > componentWiseProduct;
        for (size_t i = 0; i < N; ++i)
        {
            componentWiseProduct[i] = v[i] * w.v[i];
        }
        return componentWiseProduct;
    }

    /**
     * @brief operator /=
     * @param a
     * @return
     */
    Vec< N, T > operator/= (T a)
    {
        for (size_t i = 0; i < N; ++i)
        {
            v[i] /= a;
        }
        return *this;
    }

    /**
     * @brief operator /
     * @param a
     * @return
     */
    Vec< N, T > operator/(T a) const
    {
        Vec< N, T > w(*this);
        w /= a;
        return w;
    }
};

// Vectors
typedef Vec< 2, int64_t >       Vec2i_64;
typedef Vec< 2, uint64_t >      Vec2ui_64;
typedef Vec< 2, double >        Vec2f_64;
typedef Vec< 2, float >         Vec2f_32;
typedef Vec< 2, int32_t >       Vec2i_32;
typedef Vec< 2, uint32_t >      Vec2ui_32;
typedef Vec< 2, int16_t >       Vec2i_16;
typedef Vec< 2, uint16_t >      Vec2ui_16;
typedef Vec< 2, int8_t >        Vec2i_8;
typedef Vec< 2, uint8_t >       Vec2ui_8;

typedef Vec< 3, int64_t >       Vec3i_64;
typedef Vec< 3, uint64_t >      Vec3ui_64;
typedef Vec< 3, double >        Vec3f_64;
typedef Vec< 3, float >         Vec3f_32;
typedef Vec< 3, int32_t >       Vec3i_32;
typedef Vec< 3, uint32_t >      Vec3ui_32;
typedef Vec< 3, int16_t >       Vec3i_16;
typedef Vec< 3, uint16_t >      Vec3ui_16;
typedef Vec< 3, int8_t >        Vec3i_8;
typedef Vec< 3, uint8_t >       Vec3ui_8;

typedef Vec< 4, int64_t >       Vec4i_64;
typedef Vec< 4, uint64_t >      Vec4ui_64;
typedef Vec< 4, double >        Vec4f_64;
typedef Vec< 4, float >         Vec4f_32;
typedef Vec< 4, int32_t >       Vec4i_32;
typedef Vec< 4, uint32_t >      Vec4ui_32;
typedef Vec< 4, int16_t >       Vec4i_16;
typedef Vec< 4, uint16_t >      Vec4ui_16;
typedef Vec< 4, int8_t >        Vec4i_8;
typedef Vec< 4, uint8_t >       Vec4ui_8;

typedef Vec< 6, int64_t >       Vec6i_64;
typedef Vec< 6, uint64_t >      Vec6ui_64;
typedef Vec< 6, double >        Vec6f_64;
typedef Vec< 6, float  >        Vec6f_32;
typedef Vec< 6, uint32_t >      Vec6ui_32;
typedef Vec< 6, int32_t >       Vec6i_32;
typedef Vec< 6, int16_t >       Vec6i_16;
typedef Vec< 6, uint16_t >      Vec6ui_16;
typedef Vec< 6, int8_t >        Vec6i_8;
typedef Vec< 6, uint8_t >       Vec6ui_8;


/**
 * @brief mag2
 * @param a
 * @return
 */
template< uint64_t N, class T >
T mag2(const Vec< N, T >& a)
{
    T l = sqr(a.v[0]);
    for (size_t i = 1; i < N; ++i)
    {
        l += sqr(a.v[i]);
    }
    return l;
}

/**
 * @brief mag
 * @param a
 * @return
 */
template< uint64_t N, class T >
T mag(const Vec< N, T >& a)
{
    return sqrt(mag2(a));
}

/**
 * @brief dist2
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline T dist2(const Vec< N, T >& a, const Vec< N, T >& b)
{ 
    T d = sqr(a.v[0] - b.v[0]);
    for (size_t i = 1; i < N; ++i)
    {
        d += sqr(a.v[i] - b.v[i]);
    }
    return d;
}

/**
 * @brief dist
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline T dist(const Vec< N, T >& a, const Vec< N, T >& b)
{
    return std::sqrt(dist2(a, b));
}

/**
 * @brief normalize
 * @param a
 */
template< uint64_t N, class T >
inline void normalize(Vec< N, T >& a)
{
    a /= mag(a);
}

/**
 * @brief normalized
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, T > normalized(const Vec< N, T >& a)
{
    return a / mag(a);
}

/**
 * @brief infnorm
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline T infnorm(const Vec< N, T >& a)
{
    T d = std::fabs(a.v[0]);
    for (size_t i = 1; i < N; ++i)
    {
        d = max(std::fabs(a.v[i]), d);
    }
    return d;
}

/**
 * @brief zero
 * @param a
 */
template< uint64_t N, class T >
void zero(Vec< N, T >& a)
{ 
    for (size_t i = 0; i < N; ++i)
    {
        a.v[i] = 0;
    }
}

/**
 * @brief operator <<
 * @param out
 * @param v
 * @return
 */
template< uint64_t N, class T >
std::ostream &operator<<(std::ostream& out, const Vec< N, T >& v)
{
    out << v.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        out << ' ' << v.v[i];
    }
    return out;
}

/**
 * @brief operator >>
 * @param in
 * @param v
 * @return
 */
template< uint64_t N, class T >
std::istream &operator>>(std::istream& in, Vec< N, T >& v)
{
    in >> v.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        in >> v.v[i];
    }
    return in;
}

/**
 * @brief operator ==
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline bool operator == (const Vec< N, T >& a, const Vec< N, T >& b)
{ 
    bool t = (a.v[0] ==  b.v[0]);
    uint32_t i = 1;
    while (i < N && t)
    {
        t = t && (a.v[i] == b.v[i]);
        ++i;
    }
    return t;
}

/**
 * @brief operator !=
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline bool operator!= (const Vec< N, T >& a, const Vec< N, T >& b)
{ 
    bool t = (a.v[0] != b.v[0]);
    uint32_t i = 1;
    while (i < N && !t)
    {
        t = t || (a.v[i] != b.v[i]);
        ++i;
    }
    return t;
}

/**
 * @brief operator *
 * @param a
 * @param v
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, T > operator*(T a, const Vec< N, T >& v)
{
    Vec< N, T > w(v);
    w *= a;
    return w;
}

/**
 * @brief min
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline T min(const Vec< N, T >& a)
{
    T m = a.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        if(a.v[i] < m)
        {
            m = a.v[i];
        }
    }
    return m;
}

/**
 * @brief minUnion
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, T > minUnion(const Vec< N, T >& a, const Vec< N, T >& b)
{
    Vec< N, T > m;
    for (size_t i = 0; i < N; ++i)
    {
        (a.v[i] < b.v[i]) ? m.v[i] = a.v[i] : m.v[i] = b.v[i];
    }
    return m;
}

/**
 * @brief maxUnion
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, T > maxUnion(const Vec< N, T >& a, const Vec< N, T >& b)
{
    Vec< N, T > m;
    for (size_t i = 0; i < N; ++i)
    {
        (a.v[i] > b.v[i]) ? m.v[i] = a.v[i] : m.v[i] = b.v[i];
    }
    return m;
}

/**
 * @brief max
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline T max(const Vec< N, T >& a)
{
    T m = a.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        if(a.v[i] > m)
        {
            m = a.v[i];
        }
    }
    return m;
}

/**
 * @brief dot
 * @param a
 * @param b
 * @return
 */
template< uint64_t N, class T >
inline T dot(const Vec< N, T >& a, const Vec< N, T >& b)
{
    T d = a.v[0] * b.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        d += a.v[i] * b.v[i];
    }
    return d;
}

/**
 * @brief rotate
 * Counter-clockwise rotation.
 * @param a
 * @param angle
 * @return
 */
template< class T >
inline Vec< 2, T > rotate(const Vec< 2, T >& a, float angle)
{
    T c = cos(angle);
    T s = sin(angle);
    return Vec< 2, T >(c * a[0] - s * a[1], s * a[0] + c * a[1]);
}

/**
 * @brief perp
 * Counter-clockwise rotation by 90 degrees.
 * @param a
 * @return
 */
template< class T >
inline Vec< 2, T > perp(const Vec< 2, T >& a)
{
    return Vec< 2, T >(-a.v[1], a.v[0]);
}

/**
 * @brief cross
 * @param a
 * @param b
 * @return
 */
template< class T >
inline T cross(const Vec< 2, T >& a, const Vec< 2, T >& b)
{
    return a.v[0] * b.v[1] - a.v[1] * b.v[0];
}

/**
 * @brief cross
 * @param a
 * @param b
 * @return
 */
template< class T >
inline Vec< 3, T > cross(const Vec< 3, T >& a, const Vec< 3, T >& b)
{
    return Vec< 3, T >(a.v[1] * b.v[2] - a.v[2] * b.v[1],
                       a.v[2] * b.v[0] - a.v[0] * b.v[2],
                       a.v[0] * b.v[1] - a.v[1] * b.v[0]);
}

/**
 * @brief sdfTriple
 * @param a
 * @param b
 * @param c
 * @return
 */
template< class T >
inline T sdfTriple(const Vec< 3, T >& a,
                   const Vec< 3, T >& b,
                   const Vec< 3, T >& c)
{
    return  a.v[0] * (b.v[1] * c.v[2] - b.v[2] * c.v[1]) +
            a.v[1] * (b.v[2] * c.v[0] - b.v[0] * c.v[2]) +
            a.v[2] * (b.v[0] * c.v[1] - b.v[1] * c.v[0]);
}

/**
 * @brief sdfhash
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline uint32_t sdfHash(const Vec< N, T >& a)
{
    uint32_t h = a.v[0];
    for (size_t i = 1; i < N; ++i)
    {
        h = sdfHash(h ^ a.v[i]);
    }
    return h;
}

/**
 * @brief assign
 * @param a
 * @param a0
 * @param a1
 */
template< uint64_t N, class T >
inline void assign(const Vec< N, T >& a, T& a0, T& a1)
{ 
    assert(N == 2);
    a0 = a.v[0];
    a1 = a.v[1];
}

/**
 * @brief assign
 * @param a
 * @param a0
 * @param a1
 * @param a2
 */
template< uint64_t N, class T >
inline void assign(const Vec< N, T >& a, T& a0, T& a1, T& a2)
{ 
    assert(N == 3);
    a0 = a.v[0];
    a1 = a.v[1];
    a2 = a.v[2];
}

/**
 * @brief assign
 * @param a
 * @param a0
 * @param a1
 * @param a2
 * @param a3
 */
template< uint64_t N, class T >
inline void assign(const Vec< N, T >& a, T& a0, T& a1, T& a2, T& a3)
{ 
    assert(N == 4);
    a0 = a.v[0];
    a1 = a.v[1];
    a2 = a.v[2];
    a3 = a.v[3];
}

/**
 * @brief assign
 * @param a
 * @param a0
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 */
template< uint64_t N, class T >
inline void assign(const Vec< N, T >& a,
                   T& a0, T& a1, T& a2, T& a3, T& a4, T& a5)
{ 
    assert(N == 6);
    a0 = a.v[0];
    a1 = a.v[1];
    a2 = a.v[2];
    a3 = a.v[3];
    a4 = a.v[4];
    a5 = a.v[5];
}

/**
 * @brief round
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, int32_t > round(const Vec< N, T >& a)
{ 
    Vec< N, int32_t > rounded;
    for (size_t i = 0; i < N; ++i)
    {
        rounded.v[i] = lround(a.v[i]);
    }
    return rounded;
}

/**
 * @brief floor
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, int32_t > floor(const Vec< N, T >& a)
{ 
    Vec< N, int32_t > rounded;
    for (size_t i = 0; i < N; ++i)
    {
        rounded.v[i] = static_cast<int32_t>(floor(a.v[i]));
    }
    return rounded;
}

/**
 * @brief ceil
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, int32_t > ceil(const Vec< N, T >& a)
{ 
    Vec< N, int32_t > rounded;
    for (size_t i = 0; i < N; ++i)
    {
        rounded.v[i] = static_cast<int32_t>(ceil(a.v[i]));
    }
    return rounded;
}

/**
 * @brief fabs
 * @param a
 * @return
 */
template< uint64_t N, class T >
inline Vec< N, T > fabs(const Vec< N, T >& a)
{ 
    Vec< N, T > result;
    for (size_t i = 0; i < N; ++i)
    {
        result.v[i] = fabs(a.v[i]);
    }
    return result;
}

/**
 * @brief minMax
 * @param x0
 * @param x1
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void minMax(const Vec< N, T >& x0,
                   const Vec< N, T >& x1,
                   Vec< N, T >& xMin,
                   Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
        minMax(x0.v[i], x1.v[i], xMin.v[i], xMax.v[i]);
}

/**
 * @brief minMax
 * @param x0
 * @param x1
 * @param x2
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void minMax(const Vec< N, T >& x0,
                   const Vec< N, T >& x1,
                   const Vec< N, T >& x2,
                   Vec< N, T >& xMin, Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
    {
        minMax(x0.v[i], x1.v[i], x2.v[i], xMin.v[i], xMax.v[i]);
    }
}

/**
 * @brief minMax
 * @param x0
 * @param x1
 * @param x2
 * @param x3
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void minMax(const Vec< N, T >& x0,
                   const Vec< N, T >& x1,
                   const Vec< N, T >& x2,
                   const Vec< N, T >& x3,
                   Vec< N, T >& xMin, Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
    {
        minMax(x0.v[i], x1.v[i], x2.v[i], x3.v[i], xMin.v[i], xMax.v[i]);
    }
}

/**
 * @brief minMax
 * @param x0
 * @param x1
 * @param x2
 * @param x3
 * @param x4
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void minMax(const Vec< N, T >& x0,
                   const Vec< N, T >& x1,
                   const Vec< N, T >& x2,
                   const Vec< N, T >& x3,
                   const Vec< N, T >& x4,
                   Vec< N, T >& xMin,
                   Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
    {
        minMax(x0.v[i], x1.v[i], x2.v[i], x3.v[i], x4.v[i],
               xMin.v[i], xMax.v[i]);
    }
}

/**
 * @brief minMax
 * @param x0
 * @param x1
 * @param x2
 * @param x3
 * @param x4
 * @param x5
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void minMax(const Vec< N, T >& x0,
                   const Vec< N, T >& x1,
                   const Vec< N, T >& x2,
                   const Vec< N, T >& x3,
                   const Vec< N, T >& x4,
                   const Vec< N, T >& x5,
                   Vec< N, T >& xMin,
                   Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
    {
        minMax(x0.v[i], x1.v[i], x2.v[i],
               x3.v[i], x4.v[i], x5.v[i],
               xMin.v[i], xMax.v[i]);
    }
}

/**
 * @brief updateMinMax
 * @param x
 * @param xMin
 * @param xMax
 */
template< uint64_t N, class T >
inline void updateMinMax(const Vec< N, T >& x,
                         Vec< N, T >& xMin,
                         Vec< N, T >& xMax)
{
    for (size_t i = 0; i < N; ++i)
        updateMinMax(x[i], xMin[i], xMax[i]);
}

}
