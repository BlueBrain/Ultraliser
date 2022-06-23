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

#include <common/Common.h>
#include <utilities/TypeConversion.h>
#include <math.h>

namespace Ultraliser
{

/**
 * @brief sqr
 * @param x
 * @return
 */
template< class T >
inline T sqr(const T& x)
{
    return x * x;
}

template< class T >
inline T clamp(T a, T lower, T upper)
{
    if (a < lower)
    {
        return lower;
    }
    else if (a > upper)
    {
        return upper;
    }
    else
    {
        return a;
    }
}

/**
 * @brief inbound
 * @param ii
 * @param jj
 * @param array
 * @return
 */
template <typename T>
bool inbound(int32_t ii, int32_t jj,
             const std::vector<std::vector< T > >* array)
{
    bool ret = (ii >= 0) && (ii < array->size()) && (jj >= 0)
               && (jj < (*array)[0].size());
    return ret;
}

/**
 * @brief inBound
 * @param a
 * @param lower
 * @param upper
 * @return
 */
template< class T >
inline bool inBound(T a, T lower, T upper)
{
    if(a < lower)
    {
        return false;
    }
    else if(a >= upper)
    {
        return false;
    }
    else
    {
        return true;
    }
}

/**
 * @brief cube
 * @param x
 * @return
 */
template< class T >
inline T cube(const T& x)
{
    return x * x  * x;
}

/**
 * @brief min
 * @param a1
 * @param a2
 * @param a3
 * @return
 */
template< class T >
inline T min(T a1, T a2, T a3)
{
    return min(a1, min(a2, a3));
}

/**
 * @brief min
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @return
 */
template< class T >
inline T min(T a1, T a2, T a3, T a4)
{
    return min(min(a1, a2), min(a3, a4));
}

/**
 * @brief min
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @return
 */
template< class T >
inline T min(T a1, T a2, T a3, T a4, T a5)
{
    return min(min(a1, a2), min(a3, a4), a5);
}

/**
 * @brief min
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @param a6
 * @return
 */
template< class T >
inline T min(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return min(min(a1, a2), min(a3, a4), min(a5, a6));
}

/**
 * @brief max
 * @param a1
 * @param a2
 * @param a3
 * @return
 */
template< class T >
inline T max(T a1, T a2, T a3)
{
    return max(a1, max(a2, a3));
}

/**
 * @brief max
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @return
 */
template< class T >
inline T max(T a1, T a2, T a3, T a4)
{
    return max(max(a1, a2), max(a3, a4));
}

/**
 * @brief max
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @return
 */
template< class T >
inline T max(T a1, T a2, T a3, T a4, T a5)
{
    return max(max(a1, a2), max(a3, a4),  a5);
}

/**
 * @brief max
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @param a6
 * @return
 */
template< class T >
inline T max(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return max(max(a1, a2), max(a3, a4),  max(a5, a6));
}

/**
 * @brief minMax
 * @param a1
 * @param a2
 * @param aMin
 * @param aMax
 */
template< class T >
inline void minMax(T a1, T a2, T& aMin, T& aMax)
{
    if(a1 < a2)
    {
        aMin = a1;
        aMax = a2;
    }
    else
    {
        aMin = a2;
        aMax = a1;
    }
}

/**
 * @brief minMax
 * @param a1
 * @param a2
 * @param a3
 * @param aMin
 * @param aMax
 */
template< class T >
inline void minMax(T a1, T a2, T a3, T& aMin, T& aMax)
{
    if(a1 < a2)
    {
        if(a1 < a3)
        {
            aMin = a1;
            if(a2 < a3) aMax = a3;
            else aMax = a2;
        }
        else
        {
            aMin = a3;
            if(a1 < a2) aMax = a2;
            else aMax = a1;
        }
    }
    else
    {
        if(a2 < a3)
        {
            aMin = a2;
            if(a1 < a3) aMax = a3;
            else aMax = a1;
        }
        else
        {
            aMin = a3;
            aMax = a1;
        }
    }
}

/**
 * @brief minMax
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param aMin
 * @param aMax
 */
template< class T >
inline void minMax(T a1, T a2, T a3, T a4, T& aMin, T& aMax)
{
    if(a1 < a2)
    {
        if(a3 < a4)
        {
            aMin = min(a1, a3);
            aMax = max(a2, a4);
        }
        else
        {
            aMin = min(a1, a4);
            aMax = max(a2, a3);
        }
    }
    else
    {
        if(a3 < a4){
            aMin = min(a2, a3);
            aMax = max(a1, a4);
        }
        else
        {
            aMin = min(a2, a4);
            aMax = max(a1, a3);
        }
    }
}

/**
 * @brief minMax
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @param aMin
 * @param aMax
 */
template< class T >
inline void minMax(T a1, T a2, T a3, T a4, T a5, T& aMin, T& aMax)
{
    aMin = min(a1, a2, a3, a4, a5);
    aMax = max(a1, a2, a3, a4, a5);
}

/**
 * @brief minMax
 * @param a1
 * @param a2
 * @param a3
 * @param a4
 * @param a5
 * @param a6
 * @param aMin
 * @param aMax
 */
template< class T >
inline void minMax(T a1, T a2, T a3, T a4, T a5, T a6, T& aMin, T& aMax)
{
    aMin = min(a1, a2, a3, a4, a5, a6);
    aMax = max(a1, a2, a3, a4, a5, a6);
}

/**
 * @brief update_minmax
 * @param a1
 * @param aMin
 * @param aMax
 */
template< class T >
inline void update_minmax(T a1, T& aMin, T& aMax)
{
    if(a1 < aMin)
        aMin = a1;
    else if(a1 > aMax)
        aMax = a1;
}

/**
 * @brief sort
 * @param a
 * @param b
 * @param c
 */
template< class T >
inline void sort(T &a, T &b, T &c)
{
    T temp;
    if(a < b)
    {
        if(a < c)
        {
            if(c < b)
            {
                temp = c;
                c = b;
                b = temp;
            }
        }
        else
        {
            temp = c;
            c = b;
            b = a;
            a = temp;
        }
    }
    else
    {
        if(b < c)
        {
            if(a < c)
            {
                temp = b;
                b = a;
                a = temp;
            }
            else
            {
                temp = b;
                b = c;
                c = a;
                a = temp;
            }
        }
        else
        {
            temp = c;
            c = a ;
            a = temp;
        }
    }
}

/**
 * @brief smoothStep
 * @param r
 * @return
 */
template< class T >
inline T smoothStep (T r)
{
    if(r < 0) return 0;
    else if(r > 1) return 1;
    return r * r * r * (10 + r * (-15 + r  * 6));
}

/**
 * @brief smoothStep
 * @param r
 * @param rLower
 * @param rUpper
 * @param valueLower
 * @param valueUpper
 * @return
 */
template< class T >
inline T smoothStep(T r, T rLower, T rUpper, T valueLower, T valueUpper)
{
    return valueLower + smoothStep((r - rLower) / (rUpper - rLower)) *
                                                  (valueUpper - valueLower);
}

/**
 * @brief ramp
 * @param r
 * @return
 */
template< class T >
inline T ramp(T r)
{
    return smoothStep((r + 1) / 2) * 2 - 1;
}

/**
 * @brief roundUpToPowerOfTwo
 * @param n
 * @return
 */
inline uint64_t roundUpToPowerOfTwo(uint64_t n)
{
    int32_t exponent = 0;
    --n;
    while(n)
    {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

/**
 * @brief roundDownToPowerOfTwo
 * @param n
 * @return
 */
inline uint64_t roundDownToPowerOfTwo(uint64_t n)
{
    int32_t exponent=0;
    while(n > 1)
    {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

/**
 * @brief randHash
 * Transforms even the sequence 0,1,2,3,... into reasonably good random numbers
 * @param seed
 * @return
 */
inline uint64_t randHash(uint64_t seed)
{
    uint64_t i = (seed^0xA3C59AC3u) * 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    return i;
}

/**
 * @brief unhash
 * The inverse of randhash.
 * @param h
 * @return
 */
inline uint64_t unhash(uint64_t h)
{
    h *= 340573321u;
    h ^= (h>>16);
    h *= 340573321u;
    h ^= (h>>16);
    h *= 340573321u;
    h ^= 0xA3C59AC3u;
    return h;
}

/**
 * @brief randHash_double
 * @param seed
 * @return repeatable stateless pseudo-random number in [0,1]
 */
inline double randHash_double(uint64_t seed)
{
    return randHash(seed) / I2D(UINT_MAX);
}

/**
 * @brief randhash_float
 * @param seed
 * @return
 */
inline float randhash_float(uint64_t seed)
{
    return randHash(seed) / I2F(UINT_MAX);
}

/**
 * @brief randhashd
 * @param seed
 * @param a
 * @param b
 * @return repeatable stateless pseudo-random number in [a,b]
 */
inline double randHash_double(uint64_t seed, double a, double b)
{
    return (b - a) * randHash(seed) / I2D(UINT_MAX) + a;
}

/**
 * @brief randhash_float
 * @param seed
 * @param a
 * @param b
 * @return
 */
inline float randhash_float(uint64_t seed, float a, float b)
{
    return ((b - a) * randHash(seed) / I2F(UINT_MAX) + a);
}

/**
 * @brief intlog2
 * @param x
 * @return
 */
inline int32_t intlog2(int32_t x)
{
    int32_t exp = -1;
    while(x)
    {
        x >>= 1;
        ++exp;
    }

    return exp;
}

/**
 * @brief getBaryCentric
 * @param x
 * @param i
 * @param f
 * @param iLow
 * @param iHigh
 */
template< class T >
inline void getBaryCentric(T x, int32_t& i, T& f, int32_t iLow, int32_t iHigh)
{
    T s = std::floor(x);
    i = static_cast<int32_t>(s);
    if(i < iLow)
    {
        i = iLow;
        f = 0;
    }
    else if(i > iHigh - 2)
    {
        i = iHigh - 2;
        f = 1;
    }
    else
    {
        f = static_cast<T>(x - s);
    }
}

/**
 * @brief lerp
 * @param value0
 * @param value1
 * @param f
 * @return
 */
template< class S, class T >
inline S lerp(const S& value0, const S& value1, T f)
{
    return (1 - f) * value0 + f * value1;
}

/**
 * @brief biLerp
 * @param v00
 * @param v10
 * @param v01
 * @param v11
 * @param fx
 * @param fy
 * @return
 */
template< class S, class T >
inline S biLerp(const S& v00, const S& v10,
                const S& v01, const S& v11,
                T fx, T fy)
{ 
    return lerp(lerp(v00, v10, fx), lerp(v01, v11, fx), fy);
}

/**
 * @brief trilerp
 * @param v000
 * @param v100
 * @param v010
 * @param v110
 * @param v001
 * @param v101
 * @param v011
 * @param v111
 * @param fx
 * @param fy
 * @param fz
 * @return
 */
template< class S, class T >
inline S trilerp(const S& v000, const S& v100,
                 const S& v010, const S& v110,
                 const S& v001, const S& v101,
                 const S& v011, const S& v111,
                 T fx, T fy, T fz)
{
    return lerp(biLerp(v000, v100, v010, v110, fx, fy),
                biLerp(v001, v101, v011, v111, fx, fy),
                fz);
}

/**
 * @brief quadLerp
 * @param v0000
 * @param v1000
 * @param v0100
 * @param v1100
 * @param v0010
 * @param v1010
 * @param v0110
 * @param v1110
 * @param v0001
 * @param v1001
 * @param v0101
 * @param v1101
 * @param v0011
 * @param v1011
 * @param v0111
 * @param v1111
 * @param fx
 * @param fy
 * @param fz
 * @param ft
 * @return
 */
template< class S, class T >
inline S quadLerp(const S& v0000, const S& v1000,
                  const S& v0100, const S& v1100,
                  const S& v0010, const S& v1010,
                  const S& v0110, const S& v1110,
                  const S& v0001, const S& v1001,
                  const S& v0101, const S& v1101,
                  const S& v0011, const S& v1011,
                  const S& v0111, const S& v1111,
                  T fx, T fy, T fz, T ft)
{
    return lerp(trilerp(v0000, v1000, v0100, v1100,
                        v0010, v1010, v0110, v1110,
                        fx, fy, fz),
                trilerp(v0001, v1001, v0101, v1101,
                        v0011, v1011, v0111, v1111,
                        fx, fy, fz), ft);
}

/**
 * @brief quadraticBsplineWeights
 * @param f
 *The value of f should be between 0 and 1, with f=0.5 corresponding to balanced
 * weighting between w0 and w2.
 * @param w0
 * @param w1
 * @param w2
 */
template< class T >
inline void quadraticBsplineWeights(T f, T& w0, T& w1, T& w2)
{
    w0 = T(0.50) * sqr(f - 1);
    w1 = T(0.75) - sqr(f - T(0.5));
    w2 = T(0.50) * sqr(f);
}

/**
 * @brief cubicInterpolationWeights
 * @param f
 * @param wneg1
 * @param w0
 * @param w1
 * @param w2
 */
template< class T >
inline void cubicInterpolationWeights(T f, T& wneg1, T& w0, T& w1, T& w2)
{
    T f2(f * f), f3(f2 * f);
    wneg1 = -T(1.0 / 3) * f + T(1.0/ 2) * f2 - T(1.0 / 6) * f3;
    w0 = 1 - f2 + T(1.0 / 2) * (f3 - f);
    w1 = f + T(1.0 / 2) * (f2 - f3);
    w2 = T(1.0 / 6) * (f3 -f);
}

/**
 * @brief cubicInterpolation
 * @param valueNeg1
 * @param value0
 * @param value1
 * @param value2
 * @param f
 * @return
 */
template< class S, class T >
inline S cubicInterpolation(const S& valueNeg1,
                            const S& value0,
                            const S& value1,
                            const S& value2, T f)
{
    T wneg1, w0, w1, w2;
    cubicInterpolationWeights(f, wneg1, w0, w1, w2);
    return wneg1 * valueNeg1 + w0 * value0 + w1 * value1 + w2 * value2;
}

/**
 * @brief zero
 * @param v
 */
template< class T >
void zero(std::vector< T >& v)
{
    for(int64_t i = v.size() - 1; i >= 0; --i)
    {
        v[i] = 0;
    }
}

/**
 * @brief absMax
 * @param v
 * @return
 */
template< class T >
T absMax(const std::vector< T >& v)
{
    T m = 0;
    for(int64_t i = v.size() - 1; i >= 0; --i)
    {
        if(std::fabs(v[i]) > m)
        {
            m = std::fabs(v[i]);
        }
    }

    return m;
}

/**
 * @brief contains
 * @param a
 * @param e
 * @return
 */
template< class T >
bool contains(const std::vector< T >& a, T e)
{
    for(uint64_t i = 0; i < a.size(); ++i)
    {
        if(a[i] == e)
        {
            return true;
        }
    }

    return false;
}

/**
 * @brief addunique
 * @param a
 * @param e
 */
template< class T >
void addunique(std::vector< T >& a, T e)
{
    for(uint64_t i = 0; i < a.size(); ++i)
    {
        if(a[i] == e) return;
    }

    a.push_back(e);
}

/**
 * @brief insert
 * @param a
 * @param index
 * @param e
 */
template< class T >
void insert(std::vector< T >& a, uint64_t index, T e)
{
    a.push_back(a.back());

    for(uint64_t i = a.size() - 1; i > index; --i)
    {
        a[i] = a[i - 1];
    }

    a[index] = e;
}

/**
 * @brief erase
 * @param a
 * @param index
 */
template< class T >
void erase(std::vector< T >& a, uint64_t index)
{
    for(uint64_t i = index; i < a.size() - 1; ++i)
    {
        a[i] = a[i + 1];
    }

    a.pop_back();
}

/**
 * @brief erase_swap
 * @param a
 * @param index
 */
template< class T >
void erase_swap(std::vector< T >& a, uint64_t index)
{
    for(uint64_t i = index; i < a.size() - 1; ++i)
        swap(a[i], a[i + 1]);
    a.pop_back();
}

/**
 * @brief eraseUnordered
 * @param a
 * @param index
 */
template< class T >
void eraseUnordered(std::vector< T >& a, uint64_t index)
{
    a[index]=a.back();
    a.pop_back();
}

/**
 * @brief eraseUnorderedSwap
 * @param a
 * @param index
 */
template< class T >
void eraseUnorderedSwap(std::vector< T >& a, uint64_t index)
{
    swap(a[index], a.back());
    a.pop_back();
}

/**
 * @brief findAndEraseUnordered
 * @param a
 * @param doomedElement
 */
template< class T >
void findAndEraseUnordered(std::vector< T >& a, const T& doomedElement)
{
    for(uint64_t i = 0; i < a.size(); ++i)
    {
        if(a[i] == doomedElement)
        {
            eraseUnordered(a, i);
            return;
        }
    }
}

/**
 * @brief replaceOnce
 * @param a
 * @param oldElement
 * @param newElement
 */
template< class T >
void replaceOnce(std::vector< T >& a,
                 const T& oldElement,
                 const T& newElement)
{
    for(uint64_t i = 0; i < a.size(); ++i)
    {
        if(a[i] == oldElement)
        {
            a[i] = newElement;
            return;
        }
    }
}

/**
 * @brief isZero
 * Checks if a given floating point number is zero or not.
 * @param number
 * @return
 */
template< class T >
bool isZero(T number)
{
    if (number < std::numeric_limits<T>::epsilon())
        return true;
    return false;
}

/**
 * @brief isEqual
 * Checks if two floating point numbers are equal or not.
 * @param a
 * @param b
 * @return
 */
template< class T >
bool isEqual(T a, T b)
{
    T difference = std::abs(a - b);
    if (isZero(difference))
        return true;
    return false;
}

}

#endif // ULTRALISER_MATH_UTILITIES_H
