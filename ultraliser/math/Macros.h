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

#include <math/Constants.h>
#include <utilities/TypeConversion.h>

// Random Number
#define RAND (rand() / I2F(RAND_MAX))

// Radians to degrees and viceversa
#define DEG2RAD(DEGREE) (DEGREE * ULTRALISER_PIF / 180.f)
#define RAD2DEG(RADIANS) (RADIANS * 180.f / ULTRALISER_PIF)

// Vector cross product
#define CROSS(result, v1, v2)                                                                      \
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];                                                     \
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];                                                     \
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];

// Vector dot product
#define DOT(v1, v2)                                                                                \
    (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

// Vector subtraction, for translation
#define SUB(result, v1, v2)                                                                        \
    result[0] = v1[0] - v2[0];                                                                     \
    result[1] = v1[1] - v2[1];                                                                     \
    result[2] = v1[2] - v2[2];

// Min/Max finder
#define FIND_MIN_MAX(x0, x1, x2 , min, max)                                                        \
    min = max = x0;                                                                                \
    if (x1 < min) min = x1;                                                                        \
    if (x1 > max) max = x1;                                                                        \
    if (x2 < min) min = x2;                                                                        \
    if (x2 > max) max = x2;

// X-TESTS
#define AXIS_TEST_X01(a, b, fa, fb)                                                                \
    p0 = a * v0[Y] - b * v0[Z];                                                                    \
    p2 = a * v2[Y] - b * v2[Z];                                                                    \
    if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }                              \
    rad = fa * boxHalfSize[Y] + fb * boxHalfSize[Z];                                               \
    if (min > rad || max < -rad) return 0;

#define AXIS_TEST_X2(a, b, fa, fb)                                                                 \
    p0 = a * v0[Y] - b * v0[Z];                                                                    \
    p1 = a * v1[Y] - b * v1[Z];                                                                    \
    if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }                              \
    rad = fa * boxHalfSize[Y] + fb * boxHalfSize[Z];                                               \
    if (min > rad || max < -rad) return 0;

// Y-TESTS
#define AXIS_TEST_Y02(a, b, fa, fb)                                                                \
    p0 = -a * v0[X] + b * v0[Z];                                                                   \
    p2 = -a * v2[X] + b * v2[Z];                                                                   \
    if (p0 < p2) { min = p0; max = p2; } else { min = p2; max = p0; }                              \
    rad = fa * boxHalfSize[X] + fb * boxHalfSize[Z];                                               \
    if (min > rad || max < -rad) return 0;

#define AXIS_TEST_Y1(a, b, fa, fb)                                                                 \
    p0 = -a * v0[X] + b * v0[Z];                                                                   \
    p1 = -a * v1[X] + b * v1[Z];                                                                   \
    if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }                              \
    rad = fa * boxHalfSize[X] + fb * boxHalfSize[Z];                                               \
    if (min > rad || max <- rad) return 0;

// Z-TESTS
#define AXIS_TEST_Z12(a, b, fa, fb)                                                                \
    p1 = a * v1[X] - b * v1[Y];                                                                    \
    p2 = a * v2[X] - b * v2[Y];                                                                    \
    if (p2 < p1) { min = p2; max = p1; } else { min = p1; max = p2; }                              \
    rad = fa * boxHalfSize[X] + fb * boxHalfSize[Y];                                               \
    if (min > rad || max <- rad) return 0;

#define AXIS_TEST_Z0(a, b, fa, fb)                                                                 \
    p0 = a*v0[X] - b*v0[Y];                                                                        \
    p1 = a*v1[X] - b*v1[Y];                                                                        \
    if (p0 < p1) { min = p0; max = p1; } else { min = p1; max = p0; }                              \
    rad = fa * boxHalfSize[X] + fb * boxHalfSize[Y];                                               \
    if (min > rad || max < -rad) return 0;
