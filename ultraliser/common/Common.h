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

#ifndef ULTRALISER_COMMON_COMMON_H
#define ULTRALISER_COMMON_COMMON_H

#include <common/Headers.hh>
#include <common/Defines.h>
#include <common/Progress.h>
#include <common/Logging.h>
#include <common/Enums.h>
#include <common/OpenMP.hh>

static int64_t AXES[ AXES_COUNT ][ 3 ] =
{
    { -1,  0,  0 }, // -X
    {  1,  0,  0 }, // +X
    {  0, -1,  0 }, // -Y
    {  0,  1,  0 }, // +Y
    {  0,  0, -1 }, // -Z
    {  0,  0,  1 }  // +Z
};

inline void p_swap(void **a, void **b) {void *t = *a; *a = *b; *b = t;}

typedef unsigned char	UBYTE;
typedef   signed char	 BYTE;
typedef unsigned short UINT16;
typedef   signed short	INT16;
#endif // ULTRALISER_COMMON_COMMON_H
