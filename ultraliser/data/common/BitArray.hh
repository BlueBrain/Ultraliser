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

/// Make CHAR_BIT 8 if it's not defined in limits.h
#ifndef CHAR_BIT
#define CHAR_BIT 8
#endif

/// Position of bit within character
#define BIT_CHAR(bit)         ((bit) / CHAR_BIT)

/// Array index for character containing bit
#define BIT_IN_CHAR(bit)      (1 << (CHAR_BIT - 1 - ((bit)  % CHAR_BIT)))

/// Number of characters required to contain number of bits
#define BITS_TO_CHARS(bits)   ((((bits) - 1) / CHAR_BIT) + 1)

/// Most significant bit in a character
#define MS_BIT                (1 << (CHAR_BIT - 1))
