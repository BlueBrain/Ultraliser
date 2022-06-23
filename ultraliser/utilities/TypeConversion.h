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

namespace Ultraliser
{
namespace Utilities
{

/**
 * @brief int2float
 * @param input
 * @return
 */
template <class T >
float int2float(T input);

/**
 * @brief uint2float
 * @param input
 * @return
 */
template <class T >
float uint2float(T input);

/**
 * @brief int2double
 * @param input
 * @return
 */
template <class T >
double int2double(T input);

/**
 * @brief uint2double
 * @param input
 * @return
 */
template <class T >
double uint2double(T input);

/**
 * @brief double2float
 * @param input
 * @return
 */
float double2float(double input);

/**
 * @brief float2double
 * @param input
 * @return
 */
double float2double(float input);

/**
 * @brief string2int8
 * @param input
 * @return
 */
template < class T >
T string2int(std::string input);

/**
 * @brief string2uint
 * @param input
 * @return
 */
template < class T >
T string2uint(std::string input);

/**
 * @brief string2float
 * @param input
 * @return
 */
float string2float(std::string input);

/**
 * @brief string2double
 * @param input
 * @return
 */
double string2double(std::string input);

// Specialization
template<> float int2float< int >(int);
template<> float int2float< int8_t >(int8_t);
template<> float int2float< int16_t >(int16_t);
template<> float int2float< int32_t >(int32_t);
template<> float int2float< int64_t >(int64_t);

template<> float uint2float< uint32_t >(uint32_t);
template<> float uint2float< uint8_t >(uint8_t);
template<> float uint2float< uint16_t >(uint16_t);
template<> float uint2float< uint32_t >(uint32_t);
template<> float uint2float< uint64_t >(uint64_t);

template<> double int2double< int >(int);
template<> double int2double< int8_t >(int8_t);
template<> double int2double< int16_t >(int16_t);
template<> double int2double< int32_t >(int32_t);
template<> double int2double< int64_t >(int64_t);

template<> double int2double< uint32_t >(uint32_t);
template<> double int2double< uint8_t >(uint8_t);
template<> double int2double< uint16_t >(uint16_t);
template<> double int2double< uint32_t >(uint32_t);
template<> double int2double< uint64_t >(uint64_t);

template<> int string2int< int >(std::string);
template<> int8_t string2int< int8_t >(std::string);
template<> int16_t string2int< int16_t >(std::string);
template<> int32_t string2int< int32_t >(std::string);
template<> int64_t string2int< int64_t >(std::string);

template<> uint32_t string2uint< uint32_t >(std::string);
template<> uint8_t string2uint< uint8_t >(std::string);
template<> uint16_t string2uint< uint16_t >(std::string);
template<> uint32_t string2uint< uint32_t >(std::string);
template<> uint64_t string2uint< uint64_t >(std::string);
template<> size_t string2uint< size_t >(std::string);

}
}
