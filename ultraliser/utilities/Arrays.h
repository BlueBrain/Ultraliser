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

#include <data/common/BitArray.h>

namespace Ultraliser
{
namespace Array
{

/**
 * @brief convertStringTo8UIArray
 * Converts a string array to uint8_t array.
 * @param stringArray
 * Input string array that contains binary data.
 * @return
 * uint8_t array
 */
uint8_t* convertStringTo8UIArray(const std::string &stringArray);

/**
 * @brief convertStringTo16UIArray
 * Converts a string array to uint8_t array.
 * @param stringArray
 * Input string array that contains binary data.
 * @return
 * An uint16_t array
 */
uint16_t* convertStringTo16UIArray(const std::string &stringArray);

/**
 * @brief convertStringTo32UIArray
 * Converts a string array to uint16_t array.
 * @param stringArray
 * Input string array that contains binary data.
 * @return
 * An uint32_t array
 */
uint32_t* convertStringTo32UIArray(const std::string &stringArray);

/**
 * @brief convertStringTo64UIArray
 * Converts a string array to uint32_t array.
 * @param stringArray
 * Input string array that contains binary data.
 * @return
 * An uint64_t array
 */
uint64_t* convertStringTo64UIArray(const std::string &stringArray);

/**
 * @brief convertStringToBitArray
 * Converts a string array to uint64_t array.
 * @param stringArray
 * Input string array that contains binary data.
 * @return
 * A pointer to the resulting BitArray.
 */
BitArray* convertStringToBitArray(const std::string &stringArray, const size_t &numberBits);

}
}
