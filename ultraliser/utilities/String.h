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
namespace String
{

/**
 * @brief removeSubstring
 * Removes a substring from the given string.
 * @param mainString
 * The main string.
 * @param subString
 * A substring that will be removed from the original string.
 */
void removeSubstring(std::string& mainString, const std::string subString);

/**
 * @brief replaceSubstring
 * Replaces a substring in a given string with another one.
 * @param mainString
 * The main string.
 * @param from
 * The original substring.
 * @param to
 * The replaced substring.
 * @return
 * True or False to indicate the correction of the function.
 */
bool replaceSubstring(std::string& mainString,
                      const std::string &from,
                      const std::string &to);

/**
 * @brief subStringFound
 * Checks if a substring is found in a string or not.
 * @param mainString
 * The main string.
 * @param subString
 * The substring.
 * @return
 * True or False.
 */
bool subStringFound(std::string& mainString, const std::string &subString);

/**
 * @brief split
 * @param string
 * @param delimiter
 * @return
 */
std::vector< std::string > split(std::string string, char delimiter);

/**
 * @brief formatStringToMultiLine
 * @param input
 * @param limit
 * @return
 */
std::string formatStringToMultiLine(const std::string &input,
                                    const uint64_t limit = 60,
                                    const bool &aligned = false);

bool sameString(char *a, char *b);

bool sameString(std::string s1, std::string s2);

/**
 * @brief removeExtraSpaces
 * Removes extra spaces from a given string.
 * @param inputString
 * String to be cleaned.
 */
void removeExtraSpaces(std::string &inputString);

/**
 * @brief toLower
 * @param inputString
 */
void toLower(std::string &inputString);

/**
 * @brief toUpper
 * @param inputString
 */
void toUpper(std::string &inputString);

}
}
