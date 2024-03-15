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

#include <common/Logging.hh>
#include <string>

// Safe conversion of any integert type to float
#define I2F(NUMBER) static_cast< float >(NUMBER)

// Safe conversion of any unsigned integer number to float
#define UI2F(NUMBER) static_cast< float >(NUMBER)

// Safe conversion of any integert type to double
#define I2D(NUMBER) static_cast< double >(NUMBER)

// Safe conversion of any unsigned integer number to double
#define UI2D(NUMBER) static_cast< double >(NUMBER)

// fSafe conversion of double to float
#define D2F(NUMBER) static_cast< float >(NUMBER)

// Safe conversion of float to double
#define F2D(NUMBER) static_cast< double >(NUMBER)

// Safe conversion of string to any integer number
#define S2I(STRING_NUMBER) static_cast< int64_t >(atoi(STRING_NUMBER.c_str()))

// Safe conversion of string to any integer number
#define S2UI(STRING_NUMBER) static_cast< uint64_t >(atoi(STRING_NUMBER.c_str()))

// Safe conversion of string to any float number
#define S2F(STRING_NUMBER) static_cast< float >(atof(STRING_NUMBER.c_str()))

// Safe conversion of string to any double number
#define S2D(STRING_NUMBER) static_cast< double >(atof(STRING_NUMBER.c_str()))

// Safe conversion of any interger to int
#define I2I8(NUMBER) static_cast< int8_t >(NUMBER)
#define I2I16(NUMBER) static_cast< int16_t >(NUMBER)
#define I2I32(NUMBER) static_cast< int32_t >(NUMBER)
#define I2I64(NUMBER) static_cast< int64_t >(NUMBER)

// Safe conversion of any unsigned interger to int
#define UI2I8(NUMBER) static_cast< int8_t >(NUMBER)
#define UI2I16(NUMBER) static_cast< int16_t >(NUMBER)
#define UI2I32(NUMBER) static_cast< int32_t >(NUMBER)
#define UI2I64(NUMBER) static_cast< int64_t >(NUMBER)

// Safe conversion of any unsigned interger to unsigned int
#define UI2UI8(NUMBER) static_cast< uint8_t >(NUMBER)
#define UI2UI16(NUMBER) static_cast< uint16_t >(NUMBER)
#define UI2UI32(NUMBER) static_cast< uint32_t >(NUMBER)
#define UI2UI64(NUMBER) static_cast< uint64_t >(NUMBER)

// Safe conversion of any interger to unsigned int
#define I2UI8(NUMBER) static_cast< uint8_t >(NUMBER)
#define I2UI16(NUMBER) static_cast< uint16_t >(NUMBER)
#define I2UI32(NUMBER) static_cast< uint32_t >(NUMBER)
#define I2UI64(NUMBER) static_cast< uint64_t >(NUMBER)

// Safe conversion of float to any integer
#define F2I8(NUMBER) static_cast< int8_t >(NUMBER)
#define F2I16(NUMBER) static_cast< int16_t >(NUMBER)
#define F2I32(NUMBER) static_cast< int32_t >(NUMBER)
#define F2I64(NUMBER) static_cast< int64_t >(NUMBER)

// Safe conversion of float to any unsigned integer
#define F2UI8(NUMBER) static_cast< uint8_t >(NUMBER)
#define F2UI16(NUMBER) static_cast< uint16_t >(NUMBER)
#define F2UI32(NUMBER) static_cast< uint32_t >(NUMBER)
#define F2UI64(NUMBER) static_cast< uint64_t >(NUMBER)

// Safe conversion of double to any integer
#define D2I8(NUMBER) static_cast< int8_t >(NUMBER)
#define D2I16(NUMBER) static_cast< int16_t >(NUMBER)
#define D2I32(NUMBER) static_cast< int32_t >(NUMBER)
#define D2I64(NUMBER) static_cast< int64_t >(NUMBER)

// Safe conversion of double to any unsigned integer
#define D2UI8(NUMBER) static_cast< uint8_t >(NUMBER)
#define D2UI16(NUMBER) static_cast< uint16_t >(NUMBER)
#define D2UI32(NUMBER) static_cast< uint32_t >(NUMBER)
#define D2UI64(NUMBER) static_cast< uint64_t >(NUMBER)

/**
 * @brief log
 * @param logLevel
 * @param filePath
 * @param lineNumber
 * @param functionName
 * @param string
 */
void log(const LOG_LEVEL& logLevel,
          const std::string &filePath,
          const int& lineNumber,
          const std::string &functionName,
          const char* lstring, ...);


#define LOG_HEADER(ARG...)                                                                          \
    log(LOG_LEVEL::HEADER, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_STATUS(ARG...)                                                                          \
    log(LOG_LEVEL::STATUS, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_STATUS_IMPORTANT(ARG...)                                                                \
    log(LOG_LEVEL::GREEN_STATUS, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_DETAIL(ARG...)                                                                          \
    log(LOG_LEVEL::DETAIL, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_INFO(ARG...)                                                                            \
    log(LOG_LEVEL::INFO, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_DEBUG(ARG...)                                                                           \
    log(LOG_LEVEL::DEBUG, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_WARNING(ARG...)                                                                         \
    log(LOG_LEVEL::WARNING, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_SUCCESS(ARG...)                                                                         \
    log(LOG_LEVEL::SUCCESS, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_ERROR(ARG...)                                                                           \
    log(LOG_LEVEL::ERROR, __FILE__, __LINE__, __FUNCTION__, ARG)

#define LOG_EXIT(ARG...)                                                                            \
    log(LOG_LEVEL::EXIT, __FILE__, __LINE__, __FUNCTION__, ARG)

// Print the statistics
#define LOG_STATS(TIME) std::cout << "\t* Time     | " << TIME << " Seconds \n" << std::endl

/**
 * @brief logTitle
 * Logs the title of a specific section.
 * @param title
 * Section title.
 * @param titleLength
 * Title length normally, 80 chars.
 */
static inline void logTitle(std::string title, const uint64_t& titleLength = 80)
{
    int size = I2I32(std::string(title).size());
    int spaces = titleLength - size;
    std::string bar = "";
    for (size_t i = 0; i < titleLength + 2; i++) bar += "-";
    printf("\n%s\n", bar.c_str()); bar = " ";
    for (size_t i = 0; i < uint64_t(0.5 * spaces) - 2; i++) bar += "-";
    bar += "| " + std::string(title) + " |";
    for (size_t i = 0; i < uint64_t(0.5 * spaces) - 2; i++) bar += "-";
    printf("%s\n", bar.c_str());
    bar = "" ; for (size_t i = 0; i < titleLength + 2; i++) bar += "-";
    printf("%s\n\n", bar.c_str()); bar = "";
    fflush(stdout);
}

// Print the title
#define LOG_TITLE(STRING) logTitle(STRING)

// Print a conditioned message
#define VERBOSE_LOG(LOG_MESSAGE, CONDITION) ({ if (CONDITION) { LOG_MESSAGE; }})

