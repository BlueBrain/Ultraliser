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

#include "ArgumentParser.h"

namespace Ultraliser
{

/**
 * @brief The Args class
 */
class Args
{
public:

    /**
     * @brief Args
     * Constructor
     * @param argc
     * Arguments count.
     * @param argv
     * Arguments
     */
    Args(const int argc, const char** argv, const std::string &help = "");
    ~Args();

    /**
     * @brief addArgument
     * Adds a new argument.
     * @param argument
     * A new argument.
     */
    void addArgument(Ultraliser::Argument* argument);

    /**
     * @brief getStringValue
     * Returns the string value that corresponds to the argument.
     * @param argument
     * @return
     */
    std::string getStringValue(Ultraliser::Argument* argument);

    /**
     * @brief getIntegrValue
     * @param argument
     * @return
     */
    int32_t getIntegrValue(Ultraliser::Argument* argument);

    /**
     * @brief getUnsignedIntegrValue
     * @param argument
     * @return
     */
    size_t getUnsignedIntegrValue(Ultraliser::Argument* argument);

    /**
     * @brief getFloatValue
     * @param argument
     * @return
     */
    float getFloatValue(Ultraliser::Argument* argument);

    /**
     * @brief getBoolValue
     * @param argument
     * @return
     */
    bool getBoolValue(Ultraliser::Argument* argument);

    /**
     * @brief parse
     * Parse the command line options.
     */
    void parse();

private:

    /**
     * @brief _parser
     * Command lines arguments parser.
     */
     ArgumentParser* _parser;
};

}
