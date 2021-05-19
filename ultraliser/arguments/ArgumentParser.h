/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_ARGUMENTS_ARGUMENT_PARSER_H
#define ULTRALISER_ARGUMENTS_ARGUMENT_PARSER_H

#include "Argument.h"

namespace Ultraliser
{
class ArgumentParser
{
public:

    /**
     * @brief ArgumentParser
     * @param argc
     * @param argv
     * @param applicationHelp
     */
    ArgumentParser(const int argc,
                   const char** argv,
                   const std::string &applicationHelp = "");
    ~ArgumentParser();

    /**
     * @brief parse
     */
    void parse();

    /**
     * @brief updateArguments
     */
    void updateArguments();

    /**
     * @brief displayValues
     */
    void displayValues();

    /**
     * @brief displayArgumentsString
     */
    void displayArgumentsString();

    /**
     * @brief addArgument
     * @param argument
     */
    void addArgument(Argument* argument);

    /**
     * @brief showHelpIfNeeded
     */
    void showHelpIfNeeded();

    /**
     * @brief getStringValue
     * @param argument
     * @return
     */
    std::string getStringValue(Argument* argument);

    /**
     * @brief getIntegrValue
     * @param argument
     * @return
     */
    int getIntegrValue(Argument* argument);

    /**
     * @brief getUnsignedIntegrValue
     * @param argument
     * @return
     */
    uint32_t getUnsignedIntegrValue(Argument* argument);

    /**
     * @brief getFloatValue
     * @param argument
     * @return
     */
    float getFloatValue(Argument* argument);

    /**
     * @brief getBoolValue
     * @param argument
     * @return
     */
    bool getBoolValue(Argument* argument);

private:

    /**
     * @brief _cleanInputArguments
     */
    void _cleanInputArguments();

    /**
     * @brief _validateArguments
     */
    void _validateArguments();

    /**
     * @brief _validateCurrentArguments
     */
    void _validateCurrentArguments();

    /**
     * @brief _evaluateArguments
     */
    void _evaluateArguments();

private:

    /**
     * @brief _argc
     */
    const int _argc;

    /**
     * @brief _argv
     */
    const char** _argv;

    /**
     * @brief _argvString
     */
    std::string _argvString;

    /**
     * @brief _arguments
     */
    Arguments _arguments;

    /**
     * @brief _applicationHelp
     */
    std::string _applicationHelp;

    bool _showHelp;
};

}

#endif // ULTRALISER_ARGUMENTS_ARGUMENT_PARSER_H
