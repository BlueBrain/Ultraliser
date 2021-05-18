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

#include "ArgumentParser.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

ArgumentParser::ArgumentParser(const int argc, const char** argv,
                               const std::string &applicationHelp) :
    _argc(argc), _argv(argv), _applicationHelp(applicationHelp)
{
    // Construct a string from argv
    for (int i = 1; i < _argc; i++)
    {
        // Just append it to the string
        _argvString += _argv[i];

        // Append a space for clarity
        _argvString += " ";
    }
}

ArgumentParser::~ArgumentParser()
{
    for (auto& argument : _arguments)
    {
        delete argument;
    }

    // Clean the std::vector
    _arguments.clear();
    _arguments.shrink_to_fit();
}

bool areBothSpaces(char lhs, char rhs)
{
    return (lhs == rhs) && (lhs == ' ');
}

void ArgumentParser::_cleanInputArguments()
{
    // Replace all the = with spaces for generality
    std::replace(_argvString.begin(), _argvString.end(), '=', ' ');

    // Replace multiple spaces with a single one
    std::string::iterator newEnd = std::unique(_argvString.begin(),
                                               _argvString.end(),
                                               areBothSpaces);
    _argvString.erase(newEnd, _argvString.end());
}

void ArgumentParser::displayArgumentsString()
{
    std::cout << _argvString << std::endl;
}

void ArgumentParser::addArgument(Argument* argument)
{
    _arguments.push_back(argument);
}

void ArgumentParser::_validateArguments()
{
    // Validate the arguments for replication
    for (auto argument: _arguments)
    {
        if (!argument->isValid(_argvString))
            LOG_ERROR("Invalid argument: %s", argument->getName().c_str());
    }

    // Split all the arguments into a list to make it easy to validate them
    std::vector< std::string> argvStringList = String::split(_argvString, ' ');

    // Validate the arguments for existence
    for (auto argumentString : argvStringList)
    {
        // If the argument string starts with --, then this is a given argument
        if (argumentString.find("--") == 0)
        {
            // Check if this given argument exists in the defined ones or not
            bool exists = false;
            for (auto argument: _arguments)
            {
                if (argumentString == argument->getName())
                {
                    exists = true;
                }
            }

            if (!exists)
            {
                // Invalid argument
                LOG_ERROR("Argument: %s is not defined", argumentString.c_str());
            }
        }
    }
}

void ArgumentParser::_evaluateArguments()
{
    for (auto argument: _arguments)
        argument->evaluate(_argvString);
}

void ArgumentParser::parse()
{
    _cleanInputArguments();
    _validateArguments();
    _evaluateArguments();
}

void ArgumentParser::displayValues()
{
    for (const auto& argument:_arguments)
        argument->displayValue();
}

bool subStringFoundX(std::string& str,
                     const std::string &subString)
{
    size_t start_pos = str.find(subString);
    if (start_pos == std::string::npos)
        return false;
    return true;
}

void ArgumentParser::showHelpIfNeeded()
{
    if (_argc == 1                              ||
        subStringFoundX(_argvString, "--help ") ||
        subStringFoundX(_argvString, "-h "))
    {
        // Show tha application help
        printf("\n\t* About %s \n\n %s \n\n",
               _argv[0], String::formatStringToMultiLine(_applicationHelp, 60, true).c_str());

        printf("\t* Arguments \n\n");
        for (const auto& argument: _arguments)
            argument->showHelp();
        exit(0);
    }
}

std::string ArgumentParser::getStringValue(Argument* input)
{
    for (const auto& argument: _arguments)
    {
        bool valid;
        if (argument->getName() == input->getName())
            return argument->getStringValue(valid);
    }

    return EMPTY;
}

int32_t ArgumentParser::getIntegrValue(Argument* input)
{
    for (const auto& argument: _arguments)
    {
        bool valid;
        if (argument->getName() == input->getName())
            return argument->getIntegerValue(valid);
    }

    return 0;
}

uint32_t ArgumentParser::getUnsignedIntegrValue(Argument* input)
{
    for (const auto& argument: _arguments)
    {
        bool valid;
        if (argument->getName() == input->getName())
            return argument->getUnsignedIntegerValue(valid);
    }

    return 0;
}

float ArgumentParser::getFloatValue(Argument* input)
{
    for (const auto& argument: _arguments)
    {
        bool valid;
        if (argument->getName() == input->getName())
            return argument->getFloatValue(valid);
    }

    return 0.0;
}

bool ArgumentParser::getBoolValue(Argument* input)
{
    for (const auto& argument: _arguments)
    {
        bool valid;
        if (argument->getName() == input->getName())
            return argument->getBoolValue(valid);
    }

    return false;
}

}
