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

#include "Argument.h"
#include <utilities/Utilities.h>

namespace Ultraliser
{

Argument::Argument(const std::string name,
                   const ARGUMENT_TYPE type,
                   const std::string help,
                   const ARGUMENT_PRESENCE presence,
                   const std::string defaultValue)
    : _name(name)
    , _defaultValue(defaultValue)
    , _help(help)
    , _type(type)
    , _presence(presence)
{
    // Use the default value
    _value = defaultValue;
}

Argument::~Argument()
{
    /// EMPTY
}

Argument* Argument::getCopy()
{
    return new Argument(_name, _type, _help, _presence, _defaultValue);
}


size_t getArgumentCount(std::string argvString, std::string argument)
{
    size_t count = 0;
    const size_t step = argument.size();

    size_t position = 0;

    while ((position=argvString.find(argument, position)) !=std::string::npos)
    {
      position +=step;
      ++count ;
    }

    return count;
}

void Argument::validateCurrentOptionalArguments(const std::string &argvString)
{
    // Get the number of times this argument was added in the command line
    size_t count = getArgumentCount(argvString, _name);

    // If count is 0 and the argument is optional, then it is always valid
    if (count == 0 && _presence == ARGUMENT_PRESENCE::OPTIONAL)
    {
        // The optional argument is not set, use the default value
        _optionalSet = false;
    }

    // If the count is 1, then it is valid
    if (count == 1 && _presence == ARGUMENT_PRESENCE::OPTIONAL)
    {
        // The optional argument is set, use the given value
        _optionalSet = true;
    }
}

bool Argument::isValid(const std::string &argvString)
{
    // If the argument is mandatory and not given, then exit
    if (_presence == ARGUMENT_PRESENCE::MANDATORY)
    {
        if (getArgumentCount(argvString, _name) == 0)
        {
            LOG_ERROR("\t Argument %s is not optional", _name.c_str());
        }
    }

    // Get the number of times this argument was added in the command line
    size_t count = getArgumentCount(argvString, _name);

    // If count is 0 and the argument is optional, then it is always valid
    if (count == 0 && _presence == ARGUMENT_PRESENCE::OPTIONAL)
    {
        // The optional argument is not set, use the default value
        _optionalSet = false;

        // Valid
        return true;
    }

    // If the count is 1, then it is valid
    if (count == 1 && _presence == ARGUMENT_PRESENCE::OPTIONAL)
    {
        // The optional argument is set, use the given value
        _optionalSet = true;

        // Valid
        return true;
    }

    // If the argument is declared more than once
    if (count > 1)
    {
        LOG_ERROR("Argument %s was already defined.", _name.c_str());
        return false;
    }

    // The argument is valid
    return true;
}

void Argument::evaluate(const std::string &argvString)
{
    if (_type == ARGUMENT_TYPE::BOOL)
    {
        size_t count = getArgumentCount(argvString, _name);
        if (count == 1)
            _value = "TRUE";
        else
            _value = "FALSE";
    }
    else
    {
        if (_presence == ARGUMENT_PRESENCE::OPTIONAL && !_optionalSet)
            _value = _defaultValue;
        else
        {
            // Find the rest of the string that starts from the name of the
            // argument
            std::string post = argvString.substr(argvString.find(_name));

            // Convert the list to a vector
            std::stringstream stream(post);
            std::istream_iterator<std::string> begin(stream);
            std::istream_iterator<std::string> end;
            std::vector<std::string> stringVector(begin, end);

            // Ensure that the argument has a value
            if (stringVector.size() == 1)
            {
                LOG_ERROR("The values of the argument [%s] is not provided by the user!",
                          stringVector[0].c_str());
            }
            else
            {
                // Pick the second element
                _value = stringVector[1];
            }

            // Remove all the spaces
            String::removeSubstring(_value, " ");
        }
    }
}

int32_t Argument::getIntegerValue(bool& valid)
{
    if (_type == ARGUMENT_TYPE::INTEGER)
    {
        valid = true;
        return atoi(_value.c_str());
    }

    valid = false;
    return 0;
}

uint32_t Argument::getUnsignedIntegerValue(bool& valid)
{
    if (_type == ARGUMENT_TYPE::INTEGER)
    {
        valid = true;
        return static_cast< uint32_t >(atoi(_value.c_str()));
    }

    valid = false;
    return 0;
}

float Argument::getFloatValue(bool& valid)
{
    if (_type == ARGUMENT_TYPE::FLOAT)
    {
        valid = true;
        return static_cast<float>(atof(_value.c_str()));
    }

    valid = false;
    return 0.0;
}

std::string Argument::getStringValue(bool &valid)
{
    if (_type == ARGUMENT_TYPE::STRING)
    {
        valid = true;
        return _value;
    }

    valid = false;
    return _value;
}

bool Argument::getBoolValue(bool &valid)
{
    if (_type == ARGUMENT_TYPE::BOOL)
    {
        valid = true;

        if (_value == "TRUE")
            return true;
        else
            return false;
    }
    return false;
}

void Argument::displayValue()
{
    printf("[ %s ]:[ %s ] \n", _name.c_str(), _value.c_str());
}

std::string Argument::_formatHelp()
{
    return String::formatStringToMultiLine(_help, 60);
}

void Argument::showHelp()
{
    printf("\t %s: \n\t %s\n\n", _name.c_str(), _formatHelp().c_str());
}

std::string Argument::getName() const
{
    return _name;
}

}
