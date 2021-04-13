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

#ifndef ULTRALISER_ARGUMENTS_ARGUMENT_H
#define ULTRALISER_ARGUMENTS_ARGUMENT_H

#include <arguments/ArgumentEnums.hh>
#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The Argument class
 */
class Argument
{
public:

    /**
     * @brief Argument
     * @param name
     * @param defaultValue
     * @param help
     * @param type
     * @param presence
     */
    Argument(const std::string name,
             const ARGUMENT_TYPE type,
             const std::string help,
             const ARGUMENT_PRESENCE presence = ARGUMENT_PRESENCE::OPTIONAL,
             const std::string defaultValue = NO_DEFAULT_VALUE);

    /**
     * @brief isValid
     * @param argvString
     * @return
     */
    bool isValid(const std::string &argvString);

    /**
     * @brief evaluate
     * @param argvString
     */
    void evaluate(const std::string &argvString);

    /**
     * @brief getIntegerValue
     * @param valid
     * @return
     */
    int32_t getIntegerValue(bool& valid);

    /**
     * @brief getUnsignedIntegerValue
     * @param valid
     * @return
     */
    uint32_t getUnsignedIntegerValue(bool& valid);

    /**
     * @brief getFloatValue
     * @param valid
     * @return
     */
    float getFloatValue(bool& valid);

    /**
     * @brief getStringValue
     * @param valid
     * @return
     */
    std::string getStringValue(bool &valid);

    /**
     * @brief getBoolValue
     * @param valid
     * @return
     */
    bool getBoolValue(bool &valid);

    /**
     * @brief displayValue
     */
    void displayValue();

    /**
     * @brief showHelp
     */
    void showHelp();

    /**
     * @brief getCopy
     * Returns a pointer to a copy of this object.
     * @return
     */
    Argument* getCopy();

    /**
     * @brief getName
     * @return
     */
    std::string getName() const;

private:

    /**
     * @brief _formatHelp
     * @return
     */
    std::string _formatHelp();

private:

    /**
     * @brief _name
     */
    const std::string _name;

    /**
     * @brief _defaultValue
     */
    const std::string _defaultValue;

    /**
     * @brief _help
     */
    const std::string _help;

    /**
     * @brief _type    // Exit
    exit(0);
     */
    const ARGUMENT_TYPE _type;

    /**
     * @brief _presence
     */
    const ARGUMENT_PRESENCE _presence;

    /**
     * @brief _optionalSet
     * A flag to indicate that the optional arguments are set.
     */
    bool _optionalSet;

    /**
     * @brief _value
     */
    std::string _value;
};

/**
 * @brief Arguments
 */
typedef std::vector<Argument*> Arguments;

}

#endif // ULTRALISER_ARGUMENTS_ARGUMENT_H
