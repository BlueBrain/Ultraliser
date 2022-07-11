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

#include "Args.h"

namespace Ultraliser
{

Args::Args(const int argc, const char **argv, const std::string &help)
{
    // Argument parser
    _parser = new Ultraliser::ArgumentParser(argc, argv, help);
}

Args::~Args()
{
    delete _parser;
}

void Args::addArgument(Ultraliser::Argument* argument)
{
    _parser->addArgument(argument->getCopy());
}

void Args::parse()
{
    // Show help just in case there is any error before the actual parsing
    _parser->showHelpIfNeeded();

    // Parse the arguments
    _parser->parse();
}

std::string Args::getStringValue(Ultraliser::Argument* argument)
{
    return _parser->getStringValue(argument);
}

int32_t Args::getIntegrValue(Ultraliser::Argument* argument)
{
    return _parser->getIntegrValue(argument);
}

uint32_t Args::getUnsignedIntegrValue(Ultraliser::Argument* argument)
{
    return _parser->getUnsignedIntegrValue(argument);
}

float Args::getFloatValue(Ultraliser::Argument* argument)
{
    return _parser->getFloatValue(argument);
}

bool Args::getBoolValue(Ultraliser::Argument* argument)
{
    return _parser->getBoolValue(argument);
}

}
