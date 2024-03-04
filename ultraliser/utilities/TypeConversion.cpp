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

#include "TypeConversion.h"

namespace Ultraliser
{
namespace Utilities
{

template< class T >
float int2float(T input)
{
    return static_cast< float >(input);
}

template< class T >
float uint2float(T input)
{
    return static_cast< float >(input);
}

template< class T >
double int2double(T input)
{
     return static_cast< double >(input);
}

template< class T >
double uint2double(T input)
{
     return static_cast< double >(input);
}

float double2float(double input)
{
    return static_cast< float >(input);
}

double float2double(float input)
{
    return static_cast< double >(input);
}

template< class T >
T string2int(std::string input)
{
    return static_cast< T >(atoi(input.c_str()));
}

float string2float(std::string input)
{
    return static_cast< float >(atof(input.c_str()));
}

double string2double(std::string input)
{
    return static_cast< double >(atof(input.c_str()));
}

template< class T >
std::string number2string(T input)
{
    std::stringstream stream;
    stream << input;
    return stream.str();
}

}
}
