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

#ifndef ULTRALISER_UTILITIES_PARSERS_H
#define ULTRALISER_UTILITIES_PARSERS_H

#include <data/meshes/simple/primitives/Primitives.h>

namespace Ultraliser
{
namespace Parsers
{

/**
 * @brief getIntegerValue
 * Gets an integer value from a given string.
 * @param line
 * @return
 */
int getIntegerValue(const std::string line);

/**
 * @brief getStringValue
 * Gets a string value from a given line.
 * @param line
 * @return
 */
std::string getStringValue(const std::string line);

/**
 * @brief parseFloat
 * @param fileStream
 * @return
 */
float parseFloat(std::ifstream& fileStream);

/**
 * @brief parseVector3f
 * @param fileStream
 * @param v
 */
void parseVector3f(std::ifstream& fileStream, Vector3f &v);

/**
 * @brief parseVector2f
 * @param line
 * @return
 */
Vector2f parseVector2f(std::string line);

/**
 * @brief parseVector3i
 * @param line
 * @return
 */
Vec3i_64 parseVector3i(std::string line);

/**
 * @brief parseVector3f
 * @param line
 * @return
 */
Vector3f parseVector3f(std::string line);

/**
 * @brief parseFaceData
 * @param line
 * @return
 */
Vec3i_64 parseFaceData(std::string line);

/**
 * @brief parseFaces
 * @param line
 * @param actualTriangles
 */
void parseFaces(std::string line, Triangles& actualTriangles);

}
}

#endif // ULTRALISER_UTILITIES_PARSERS_H
