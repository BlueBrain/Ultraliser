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

#include <utilities/Parsers.h>
#include <common/Common.h>

namespace Ultraliser
{
namespace Parsers
{

int getIntegerValue( const std::string line )
{
    std::vector< std::string > tokens;
    std::istringstream iss( line );
    copy( std::istream_iterator< std::string >( iss ),
          std::istream_iterator< std::string >( ),
          std::back_inserter( tokens ));

    return std::atoi( tokens[1].c_str());
}

std::string getStringValue( const std::string line )
{
    std::vector< std::string > tokens;
    std::istringstream iss( line );
    copy( std::istream_iterator< std::string >( iss ),
          std::istream_iterator< std::string >( ),
          std::back_inserter( tokens ));

    return tokens[1];
}

float parseFloat(std::ifstream& fileStream)
{
    char buffer[sizeof(float)];
    fileStream.read(buffer, 4);

    float* floatPtr = (float*) buffer;
    return *floatPtr;
}

void parseVector3f(std::ifstream& fileStream, Vector3f& v)
{
    v.x() = parseFloat(fileStream);
    v.y() = parseFloat(fileStream);
    v.z() = parseFloat(fileStream);
}

Vector2f parseVector2f(std::string line)
{
    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    Vector2f value;
    value.x() = S2F(tokens[0]);
    value.y() = S2F(tokens[1]);
    return value;
}

Vec3i_64 parseVector3i(std::string line)
{
    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    Vec3i_64 value;
    value[0] = std::atoi(tokens[0].c_str());
    value[1] = std::atoi(tokens[1].c_str());
    value[2] = std::atoi(tokens[2].c_str());
    return value;
}

Vector3f parseVector3f(std::string line)
{
    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    Vector3f value;
    value.x() = S2F(tokens[0]);
    value.y() = S2F(tokens[1]);
    value.z() = S2F(tokens[2]);
    return value;
}

Vec3i_64 parseFaceData(std::string line)
{
    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    Vec3i_64 value;
    value[0] = std::atoi(tokens[1].c_str());
    value[1] = std::atoi(tokens[2].c_str());
    value[2] = std::atoi(tokens[3].c_str());
    return value;
}

void parseFaces(std::string line, Triangles& actualTriangles)
{
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));


    // Read the number of vertices forming the face
    size_t numberVertices = std::atoi(tokens[0].c_str());
    for (size_t i = 0; i < numberVertices - 2; i++)
    {
        Triangle triangle;
        triangle[0] = std::atoi(tokens[1].c_str());
        triangle[1] = std::atoi(tokens[i + 2].c_str());
        triangle[2] = std::atoi(tokens[i + 3].c_str());
        actualTriangles.push_back(triangle);
    }
}

}
}

