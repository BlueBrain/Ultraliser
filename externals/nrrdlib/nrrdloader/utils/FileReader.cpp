/***************************************************************************************************
 * Copyright (c) 2015 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * All rights reserved. Do not distribute without permission.
 * Author(s):
 *      Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of Brayns <https://github.com/BlueBrain/Brayns>
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301 USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include "FileReader.h"

#include <filesystem>
#include <fstream>
#include <stdexcept>

namespace libNRRD
{
std::string FileReader::read(const std::string &filename)
{
    if (!std::filesystem::exists(filename))
    {
        throw std::invalid_argument("File " + filename + " does not exists");
    }

    if (!std::filesystem::is_regular_file(filename))
    {
        throw std::invalid_argument("File " + filename + " is not a regular file");
    }

    std::ifstream stream(filename, std::ios::ate);
    if (!stream.is_open())
    {
        throw std::runtime_error("Cannot read file '" + filename + "'");
    }
    auto size = stream.tellg();
    std::string data;
    data.resize(size);
    stream.seekg(0);
    stream.read(&data[0], size);
    return data;
}
} // namespace brayns