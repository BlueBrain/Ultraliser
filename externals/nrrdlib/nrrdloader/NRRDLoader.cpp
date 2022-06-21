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

#include "NRRDLoader.h"

#include "data/DataParser.h"
#include "header/HeaderParser.h"
#include "utils/FileReader.h"

namespace libNRRD
{
NRRDVolume NRRDLoader::load(const std::string &path)
{
    const auto fileData = FileReader::read(path);
    auto fileContent = std::string_view(fileData);

    auto header = HeaderParser::parse(path, fileContent);

    auto data = DataParser::parse(header, fileContent);

    return NRRDVolume{std::move(header), std::move(data)};
}
}
