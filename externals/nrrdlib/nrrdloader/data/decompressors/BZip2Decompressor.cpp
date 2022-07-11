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

#include "BZip2Decompressor.h"

#include <stdexcept>

#include <bzlib.h>

namespace libNRRD
{
std::string BZip2Decompressor::decompress(std::string input) const
{
    const auto inputData = input.data();
    const auto inputSize = input.size();

    std::string result;
    result.resize(inputSize * 2);
    auto outputData = result.data();
    auto outputSize = static_cast<unsigned int>(result.capacity());

    auto bzlibResult = BZ2_bzBuffToBuffDecompress(outputData, &outputSize, inputData, inputSize, 0, 0);

    switch (bzlibResult)
    {
    case BZ_CONFIG_ERROR:
        throw std::runtime_error("Incompatibility with libbzip2 library error");
    case BZ_MEM_ERROR:
        throw std::runtime_error("Unnable to allocate enough memory for the decompression");
    case BZ_OUTBUFF_FULL:
        throw std::runtime_error("Output buffer is not big enough to hold the decompressed data");
    case BZ_PARAM_ERROR:
        throw std::runtime_error("Unknown error during decompression");
    default:
        break;
    }

    result.resize(outputSize);

    return result;
}
}
