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

#include "DecoderTable.h"

#include "AsciiDecoder.h"
#include "HexDecoder.h"
#include "RawDecoder.h"

namespace libNRRD
{
std::unique_ptr<IDecoder> DecoderTable::getDecoder(NRRDEncoding encoding)
{
    if (encoding == NRRDEncoding::Ascii)
    {
        return std::make_unique<AsciiDecoder>();
    }

    if (encoding == NRRDEncoding::Hex)
    {
        return std::make_unique<HexDecoder>();
    }

    return std::make_unique<RawDecoder>();
}
}
