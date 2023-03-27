/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Nadir Roman Guerrero < nadir.romanguerrero@epfl.ch >
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

#include "OffsetTable.h"

#include <array>
#include <cassert>

namespace
{
class OffsetTableComputer
{
public:
    static std::array<uint8_t, 2048> compute()
    {
        auto result = std::array<uint8_t, 2048>();
        result.fill(Ultraliser::OffsetTable::noEntry);

        for (size_t i = 0; i < 256; ++i)
        {
            auto idx = i << 3;
            auto offset = 0ul;

            for (size_t j = 0; j < 8; ++j)
            {
                auto mask = (1ul << j);
                if ((i & mask))
                {
                    result[idx + j] = offset;
                    offset++;
                }
            }
        }

        return result;
    }
};

struct OffsetTableStorage
{
    inline static const std::array<uint8_t, 2048> data = OffsetTableComputer::compute();
};
}

namespace Ultraliser
{
uint8_t OffsetTable::of(size_t tableIndex)
{
    assert(tableIndex < 2048);
    return OffsetTableStorage::data[tableIndex];
}
}