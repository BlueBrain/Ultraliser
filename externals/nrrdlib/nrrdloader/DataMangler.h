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

#pragma once

#include <cstdint>
#include <cstdlib>
#include <vector>

namespace libNRRD
{
class IDataMangler
{
public:
    virtual ~IDataMangler() = default;

    virtual std::vector<char> asChars() const = 0;
    virtual std::vector<uint8_t> asBytes() const = 0;
    virtual std::vector<int16_t> asShorts() const = 0;
    virtual std::vector<uint16_t> asUnsingedShorts() const = 0;
    virtual std::vector<int32_t> asIntegers() const = 0;
    virtual std::vector<uint32_t> asUnsignedIntegers() const = 0;
    virtual std::vector<int64_t> asLongs() const = 0;
    virtual std::vector<uint64_t> asUnsignedLongs() const = 0;
    virtual std::vector<float> asFloats() const = 0;
    virtual std::vector<double> asDoubles() const = 0;

    virtual bool isTypeSigned() const noexcept = 0;
    virtual size_t getTypeSize() const noexcept = 0;
    virtual size_t getNumElements() const noexcept = 0;
};

template<typename T>
class DataMangler final : public IDataMangler
{
public:
    DataMangler(std::vector<T> data)
        : _data(std::move(data))
    {
        static_assert(std::is_scalar_v<T>, "DataMangler only supports scalar values");
    }

    std::vector<char> asChars() const noexcept override
    {
        return _as<char>();
    }

    std::vector<uint8_t> asBytes() const noexcept override
    {
        return _as<uint8_t>();
    }

    std::vector<int16_t> asShorts() const noexcept override
    {
        return _as<int16_t>();
    }

    std::vector<uint16_t> asUnsingedShorts() const noexcept override
    {
        return _as<uint16_t>();
    }

    std::vector<int32_t> asIntegers() const noexcept override
    {
        return _as<int32_t>();
    }

    std::vector<uint32_t> asUnsignedIntegers() const noexcept override
    {
        return _as<uint32_t>();
    }

    std::vector<int64_t> asLongs() const noexcept override
    {
        return _as<int64_t>();
    }

    std::vector<uint64_t> asUnsignedLongs() const noexcept override
    {
        return _as<uint64_t>();
    }

    std::vector<float> asFloats() const noexcept override
    {
        return _as<float>();
    }

    std::vector<double> asDoubles() const noexcept override
    {
        return _as<double>();
    }

    bool isTypeSigned() const noexcept override
    {
        return std::is_signed_v<T>;
    }

    size_t getTypeSize() const noexcept override
    {
        return sizeof(T);
    }

    size_t getNumElements() const noexcept override
    {
        return _data.size();
    }

private:
    template<typename To>
    std::vector<To> _as() const noexcept
    {
        if constexpr (std::is_same_v<T, To>)
        {
            return _data;
        }

        std::vector<To> result(_data.size());
        for (size_t i = 0; i < _data.size(); ++i)
        {
            result[i] = static_cast<To>(_data[i]);
        }
        return result;
    }

private:
    std::vector<T> _data;
};
}
