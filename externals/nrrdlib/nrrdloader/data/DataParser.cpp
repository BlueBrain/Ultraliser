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

#include "DataParser.h"

#include "../utils/FileReader.h"

#include "decoders/DecoderTable.h"
#include "decompressors/DecompressorTable.h"

#include <filesystem>

namespace
{
class DataFilePaths
{
public:
    static std::vector<std::string> buildFixed(const libNRRD::NRRDHeader &header)
    {
        const auto filePath = std::filesystem::path(header.filePath);
        const auto basePath = filePath.parent_path();

        // Copy, manipulate, return
        auto fileList = *header.dataFiles;

        for (auto &file : fileList)
        {
            const auto filePath = basePath / std::filesystem::path(file);
            file = filePath.string();
        }

        return fileList;
    }
};
class DataFileReader
{
public:
    static std::string read(const std::vector<std::string> &fileList)
    {
        std::string result;
        for (const auto &file : fileList)
        {
            result += libNRRD::FileReader::read(file);
        }
        return result;
    }
};

class DataContentParser
{
public:
    static std::unique_ptr<libNRRD::IDataMangler> parse(
        const libNRRD::NRRDHeader &header,
        std::string content)
    {
        const auto format = header.encoding;

        const auto decompressor = libNRRD::DecompressorTable::getDecompressor(format);
        const auto decompressedContent = decompressor->decompress(std::move(content));

        const auto decoder = libNRRD::DecoderTable::getDecoder(format);
        return decoder->decode(header, decompressedContent);
    }
};
}

namespace libNRRD
{
std::unique_ptr<IDataMangler> DataParser::parse(const NRRDHeader &header, std::string_view content)
{
    if (header.dataFiles)
    {
        const auto files = DataFilePaths::buildFixed(header);
        content = DataFileReader::read(files);
    }

    return DataContentParser::parse(header, std::string(content));
}
}
