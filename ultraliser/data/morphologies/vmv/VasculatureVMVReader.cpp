/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include "VasculatureVMVReader.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

VasculatureVMVReader::VasculatureVMVReader(const std::string &vmvMorphologyFilePath)
    : _vmvMorphologyFile(vmvMorphologyFilePath)
{
    _readAttributes();

    // Read the samples
    _readSamples();

    // Read the connectivity
    _readConnectivity();
}

VasculatureMorphology* VasculatureVMVReader::getMorphology()
{
    return new VasculatureMorphology(_samples, _sections);
}

void VasculatureVMVReader::_readAttributes()
{
    std::ifstream stream;
    stream.open(_vmvMorphologyFile.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the morphology file [ %s ]", _vmvMorphologyFile.c_str());
    }

    std::vector< std::string > attributesBlock;
    bool collectData = false;

    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
    {
        if (collectData)
        {
            attributesBlock.push_back(line);
        }

        // Start collecting data
        if (Ultraliser::String::subStringFound(line, std::string("$PARAM_BEGIN")))
        {
            collectData = true;
        }

        // Break to close the file
        if (Ultraliser::String::subStringFound(line, std::string("$PARAM_END")))
        {
            collectData = false;
            break;
        }
    }

    // Close the file stream
    stream.close();

    for (std::string line: attributesBlock)
    {
        if (Ultraliser::String::subStringFound(line, std::string("NUM_VERTS")))
        {
            // Split the line into tokens split by a space
            std::vector< std::string > tokens;
            String::removeExtraSpaces(line);
            std::istringstream iss(line);
            copy(std::istream_iterator< std::string >(iss),
                 std::istream_iterator< std::string >(),
                 std::back_inserter(tokens));

            _numberVerts = S2UI(tokens[1]);
        }

        if (Ultraliser::String::subStringFound(line, std::string("NUM_STRANDS")))
        {
            // Split the line into tokens split by a space
            std::vector< std::string > tokens;
            std::istringstream iss(line);
            copy(std::istream_iterator< std::string >(iss),
                 std::istream_iterator< std::string >(),
                 std::back_inserter(tokens));

            _numberStrands = S2UI(tokens[1]);
        }

        if (Ultraliser::String::subStringFound(line, std::string("NUM_ATTRIB_PER_VERT")))
        {
            // Split the line into tokens split by a space
            std::vector< std::string > tokens;
            std::istringstream iss(line);
            copy(std::istream_iterator< std::string >(iss),
                 std::istream_iterator< std::string >(),
                 std::back_inserter(tokens));

            _numberAttributesPerVertex = S2UI(tokens[1]);
        }
    }
}

void VasculatureVMVReader::_readSamples()
{
    std::ifstream stream;
    stream.open(_vmvMorphologyFile.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the morphology file [ %s ]", _vmvMorphologyFile.c_str());
    }

    std::vector< std::string > verticesBlock;
    bool collectData = false;

    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
    {
        // Break to close the file
        if (Ultraliser::String::subStringFound(line, std::string("$VERT_LIST_END")))
        {
            break;
        }

        // We can safely collect the data
        if (collectData)
        {
            verticesBlock.push_back(line);
        }

        // Start collecting data
        if (Ultraliser::String::subStringFound(line, std::string("$VERT_LIST_BEGIN")))
        {
            collectData = true;
        }
    }

    // Close the file stream
    stream.close();

    // Reads the vertices from the vertex block
    for (std::string line: verticesBlock)
    {
        // Split the line into tokens split by a space
        std::vector< std::string > tokens;
        std::istringstream iss(line);
        copy(std::istream_iterator< std::string >(iss),
             std::istream_iterator< std::string >(),
             std::back_inserter(tokens));

        // Vertex index is ignored
        const size_t vertexIndex = static_cast< size_t >(S2UI(tokens[0]));

        // Sample attributes
        const float x = S2F(tokens[1]);
        const float y = S2F(tokens[2]);
        const float z = S2F(tokens[3]);
        float r = S2F(tokens[4]);
        if (r < 0.1)
            r = 0.1;

        r *= 0.8;

        // Construct the sample
        Sample* sample = new Sample(Vector3f(x, y, z), r, PROCESS_TYPE::VASCULATURE, vertexIndex);
        _samples.push_back(sample);
    }
}

void VasculatureVMVReader::_readConnectivity()
{
    std::ifstream stream;
    stream.open(_vmvMorphologyFile.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the morphology file [ %s ]", _vmvMorphologyFile.c_str());
    }

    std::vector< std::string > strandsBlock;
    bool collectData = false;

    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
    {
        // Break to close the file
        if (Ultraliser::String::subStringFound(line, std::string("$STRANDS_LIST_END")))
        {
            break;
        }

        // It is safe to collect the data
        if (collectData)
        {
            strandsBlock.push_back(line);
        }

        // Start collecting data
        if (Ultraliser::String::subStringFound(line, std::string("$STRANDS_LIST_BEGIN")))
        {
            collectData = true;
        }
    }

    // Close the file stream
    stream.close();

    // Reads the vertices from the vertex block
    for (std::string line: strandsBlock)
    {
        // Split the line into tokens split by a space
        std::vector< std::string > tokens;
        String::removeExtraSpaces(line);
        std::istringstream iss(line);
        copy(std::istream_iterator< std::string >(iss),
             std::istream_iterator< std::string >(),
             std::back_inserter(tokens));

        // Section index is at index 0 of the tokens
        size_t sectionIndex = static_cast<size_t> (S2UI(tokens[0]));
        Section* section = new Section(sectionIndex, PROCESS_TYPE::VASCULATURE);

        // Fill the section with its sampless
        for (size_t i = 1; i < tokens.size(); ++i)
        {
            // Get the section sample and construct it
            const size_t sampleIndex = static_cast< size_t >(S2UI(tokens[i]));
            section->addSample(_samples[sampleIndex - 1]);
        }

        // Add the section to the list of the sections in the morphology
        _sections.push_back(section);
    }
}

}
