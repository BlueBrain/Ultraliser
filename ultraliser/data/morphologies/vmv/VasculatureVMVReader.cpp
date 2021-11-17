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
        // const uint64_t vertexIndex = S2UI(tokens[0]);

        // Sample attributes
        const float x = S2F(tokens[1]);
        const float y = S2F(tokens[2]);
        const float z = S2F(tokens[3]);
        float r = S2F(tokens[4]);
        if (r < 0.0001)
            r = 0.1;

        // Construct the sample
        Sample* sample = new Sample(Vector3f(x, y, z), r);
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
        // LOG_INFO("%s", cleanLine.c_str());

        // Section index is at index 0 of the tokens
        uint64_t sectionIndex = S2UI(tokens[0]);
        Section* section = new Section(sectionIndex);

        // Fill the section with its sampless
        for (uint64_t i = 1; i < tokens.size(); ++i)
        {
            //LOG_INFO("%s", tokens[i].c_str());
            // Get the section sample and construct it
            const uint64_t sampleIndex = S2UI(tokens[i]);
            section->addSample(_samples[sampleIndex - 1]);
        }

        // Add the section to the list of the sections in the morphology
        _sections.push_back(section);
    }
}

}
