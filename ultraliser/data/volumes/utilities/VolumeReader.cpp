#include "VolumeReader.h"
#include <common/Headers.hh>
#include <utilities/String.h>
#include <nrrdloader/NRRDLoader.h>

namespace Ultraliser
{

VolumeData* readHeaderFile(const std::string& filePath)
{
    // Store the data in a VolumeData object
    auto volumeData = new VolumeData();

    // Open the header file
    std::ifstream hdrFileStream(filePath.c_str());

    std::string format;
    hdrFileStream >> format;

    if (String::subStringFound(format, FORMAT_BIT))
    {
        volumeData-> type = VOLUME_TYPE::BIT;
    }

    else if (String::subStringFound(format, FORMAT_8UI))
    {
        volumeData-> type = VOLUME_TYPE::UI8;
    }

    else if (String::subStringFound(format, FORMAT_16UI))
    {
        volumeData-> type = VOLUME_TYPE::UI16;
    }

    else if (String::subStringFound(format, FORMAT_32UI))
    {
        volumeData-> type = VOLUME_TYPE::UI32;
    }

    else if (String::subStringFound(format, FORMAT_64UI))
    {
        volumeData-> type = VOLUME_TYPE::UI64;
    }

    else if (String::subStringFound(format, FORMAT_F32))
    {
        volumeData-> type = VOLUME_TYPE::F32;
    }

    else if (String::subStringFound(format, FORMAT_F64))
    {
        volumeData-> type = VOLUME_TYPE::F64;
    }
    else
    {
        LOG_ERROR("Unknown volume format!");
    }

    // Volume dimensions
    hdrFileStream >> volumeData->width;
    hdrFileStream >> volumeData->height;
    hdrFileStream >> volumeData->depth;

    // Close the stream
    hdrFileStream.close();

    // Return the volume data
    return volumeData;
}

std::vector< uint8_t > readRawFileToByteVector(const std::string& filePath)
{
    // Input file stream
    std::fstream fileStream(filePath.c_str());

    // If the file cannot be opened
    if(!fileStream)
    {
        LOG_ERROR("The file [%s] cannot be loaded", filePath.c_str());
    }

    return std::vector< uint8_t >(std::istreambuf_iterator<char>(fileStream),
                                  std::istreambuf_iterator<char>());
}

std::string readRawFile(const std::string& filePath)
{
    // Input file stream
    std::fstream fileStream(filePath.c_str());

    // A string stream that will be used to parse the volume file
    std::stringstream dataStream;
    dataStream << fileStream.rdbuf();

    // Parse till the end of the stream
    std::string token = dataStream.str();

    // Close the file
    fileStream.close();

    // Return the string
    return token;
}

NRRDVolumeData* readNRRDVolumeFile(const std::string& filePath)
{
    // Load the volume
    auto nrrdVolume = libNRRD::NRRDLoader::load(filePath);

    // Reference to the header
    auto &header = nrrdVolume.header;

    // Store the data in a VolumeData object
    auto volumeData = new NRRDVolumeData();

    // Volume dimensions
    volumeData->width = header.sizes[0];
    volumeData->height = header.sizes[1];
    volumeData->depth = header.sizes[2];

    // Volume center
    const auto spaceOrigin = header.spaceOrigin;
    volumeData->center.x() = spaceOrigin->at(0);
    volumeData->center.y() = spaceOrigin->at(1);
    volumeData->center.z() = spaceOrigin->at(2);

    // Volume scale
    const auto spaceDirections = header.spaceDirections;

    volumeData->scale.x() = spaceDirections->at(0).at(0);
    volumeData->scale.y() = spaceDirections->at(1).at(1);
    volumeData->scale.z() = spaceDirections->at(2).at(2);

    // Volume type
    if (header.type == libNRRD::NRRDType::UnsignedChar)
    {
        volumeData-> type = VOLUME_TYPE::UI8;
    }
    else if (header.type == libNRRD::NRRDType::UnsignedShort)
    {
        volumeData-> type = VOLUME_TYPE::UI16;
    }
    else if (header.type == libNRRD::NRRDType::UnsignedInt)
    {
        volumeData-> type = VOLUME_TYPE::UI32;
    }
    else if (header.type == libNRRD::NRRDType::UnsignedLong)
    {
        volumeData-> type = VOLUME_TYPE::UI64;
    }
    else if (header.type == libNRRD::NRRDType::Float)
    {
        volumeData-> type = VOLUME_TYPE::F32;
    }
    else if (header.type == libNRRD::NRRDType::Double)
    {
        volumeData-> type = VOLUME_TYPE::F64;
    }
    else
    {
        LOG_ERROR("Undefined volume format");
    }

    // Store the data
    volumeData->data = nrrdVolume.data.release();

    // Return the volume data
    return volumeData;
}

UltraliserVolumeData *readUltraliserVolumeFile(const std::string &filePath)
{
    // Input file stream
    std::fstream fileStream(filePath.c_str());

    // A string stream that will be used to parse the volume file
    std::stringstream dataStream;
    dataStream << fileStream.rdbuf();

    // Parsing the file token-by-token
    std::string token;

    // Verify if reading the header is done or not
    auto isHeaderDone = false;

    // Store the data in a VolumeData object
    auto volumeData = new UltraliserVolumeData();

    // Parse the file
    while(1)
    {
        // Parse till the end of the line
        std::getline(dataStream, token, '\n');

        // Parse the type of the file
        if (String::subStringFound(token, std::string("format")))
        {
            if (String::subStringFound(token, FORMAT_BIT))
            {
                volumeData-> type = VOLUME_TYPE::BIT;
            }

            else if (String::subStringFound(token, FORMAT_8UI))
            {
                volumeData-> type = VOLUME_TYPE::UI8;
            }

            else if (String::subStringFound(token, FORMAT_16UI))
            {
                volumeData-> type = VOLUME_TYPE::UI16;
            }

            else if (String::subStringFound(token, FORMAT_32UI))
            {
                volumeData-> type = VOLUME_TYPE::UI32;
            }

            else if (String::subStringFound(token, FORMAT_64UI))
            {
                volumeData-> type = VOLUME_TYPE::UI64;
            }

            else if (String::subStringFound(token, FORMAT_F32))
            {
                volumeData-> type = VOLUME_TYPE::F32;
            }

            else if (String::subStringFound(token, FORMAT_F64))
            {
                volumeData-> type = VOLUME_TYPE::F64;
            }
            else
            {
                LOG_ERROR("Unknown volume format!");
            }
        }

        // Parse the type of the file
        if (String::subStringFound(token, std::string("sizes")))
        {
            // Get the volume dimensions
            std::vector< std::string > strings = String::split(token, ':');
            std::string dimensions = strings[1];
            strings.clear();
            strings = String::split(dimensions, 'x');
            volumeData->width = S2UI(strings[0]);
            volumeData->height= S2UI(strings[1]);
            volumeData->depth = S2UI(strings[2]);
        }

        // Check if the header is completely parsed to be able to parse the volume data
        if (String::subStringFound(token, std::string("HEADER_DONE")))
        {
            isHeaderDone = true;
        }

        // Read the binary chunk
        if (isHeaderDone)
        {
            // Parse till the end of the stream
            std::getline(dataStream, token);

            // Update the stream
            volumeData->data= token;

            // The string stream is ready
            break;
        }
    }

    // Return the volume data
    return volumeData;
}

}
