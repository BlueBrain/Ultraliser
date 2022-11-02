/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
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

#include <utilities/File.h>
#include <utilities/Parsers.h>
#include <utilities/String.h>
#include <common/Common.h>

namespace Ultraliser
{
namespace File
{

std::vector< Vector4f >parseColorMap(const std::string &filePath, const size_t numberTags)
{
    std::vector< Vector4f > colormap;
    colormap.resize(numberTags);

    std::ifstream stream;
    stream.open(filePath.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the colormap file");
    }

    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
    {
        // Split the line into 5 float values split by a space
        std::vector< std::string > tokens;
        std::istringstream iss(line);
        copy(std::istream_iterator< std::string >(iss),
              std::istream_iterator< std::string >(),
              std::back_inserter(tokens));

        size_t index = I2UI64(std::atoi(tokens[0].c_str()));
        Vector4f color;
        color.x() = D2F(std::atof(tokens[1].c_str()) / 256.0);
        color.y() = D2F(std::atof(tokens[2].c_str()) / 256.0);
        color.z() = D2F(std::atof(tokens[3].c_str()) / 256.0);
        color.w() = D2F(std::atof(tokens[4].c_str()) / 256.0);
        colormap[ index - 1 ] = color;
    }

    return colormap;
}

Neurons parseCircuitConfiguration(const std::string &circuitFile)
{
    std::ifstream stream;
    stream.open(circuitFile.c_str());
    if (!stream.good())
        LOG_ERROR("Invalid circuit file [ %s ]", circuitFile.c_str());

    // A list the contains all the neurons and their data
    Neurons neurons;

    // Parse the configuration file line by line into a list
    std::vector< std::string > configurationData;
    std::string line;
    while (std::getline(stream, line))
    {
        // Add the line to the list
        configurationData.push_back(line);
    }

    // Close the file
    stream.close();

    // Parse the neuron data from the configurationData list
    uint32_t index = 0;
    while (true)
    {
        // Parsing done
        if (index > configurationData.size() - 1)
            break;

        // Get a new line
        std::string line = configurationData[index];

        // Neuron found
        if (String::subStringFound(line, "NEURON"))
        {
            // Construct a new neuron object
            Neuron neuron;

            while (true)
            {
                // Get a new line
                index++;
                line = configurationData[index];

                if (String::subStringFound(line, GID_KEY))
                {
                    neuron.gid = Parsers::getStringValue(line);
                }
                if (String::subStringFound(line, TAG_KEY))
                {
                    neuron.tag = I2UI8(Parsers::getIntegerValue(line));
                }
                if (String::subStringFound(line, LAYER_KEY))
                {
                    neuron.layer = Parsers::getStringValue(line);
                }
                if (String::subStringFound(line, COLUMN_KEY))
                {
                    neuron.column = Parsers::getStringValue(line);
                }
                else if (String::subStringFound(line, MTYPE_KEY))
                {
                    neuron.morphologyType = Parsers::getStringValue(line);
                }
                else if (String::subStringFound(line, MLABEL_KEY))
                {
                    neuron.morphologyLabel = Parsers::getStringValue(line);
                }
                else if (String::subStringFound(line, POSITION_KEY))
                {
                    // Tokenize the line
                    std::vector< std::string > tokens;
                    std::istringstream iss(line);
                    copy(std::istream_iterator< std::string >(iss),
                          std::istream_iterator< std::string >(),
                          std::back_inserter(tokens));

                    Vector3f somaPosition;
                    somaPosition.x() = D2F(std::atof(tokens[1].c_str()));
                    somaPosition.y() = D2F(std::atof(tokens[2].c_str()));
                    somaPosition.z() = D2F(std::atof(tokens[3].c_str()));
                    neuron.somaPosition = somaPosition;
                }
                else if (String::subStringFound(line, ORIENTATION_KEY))
                {
                    // Tokenize the line
                    std::vector< std::string > tokens;
                    std::istringstream iss(line);
                    copy(std::istream_iterator< std::string >(iss),
                          std::istream_iterator< std::string >(),
                          std::back_inserter(tokens));
                    neuron.yOrientation = D2F(std::atof(tokens[1].c_str()));
                }
                else if (String::subStringFound(line, TRANSFORM_KEY))
                {
                    // Tokenize the line
                    std::vector< std::string > tokens;
                    std::istringstream iss(line);
                    copy(std::istream_iterator< std::string >(iss),
                          std::istream_iterator< std::string >(),
                          std::back_inserter(tokens));

                    // Parse the transformations
                    float t00 = D2F(std::atof(tokens[1].c_str()));
                    float t01 = D2F(std::atof(tokens[2].c_str()));
                    float t02 = D2F(std::atof(tokens[3].c_str()));
                    float t03 = D2F(std::atof(tokens[4].c_str()));
                    float t10 = D2F(std::atof(tokens[5].c_str()));
                    float t11 = D2F(std::atof(tokens[6].c_str()));
                    float t12 = D2F(std::atof(tokens[7].c_str()));
                    float t13 = D2F(std::atof(tokens[8].c_str()));
                    float t20 = D2F(std::atof(tokens[9].c_str()));
                    float t21 = D2F(std::atof(tokens[10].c_str()));
                    float t22 = D2F(std::atof(tokens[11].c_str()));
                    float t23 = D2F(std::atof(tokens[12].c_str()));
                    float t30 = D2F(std::atof(tokens[13].c_str()));
                    float t31 = D2F(std::atof(tokens[14].c_str()));
                    float t32 = D2F(std::atof(tokens[15].c_str()));
                    float t33 = D2F(std::atof(tokens[16].c_str()));
                    neuron.localToGlobalTransform =
                            Matrix4f(t00, t01, t02, t03,
                                     t10, t11, t12, t13,
                                     t20, t21, t22, t23,
                                     t30, t31, t32, t33);
                }
                else
                {
                    // Split the line into 5 values split by a space
                    std::vector< std::string > tokens;
                    std::istringstream iss(line);
                    copy(std::istream_iterator< std::string >(iss),
                          std::istream_iterator< std::string >(),
                          std::back_inserter(tokens));

                    // New line
                    if (!(tokens.size() > 0))
                        break;

                    // Otherwise, it is an UNRECOGNIZED parameter
                    else
                        continue;
                }
            }

            // Add the neuron to the list
            neurons.push_back(neuron);

            // Increment the index for the next neuron
            index++;
        }
        else
            break;
    }

   // Return a list of neurons
    return neurons;
}

void parseBoundsFile(std::string boundsFile, Vector3f& pMin, Vector3f& pMax)
{
    std::ifstream stream;
    stream.open(boundsFile.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the bounds file");
    }

    std::string line;
    std::getline(stream, line);

    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    pMin.x() = S2F(tokens[0]);
    pMin.y() = S2F(tokens[1]);
    pMin.z() = S2F(tokens[2]);

    pMax.x() = S2F(tokens[3]);
    pMax.y() = S2F(tokens[4]);
    pMax.z() = S2F(tokens[5]);
}

std::vector< size_t > parseIsovaluesFile(const std::string &filePath)
{
    std::ifstream stream;
    stream.open(filePath.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot parse the isovalues file [%s]", filePath.c_str());
    }

    std::string line;
    std::getline(stream, line);

    // Split the line into 5 values split by a space
    std::vector< std::string > tokens;
    std::istringstream iss(line);
    copy(std::istream_iterator< std::string >(iss),
          std::istream_iterator< std::string >(),
          std::back_inserter(tokens));

    std::vector<size_t> isovalues;
    for (size_t i = 0; i < tokens.size(); ++i)
    {
        isovalues.push_back(atoll(tokens[i].c_str()));
    }

    return isovalues;
}


std::vector< std::string > parseVolumeList(const std::string &filePath)
{
    std::vector< std::string > list;

    std::ifstream stream;
    stream.open(filePath.c_str());
    if (!stream.good())
        LOG_ERROR("Cannot parse the volume file [ %s ]", filePath.c_str());


    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
            list.push_back(std::string(line));

    return list;
}

std::vector< std::string > parsePBRTConfig(const std::string &filePath)
{
    std::vector< std::string > config;

    std::ifstream stream;
    stream.open(filePath.c_str());
    if (!stream.good())
        LOG_ERROR("Cannot parse the config file [ %s ]", filePath.c_str());

    // Parse the list line by line
    std::string line;
    while (std::getline(stream, line))
    {
        config.push_back(std::string(line));
    }

    return config;
}

bool exists(const std::string &path)
{
    struct stat buffer;
    return (stat (path.c_str(), &buffer) == 0);
}

std::string getName(std::string& filePath, bool withExtension)
{
    struct tokensBS: std::ctype<char>
    {
        tokensBS(): std::ctype<char>(get_table()) {}

        static std::ctype_base::mask const* get_table()
        {
            typedef std::ctype<char> cctype;
            static const cctype::mask *const_rc= cctype::classic_table();

            static cctype::mask rc[cctype::table_size];
            std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));

            rc['/'] = std::ctype_base::space;
            return &rc[0];
        }
    };

    struct tokensDot: std::ctype<char>
    {
        tokensDot(): std::ctype<char>(get_table()) {}

        static std::ctype_base::mask const* get_table()
        {
            typedef std::ctype<char> cctype;
            static const cctype::mask *const_rc= cctype::classic_table();

            static cctype::mask rc[cctype::table_size];
            std::memcpy(rc, const_rc, cctype::table_size * sizeof(cctype::mask));

            rc['.'] = std::ctype_base::space;
            return &rc[0];
        }
    };

    // Split the file path to a std::vector with the  back slach separator
    std::stringstream streamBS(filePath);
    streamBS.imbue(std::locale(std::locale(), new tokensBS()));
    std::istream_iterator<std::string> begin(streamBS);
    std::istream_iterator<std::string> end;
    std::vector<std::string> stringsVectorBS(begin, end);

    // Now, we have the file name with the extension
    std::string fileNameWithExtension = stringsVectorBS.at(stringsVectorBS.size() - 1);

    // With extension
    if (withExtension)
        return fileNameWithExtension;

    // Split the file name to a std::vector with the . separator
    std::stringstream streamDot(fileNameWithExtension);
    streamDot.imbue(std::locale(std::locale(), new tokensDot()));
    std::istream_iterator<std::string> beginDot(streamDot);
    std::istream_iterator<std::string> endDot;
    std::vector<std::string> stringsVectorDot(beginDot, endDot);

    // Withour extension
    return stringsVectorDot.at(0);
}

void writeFloatDistributionToFile(const std::string &filePath,
                                  const std::vector< float > &distribution)
{
    std::ofstream outputStream(filePath.c_str());
    if (!outputStream.good())
        LOG_ERROR("Cannot write file [ %s ]", filePath.c_str());

    for (size_t i = 0; i < distribution.size(); ++i)
        if (distribution[i] > 0.f)
            outputStream << i << " " << distribution[i] << NEW_LINE;

    // Close the stream
    outputStream.close();
}

void writeFloatHistogramToFile(const std::string &filePath,
                               const std::vector< size_t > &histogram,
                               const std::vector< float > &bins)
{
    std::ofstream outputStream(filePath.c_str());
    if (!outputStream.good())
        LOG_ERROR("Cannot write file [ %s ]", filePath.c_str());

    // Histogram
    for (size_t i = 0; i < histogram.size(); ++i)
        outputStream << histogram[i] << " ";
    outputStream << NEW_LINE;

    // Bins
    for (size_t i = 0; i < bins.size(); ++i)
        outputStream << bins[i] << " ";
    outputStream << NEW_LINE;

    // Close the stream
    outputStream.close();
}

void writeIntegerDistributionToFile(const std::string &filePath,
                                    std::vector<size_t> distribution)
{
    std::ofstream outputStream(filePath.c_str());
    if (!outputStream.good())
        LOG_ERROR("Cannot write file [ %s ]", filePath.c_str());

    for (size_t i = 0; i < distribution.size(); ++i)
        if (distribution[i] > 0.f)
            outputStream << i << " " << distribution[i] << NEW_LINE;

    // Close the stream
    outputStream.close();
}

size_t getNumberLinesInFile(const std::string & fileName)
{
    // NUmber of line in the file
    size_t numberLines = 0;

    // Line-by-line
    std::string line;

    // File stream
    std::ifstream fileStream(fileName.c_str());

    // Do it line by line
    while (std::getline(fileStream, line))
        ++numberLines;

    // Close the file stream
    fileStream.close();

    // Return the number of lines
    return numberLines;
}


// Swap endian-ness for four-byte elements

void endianSwapLong(unsigned char *p)
{
    unsigned char b0,b1,b2,b3;

    b0 = *p;
    b1 = *(p+1);
    b2 = *(p+2);
    b3 = *(p+3);

    *p = b3;
    *(p+1) = b2;
    *(p+2) = b1;
    *(p+3) = b0;
}

char *readLineFromFile(FILE *in, bool exit_on_eof)
{
#define MAX_READLINE_CHARS	1024
    static char line[MAX_READLINE_CHARS];
    int i=0;
    char c;

    while ((c = fgetc(in)) != '\n' && i<(MAX_READLINE_CHARS-1))
        if (c==EOF)
        {
            if (exit_on_eof) LOG_ERROR("\nUnexpected end of file!\n");
            else return nullptr;
        }
        else if (c != '\r') line[++i] = c;
    line[i] = '\0';

    if (i==MAX_READLINE_CHARS-1)
        LOG_WARNING("readLineFromFile: Line is too long. Truncated !\n");

    return line;
}

void skipCommentAndBlankLines(FILE *filePointer)
{
    long pos0;
    char *line, s[2];
    do
    {
        pos0 = ftell(filePointer);

        line = readLineFromFile(filePointer);}
    while (line[0] == '#' || line[0] == '\0' || !sscanf(line,"%1s",s));

    fseek(filePointer, pos0, SEEK_SET);
}


}
}
