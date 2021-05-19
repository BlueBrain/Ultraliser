/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_UTILITIES_FILE_H
#define ULTRALISER_UTILITIES_FILE_H

#include <common/Common.h>
#include <math/Vector3f.h>
#include <data/NeuronData.h>

namespace Ultraliser
{
namespace File
{

/**
 * @brief parseColorMap
 *
 * @param filePath
 * @param numberTags
 * @return
 */
std::vector < Vector4f > parseColorMap(const std::string &filePath, const uint64_t numberTags);

/**
 * @brief parseSomaPositions
 *
 * @param filePath
 * @return
 */
std::vector< Vector3f > parseSomaPositions(const std::string &filePath);
/**
 * @brief parseListsFile
 * @param filePath
 * @return
 */
std::vector< std::string > parseListsFile(const std::string &filePath);

/**
 * @brief parseCircuitConfiguration
 *
 * @param circuitFile
 * @return
 */
Neurons parseCircuitConfiguration(const std::string &circuitFile);

/**
 * @brief parseBoundsFile
 *
 * @param boundsFile
 * @param pMin
 * @param pMax
 */
void parseBoundsFile(std::string boundsFile, Vector3f& pMin, Vector3f& pMax);

/**
 * @brief parseVolumeList
 *
 * @param filePath
 * @return
 */
std::vector< std::string > parseVolumeList(const std::string &filePath);

/**
 * @brief parsePBRTConfig
 * Parses a given PBRT configuration.
 *
 * @param filePath
 * The file path.
 * @return
 * A vector of strings including all the data.
 */
std::vector< std::string > parsePBRTConfig(const std::string &filePath);

/**
 * @brief fileExists
 * Checks if a given file exists or not.
 *
 * @param path
 * File path.
 * @return
 * True or False.
 */
bool exists(const std::string &path);

/**
 * @brief getFileName
 * Gets the name of the file from its full path.
 *
 * @param filePath
 * The path of the file.
 * @return
 * The file name.
 */
std::string getName(std::string &filePath, bool withExtension = false);

/**
 * @brief writeFloatDistributionToFile
 *
 * @param filePath
 */
void writeFloatDistributionToFile(const std::string &filePath, std::vector<float> distribution);

/**
 * @brief writeIntegerDistributionToFile
 *
 * @param filePath
 * @param distribution
 */
void writeIntegerDistributionToFile(const std::string &filePath, std::vector<uint64_t> distribution);

/**
 * @brief getNumberLinesInFile
 *
 * @param fileName
 * @return
 */
uint64_t getNumberLinesInFile(const std::string & fileName);

/**
 * @brief readLineFromFile
 * Read one line (max 1024 chars) and exit if EOF
 *
 * @param in
 * @param exitWithEOF
 * @return
 */
char *readLineFromFile(FILE *file, bool exitWithEOF = 1);

/**
 * @brief skipCommentAndBlankLines
 * @param filePointer
 */
void skipCommentAndBlankLines(FILE *filePointer);

}
}

#endif // ULTRALISER_UTILITIES_FILE_H
