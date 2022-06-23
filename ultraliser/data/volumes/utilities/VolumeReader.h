#pragma once

#include <common/Headers.hh>
#include <data/volumes/utilities/VolumeData.hh>

namespace Ultraliser
{

/**
 * @brief readHeaderFile
 * @param filePath
 * @return
 */
VolumeData* readHeaderFile(const std::string& filePath);

/**
 * @brief readRawFile
 * @param filePath
 * @return
 */
std::string readRawFile(const std::string& filePath);

/**
 * @brief readRawFileToByteVector
 * @param filePath
 * @return
 */
std::vector< uint8_t > readRawFileToByteVector(const std::string& filePath);

/**
 * @brief readNRRDVolumeFile
 * @param filePath
 * @return
 */
NRRDVolumeData* readNRRDVolumeFile(const std::string& filePath);

/**
 * @brief UltraliserVolumeData
 * Reads
 * @param filePath
 * @return
 */
UltraliserVolumeData* readUltraliserVolumeFile(const std::string &filePath);


}
