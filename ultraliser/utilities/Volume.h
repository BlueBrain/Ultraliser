#pragma once

#include <common/Common.h>
#include <data/volumes/volumes/VolumeType.hh>

namespace Ultraliser
{

// Forward declarations
class BitVolumeGrid;
template <class T> class UnsignedVolumeGrid;
template <class T> class FloatVolumeGrid;

namespace Utils
{
namespace Volume
{

/**
 * @brief writeBitGridToNRRDFile
 * @param prefix
 * @param grid
 */
void writeBitGridToNRRDFile(const std::string &prefix, const BitVolumeGrid *grid);

/**
 * @brief writeUnsignedGridToNRRDFile
 * @param prefix
 * @param grid
 */
template< class T >
void writeUnsignedGridToNRRDFile(const std::string &prefix, const UnsignedVolumeGrid<T>* grid);

/**
 * @brief writeFloatGridToNRRDFile
 * @param prefix
 * @param grid
 */
template< class T >
void writeFloatGridToNRRDFile(const std::string &prefix, const FloatVolumeGrid<T>* grid);

/**
 * @brief writeBitGridToVOLFile
 * @param prefix
 * @param grid
 * @param oneBitPerVoxel
 */
void writeBitGridToVOLFile(const std::string &prefix, const BitVolumeGrid* grid,
                           const bool& oneBitPerVoxel = true);

/**
 * @brief writeUnsignedGridToVOLFile
 * @param prefix
 * @param grid
 */
template< class T >
void writeUnsignedGridToVOLFile(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid);

/**
 * @brief writeFloatGridToVOLFile
 * @param prefix
 * @param grid
 */
template< class T >
void writeFloatGridToVOLFile(const std::string &prefix,  const FloatVolumeGrid<T>* grid);

/**
 * @brief writeHeaderFile
 * @param prefix
 * @param type
 * @param width
 * @param height
 * @param depth
 */
void writeHeaderFile(const std::string &prefix, const VOLUME_TYPE& type,
                     const size_t& width, const size_t& height, const size_t& depth);

/**
 * @brief writeBitGridToRAWFile
 * @param prefix
 * @param grid
 */
void writeBitGridToRAWFile(const std::string &prefix, const BitVolumeGrid* grid);

/**
 * @brief writeUnsignedGridToRawFile
 * @param prefix
 * @param grid
 */
template< class T >
void writeUnsignedGridToRawFile(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid);

}
}
}
