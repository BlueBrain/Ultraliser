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
 * @brief writeNRRD
 * @param prefix
 * @param grid
 */
void writeNRRD(const std::string &prefix, const BitVolumeGrid *grid);

/**
 * @brief writeNRRD
 * @param prefix
 * @param grid
 */
template< class T >
void writeNRRD(const std::string &prefix, const UnsignedVolumeGrid<T>* grid);

/**
 * @brief writeNRRD
 * @param prefix
 * @param grid
 */
template< class T >
void writeNRRD(const std::string &prefix, const FloatVolumeGrid<T>* grid);

/**
 * @brief writeVOL
 * @param prefix
 * @param grid
 * @param oneBitPerVoxel
 */
void writeVOL(const std::string &prefix, const BitVolumeGrid* grid,
              const bool& oneBitPerVoxel = true);

/**
 * @brief writeVOL
 * @param prefix
 * @param grid
 */
template< class T >
void writeVOL(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid);

/**
 * @brief writeVOL
 * @param prefix
 * @param grid
 */
template< class T >
void writeVOL(const std::string &prefix,  const FloatVolumeGrid<T>* grid);

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
 * @brief writeRAW
 * @param prefix
 * @param grid
 */
void writeRAW(const std::string &prefix, const BitVolumeGrid* grid);

/**
 * @brief writeRAW
 * @param prefix
 * @param grid
 */
template< class T >
void writeRAW(const std::string &prefix,  const UnsignedVolumeGrid<T>* grid);

}
}
}
