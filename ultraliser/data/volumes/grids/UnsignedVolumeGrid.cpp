/***************************************************************************************************
 * Copyright (c) 2016 - 2022
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

#include <data/volumes/grids/VolumeGrid.h>
#include "UnsignedVolumeGrid.h"
#include <data/volumes/voxels/DMCVoxel.h>
#include <algorithms/FloodFiller.h>
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <data/volumes/grids/Projection.h>
#include <utilities/Volume.h>

namespace Ultraliser
{

template <class T>
UnsignedVolumeGrid<T>::UnsignedVolumeGrid(const size_t &width,
                   const size_t &height,
                   const size_t &depth,
                   T* data)
    : VolumeGrid(width, height, depth)
{
    _allocateMemory();

    // Update the data
    for (uint64_t i = 0; i < width * height * depth; ++i)
        _data[i] =  data[i];

    // _data = data;
}

template <class T>
UnsignedVolumeGrid<T>::UnsignedVolumeGrid(const uint64_t &width,
                                          const uint64_t &height,
                                          const uint64_t &depth,
                                          const bool &preAllocateMemory)
    : VolumeGrid(width, height, depth)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

template <class T>
UnsignedVolumeGrid<T>::UnsignedVolumeGrid(const uint64_t &width,
                                          const uint64_t &height,
                                          const uint64_t &depth,
                                          std::vector<T> &data)
    : VolumeGrid(width, height, depth)
{
    _allocateMemory();

    // Update the data
    for (uint64_t i = 0; i < _numberVoxels; ++i)
        _data[i] =  data[i];
}


template <class T>
UnsignedVolumeGrid<T>::UnsignedVolumeGrid(const UnsignedVolumeGrid* inputGrid)
    : VolumeGrid(*inputGrid)
{
    // Allocate the memory to be able to copy the data
    _allocateMemory();

    // Fill the _data array
    const T* inputData = inputGrid->getGridData();

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberVoxels; ++i)
    {
        _data[i] = inputData[i];
    }
}


template <class T>
uint64_t UnsignedVolumeGrid<T>::getNumberBytes() const
{
    if (typeid (T) == typeid (uint8_t))
    {
        return _numberVoxels;
    }
    else if (typeid (T) == typeid (uint16_t))
    {
        return _numberVoxels * 2;
    }
    else if (typeid (T) == typeid (uint32_t))
    {
        return _numberVoxels * 4;
    }
    else if (typeid (T) == typeid (uint64_t))
    {
        return _numberVoxels * 8;
    }
    else
    {
        LOG_ERROR("Undefined type in UnsignedVolumeGrid<T>::getNumberBytes()!");
        return 0;
    }
}

template <class T>
uint8_t UnsignedVolumeGrid<T>::getValueUI8(const uint64_t &index) const
{
    if (typeid (T) == typeid (uint8_t))
    {
        return _data[index];
    }
    else
    {
        LOG_ERROR("Precision error in UnsignedVolumeGrid<T>::getValueUI8()!");
        return 0;
    }
}

template <class T>
uint16_t UnsignedVolumeGrid<T>::getValueUI16(const uint64_t &index) const
{
    if (typeid (T) == typeid (uint8_t))
    {
        return static_cast< uint16_t >(_data[index]);
    }
    else if (typeid (T) == typeid (uint16_t))
    {
        return _data[index];
    }
    else
    {
        LOG_ERROR("Precision error in UnsignedVolumeGrid<T>::getValueUI16()!");
        return 0;
    }
}

template <class T>
uint32_t UnsignedVolumeGrid<T>::getValueUI32(const uint64_t &index) const
{
    if (typeid (T) == typeid (uint8_t) || typeid (T) == typeid (uint16_t))
    {
        return static_cast< uint32_t >(_data[index]);
    }
    else if (typeid (T) == typeid (uint32_t))
    {
        return _data[index];
    }
    else
    {
        LOG_ERROR("Precision error in UnsignedVolumeGrid<T>::getValueUI32()!");
        return 0;
    }
}

template <class T>
uint64_t UnsignedVolumeGrid<T>::getValueUI64(const uint64_t &index) const
{
    return static_cast< uint64_t >(_data[index]);
}

template <class T>
float UnsignedVolumeGrid<T>::getValueF32(const uint64_t &index) const
{
    return static_cast< float >(_data[index]);
}

template <class T>
double UnsignedVolumeGrid<T>::getValueF64(const uint64_t &index) const
{
    return static_cast< double >(_data[index]);
}

template <class T>
uint8_t UnsignedVolumeGrid<T>::getByte(uint64_t index) const
{
    if (typeid (T) == typeid (uint8_t))
    {
        return _data[index];
    }
    else
    {
        LOG_ERROR("Unimplemented function UnsignedVolumeGrid<T>::getByte()!");
        return 0;
    }
}
template <class T>

void UnsignedVolumeGrid<T>::addByte(const uint64_t &index, const uint8_t &byte)
{
    if (typeid (T) == typeid (uint8_t))
    {
        _data[index] = byte;
    }
    else
    {
        LOG_ERROR("Unimplemented function UnsignedVolumeGrid<T>::addByte()!");
    }
}

template <class T>
void UnsignedVolumeGrid<T>::clear()
{
    for (int64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = static_cast<T>(0);
}

template <class T>
void UnsignedVolumeGrid<T>::fillVoxel(const uint64_t &index)
{
    _data[index] = 255;
}

template <class T>
void UnsignedVolumeGrid<T>::clearVoxel(const uint64_t &index)
{
    _data[index] = static_cast<T>(0);
}

template <class T>
bool UnsignedVolumeGrid<T>::isFilled(const uint64_t &index) const
{
    return _data[index] > 0;
}

template <class T>
bool UnsignedVolumeGrid<T>::isEmpty(const uint64_t &index) const
{
    return _data[index] == 0;
}

template <class T>
void UnsignedVolumeGrid<T>::andWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    T* inputGridData = static_cast< UnsignedVolumeGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel] &= inputGridData[voxel];
    }
}

template <class T>
void UnsignedVolumeGrid<T>::orWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    T* inputGridData = static_cast< UnsignedVolumeGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel] |= inputGridData[voxel];
    }
}

template <class T>
void UnsignedVolumeGrid<T>::writeBitVolume(const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "1bit\n");
    fprintf(fptr,"sizes: %" PRId64 " %" PRId64 " %" PRId64 "\n",
            getWidth(), getHeight(), getDepth());

    // Create a BitArray
    auto binData = std::make_unique< BitArray >(_numberVoxels);

    // Fill the BitArray
    LOOP_STARTS("Filling the BitArray");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);
        if (_data[voxel])
        {
            binData->setBit(voxel);
        }
        else
        {
            binData->clearBit(voxel);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);


    LOOP_STARTS("Writing Voxels (1 Bit per Voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value = 0;
        for (uint64_t i = 0; i < 8; ++i)
        {
            if (binData->bit(voxel + i))
                value |= 1 << i;
        }

        // Add the byte value to the volume
        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

template <class T>
void UnsignedVolumeGrid<T>::writeUnsignedVolume(const std::string &prefix) const
{
    Utils::Volume::writeVOL(prefix, this);
}

template <class T>
void UnsignedVolumeGrid<T>::writeFloatVolume(const std::string &prefix) const
{
    LOG_WARNING("Unimplemented UnsignedVolumeGrid<T>::writeFloatVolume!");
}

template <class T>
void UnsignedVolumeGrid<T>::writeNRRDVolume(const std::string &prefix) const
{
    Utils::Volume::writeNRRD<T>(prefix, this);
}

template <class T>
void UnsignedVolumeGrid<T>::writeRAWVolume(const std::string &prefix) const
{
    Utils::Volume::writeRAW<T>(prefix, this);
}

template <class T>
void UnsignedVolumeGrid<T>::_allocateMemory()
{
    // Allocate the array
    _data = new T[_numberVoxels];
    for (uint64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = static_cast<T>(0);
}

template <class T>
void UnsignedVolumeGrid<T>::_freeMemory()
{
    delete [] _data;
}

template <class T>
UnsignedVolumeGrid<T>::~UnsignedVolumeGrid()
{
    _freeMemory();
}

template class UnsignedVolumeGrid<uint8_t>;
template class UnsignedVolumeGrid<uint16_t>;
template class UnsignedVolumeGrid<uint32_t>;
template class UnsignedVolumeGrid<uint64_t>;

}
