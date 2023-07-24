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

#include "BitVolumeGrid.h"
#include <algorithms/floodfill/FloodFiller.h>
#include <common/Common.h>
#include <data/common/ColorMap.h>
#include <data/volumes/voxels/DMCVoxel.h>
#include <data/volumes/grids/Projection.h>
#include <utilities/Utilities.h>
#include <data/volumes/utilities/VolumeWriter.h>

namespace Ultraliser
{

BitVolumeGrid::BitVolumeGrid(const size_t &width,
                             const size_t &height,
                             const size_t &depth,
                             BitArray* data)
    : VolumeGrid(width, height, depth)
{
    _data = data;
}

BitVolumeGrid::BitVolumeGrid(const int64_t &width, const int64_t &height, const int64_t &depth,
                             const bool &preAllocateMemory)
    : VolumeGrid(width, height, depth)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

BitVolumeGrid::BitVolumeGrid(const BitVolumeGrid *inputGrid) : VolumeGrid(*inputGrid)
{
    // Allocate the memory to be able to copy the data
    _allocateMemory();
    *_data = *inputGrid->getGridData();
}

size_t BitVolumeGrid::getNumberBytes() const
{
    return _data->getNumberBytes();
}

uint8_t BitVolumeGrid::getValueUI8(const size_t &index) const
{
    if (_data->bit(index))
        return uint8_t(1);
    return 0;
}

uint16_t BitVolumeGrid::getValueUI16(const size_t &index) const
{
    if (_data->bit(index))
        return uint16_t(1);
    return 0;
}

uint32_t BitVolumeGrid::getValueUI32(const size_t &index) const
{
    if (_data->bit(index))
        return uint32_t(1);
    return 0;
}

uint64_t BitVolumeGrid::getValueUI64(const size_t &index) const
{
    if (_data->bit(index))
        return uint64_t(255);
    return 0;
}

float BitVolumeGrid::getValueF32(const size_t &index) const
{
    if (_data->bit(index))
        return 1.f;
    return 0.f;
}

double BitVolumeGrid::getValueF64(const size_t &index) const
{
    if (_data->bit(index))
        return 1.0;
    return 0.0;
}

uint8_t BitVolumeGrid::getByte(size_t index) const
{
    return _data->getByte(index);
}

void BitVolumeGrid::addByte(const size_t &index, const uint8_t &byte)
{
    _data->addByte(index, byte);
}

void BitVolumeGrid::clear()
{
    _data->clearAll();
}

void BitVolumeGrid::fillVoxel(const size_t &index)
{
    _data->setBit(index);
}

void BitVolumeGrid::clearVoxel(const size_t &index)
{
    _data->clearBit(index);
}

bool BitVolumeGrid::isFilled(const size_t &index) const
{
    return _data->bit(index);
}

bool BitVolumeGrid::isEmpty(const size_t &index) const
{
    return !_data->bit(index);
}

void BitVolumeGrid::andWithAnotherGrid(VolumeGrid *anotherGrid)
{
    _data->operator&=(*(static_cast< BitVolumeGrid* >(anotherGrid)->getGridData()));
}

void BitVolumeGrid::orWithAnotherGrid(VolumeGrid *anotherGrid)
{
    _data->operator|=(*(static_cast< BitVolumeGrid* >(anotherGrid)->getGridData()));
}

void BitVolumeGrid::writeBitVolume(const std::string &prefix) const
{
    // Add the "-bit" suffix to make it understandable that this volume is a bit volume.
    std::string binaryPrefix = prefix + "-bit";
    writeVOL(binaryPrefix, this, true);
}

void BitVolumeGrid::writeUnsignedVolume(const std::string &prefix) const
{
    writeVOL(prefix, this);
}

void BitVolumeGrid::writeFloatVolume(const std::string &prefix) const
{
   LOG_WARNING("BitVolumeGrid::writeFloatVolume, BitVolumeGrids CANNOT be written to float volumes!");
}

void BitVolumeGrid::writeNRRDVolume(const std::string &prefix) const
{
    writeNRRD(prefix, this);
}

void BitVolumeGrid::writeRAWVolume(const std::string &prefix) const
{
    writeRAW(prefix, this);
}

void BitVolumeGrid::_allocateMemory()
{
    // Allocate the BIT ARRAY
    _data = new BitArray(I2UI64(_numberVoxels));

    // Initialize to ZERO
    _data->clearAll();
}

void BitVolumeGrid::_freeMemory()
{
    // Free the bit array
    delete _data;
}

BitVolumeGrid::~BitVolumeGrid()
{
    // Free the allocated memory
    _freeMemory();
}

}
