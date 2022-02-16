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

#include <data/volumes/grids/VolumeGrid.h>
#include "ByteVolumeGrid.h"
#include <data/volumes/voxels/DMCVoxel.h>
#include <algorithms/FloodFiller.h>
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <data/volumes/grids/Projection.h>

namespace Ultraliser
{

ByteVolumeGrid::ByteVolumeGrid(const int64_t &width,
                               const int64_t &height,
                               const int64_t &depth,
                               const bool &preAllocateMemory)
    : VolumeGrid(width, height, depth)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

ByteVolumeGrid::ByteVolumeGrid(const Vec3i_64& dimensions,
                               const bool &preAllocateMemory)
    : VolumeGrid(dimensions)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

ByteVolumeGrid::ByteVolumeGrid(const ByteVolumeGrid* inputGrid) : VolumeGrid(*inputGrid)
{
    // Allocate the memory to be able to copy the data
    _allocateMemory();

    // Fill the _data array
    uint8_t* inputData = inputGrid->getGridData();

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int64_t i = 0; i < _numberVoxels; ++i)
    {
        _data[i] = inputData[i];
    }
}

void ByteVolumeGrid::loadBinaryVolumeData(const std::string &prefix)
{

}

void ByteVolumeGrid::loadByteVolumeData(const std::string &prefix)
{
    // Read the volume file from the input stream
    std::string filePath = prefix + RAW_EXTENSION;
    std::ifstream imgFileStream;
    imgFileStream.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (imgFileStream.fail())
    {
        LOG_ERROR("Could not open the volume file %s!", filePath.c_str());
    }

    imgFileStream.read((char*) _data, _numberVoxels);

    // Close the stream
    imgFileStream.close();
}

uint64_t ByteVolumeGrid::getNumberBytes() const
{
    return I2UI64(_numberVoxels);
}

uint8_t ByteVolumeGrid::getValue(const uint64_t &index) const
{
    return _data[index];
}

uint8_t ByteVolumeGrid::getByte(uint64_t index) const
{
    return _data[index];
}

void ByteVolumeGrid::addByte(const uint64_t &index, const uint8_t &byte)
{
    _data[index] = byte;
}

void ByteVolumeGrid::clear()
{
    for (int64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = 0;
}

void ByteVolumeGrid::fillVoxel(const uint64_t &index)
{
    _data[index] = 255;
}

void ByteVolumeGrid::clearVoxel(const uint64_t &index)
{
    _data[index] = 0;
}

bool ByteVolumeGrid::isFilled(const uint64_t &index) const
{
    return _data[index] > 0;
}

bool ByteVolumeGrid::isEmpty(const uint64_t &index) const
{
    return _data[index] == 0;
}

void ByteVolumeGrid::andWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    uint8_t* inputGridData = static_cast< ByteVolumeGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel] &= inputGridData[voxel];
    }
}

void ByteVolumeGrid::orWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    uint8_t* inputGridData = static_cast< ByteVolumeGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel] |= inputGridData[voxel];
    }
}

void ByteVolumeGrid::_writeHeader(const std::string &prefix)
{
    std::string fileName = prefix + std::string(HEADER_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    header << getWidth() << " " << getHeight() << " " << getDepth() << std::endl;
    header.close();
}

void ByteVolumeGrid::writeRAW(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Write the image file
    std::string fileName = prefix + std::string(RAW_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

    LOOP_STARTS("Writing Voxels (1 Byte)");
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);
        image << _data[voxel];
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void ByteVolumeGrid::writeNRRD(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // File name
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    fprintf(fptr, "type: unsigned char\n");
    fprintf(fptr, "dimension: 3\n");
    fprintf(fptr,"sizes: %" PRId64 " %" PRId64 " %" PRId64 "\n",
            getWidth(), getHeight(), getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "encoding: raw\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels");
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        fputc(_data[voxel], fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void ByteVolumeGrid::writeBIN(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Create a BitArray
    auto binData = std::make_unique< BitArray >(I2UI64(_numberVoxels));

    // Fill the BitArray
    LOOP_STARTS("Filling the BitArray");
    for (int64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
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

    // Write the image file
    TIMER_RESET;
    std::string fileName = prefix + std::string(BINARY_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

    LOOP_STARTS("Writing Voxels (1 Bit)");
    for (int64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value = 0;
        for (int64_t i = 0; i < 8; ++i)
        {
            if (binData->bit(I2UI64(voxel + i)))
                value |= 1 << i;
        }
        image << value;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void ByteVolumeGrid::_allocateMemory()
{
    // Allocate the array
    _data = new uint8_t[I2UI64(_numberVoxels)];
    for (int64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = 0;
}

void ByteVolumeGrid::_freeMemory()
{
    delete [] _data;
}

ByteVolumeGrid::~ByteVolumeGrid()
{
    _freeMemory();
}

}
