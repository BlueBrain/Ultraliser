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
#include <algorithms/FloodFiller.h>
#include <common/Common.h>
#include <data/common/ColorMap.h>
#include <data/volumes/voxels/DMCVoxel.h>
#include <data/volumes/grids/Projection.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

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

BitVolumeGrid::BitVolumeGrid(const Vec3ui_64 &dimensions, const bool &preAllocateMemory)
    : VolumeGrid(dimensions)
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



struct SomeStruct{
    char format[3];
    uint64_t width, height, depth;
};

void BitVolumeGrid::readUVOLBData(const std::string &filePath)
{
//    FILE * pFile = std::fopen(filePath.c_str(), "rb" );
//    if (pFile == NULL)
//    {
//        LOG_ERROR("Could not open the volume file [ %s ]!", filePath.c_str());
//    }
//    // Read the volume file from the input stream
//    std::ifstream imgFileStream;
//    imgFileStream.open(filePath.c_str(), std::ios::in | std::ios::binary);
//    if (imgFileStream.fail())
//    {
//        LOG_ERROR("Could not open the volume file [ %s ]!", filePath.c_str());
//    }
}


void BitVolumeGrid::readUVOLData(const std::string &filePath)
{

}

void BitVolumeGrid::loadBinaryVolumeData(const std::string &prefix)
{
    // Read the volume file from the input stream
    std::string filePath = prefix + BINARY_EXTENSION;
    std::ifstream imgFileStream;
    imgFileStream.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (imgFileStream.fail())
    {
        LOG_ERROR("Could not open the volume file [ %s ]!", filePath.c_str());
    }

    imgFileStream.read((char*) _data, I2I64(_data->getNumberBytes()));

    // Close the stream
    imgFileStream.close();
}

void BitVolumeGrid::loadUnsignedVolumeData(const std::string &rawvolumepath)
{
    // Read the volume file from the input stream
    std::string filePath = rawvolumepath + BINARY_EXTENSION;
    std::ifstream imgFileStream;
    imgFileStream.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (imgFileStream.fail())
    {
        LOG_ERROR("Could not open the volume file [ %s ]!", filePath.c_str());
    }

    // Read the file in a temporary array
    uint8_t* fileData = new uint8_t[I2UI64(_numberVoxels)];
    imgFileStream.read((char*) fileData, I2I64(_data->getNumberBytes()));

    // Close the stream
    imgFileStream.close();

    // Fill the actual data array
    for (uint64_t i = 0; i < I2UI64(_numberVoxels); ++i)
    {
        if (fileData[i] > 0)
            _data->setBit(i);
        else
            _data->clearBit(i);
    }

    // Release the data loaded from file
    delete [] fileData;
}

uint64_t BitVolumeGrid::getNumberBytes() const
{
    return _data->getNumberBytes();
}

uint8_t BitVolumeGrid::getValueUI8(const uint64_t &index) const
{
    if (_data->bit(index))
        return uint8_t(1);
    return 0;
}

uint16_t BitVolumeGrid::getValueUI16(const uint64_t &index) const
{
    if (_data->bit(index))
        return uint16_t(1);
    return 0;
}

uint32_t BitVolumeGrid::getValueUI32(const uint64_t &index) const
{
    if (_data->bit(index))
        return uint32_t(1);
    return 0;
}

uint64_t BitVolumeGrid::getValueUI64(const uint64_t &index) const
{
    if (_data->bit(index))
        return uint64_t(255);
    return 0;
}

float BitVolumeGrid::getValueF32(const uint64_t &index) const
{
    if (_data->bit(index))
        return 1.f;
    return 0.f;
}

double BitVolumeGrid::getValueF64(const uint64_t &index) const
{
    if (_data->bit(index))
        return 1.0;
    return 0.0;
}

//uint8_t BitVolumeGrid::getValue(const uint64_t &index) const
//{
//    if (_data->bit(index))
//    {
//        return uint8_t(1);
//    }
//    else
//    {
//        return uint8_t(0);
//    }
//}

uint8_t BitVolumeGrid::getByte(uint64_t index) const
{
    return _data->getByte(index);
}

void BitVolumeGrid::addByte(const uint64_t &index, const uint8_t &byte)
{
    _data->addByte(index, byte);
}

void BitVolumeGrid::clear()
{
    _data->clearAll();
}

void BitVolumeGrid::fillVoxel(const uint64_t &index)
{
    _data->setBit(index);
}

void BitVolumeGrid::clearVoxel(const uint64_t &index)
{
    _data->clearBit(index);
}

bool BitVolumeGrid::isFilled(const uint64_t &index) const
{
    return _data->bit(index);
}

bool BitVolumeGrid::isEmpty(const uint64_t &index) const
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

void BitVolumeGrid::_writeHeader(const std::string &prefix)
{
    // Header path
    std::string fileName = prefix + std::string(HEADER_EXTENSION);

    // Open the file
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    LOG_STATUS("Writing Header [ %s ]", fileName.c_str());

    // Write the dimensions
    header << getWidth() << " " << getHeight() << " " << getDepth() << std::endl;

    // Close the file
    header.close();
}

void BitVolumeGrid::writeBIN(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    std::string fileName = prefix + std::string(BINARY_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Bit per voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value = 0;
        for (uint64_t i = 0; i < 8; ++i)
        {
            if (_data->bit(voxel + i))
                value |= 1 << i;
        }

        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void BitVolumeGrid::writeRAW(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Write the image file
    std::string fileName = prefix + std::string(RAW_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value;
        if (_data->bit(voxel))
            value = FILLED_VOXEL_VALUE;
        else
            value = EMPTY_VOXEL_VALUE;

        image << value;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void BitVolumeGrid::writeNRRD(const std::string &prefix)
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

    LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value;
        if (_data->bit(voxel))
            value = FILLED_VOXEL_VALUE;
        else
            value = EMPTY_VOXEL_VALUE;

        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void BitVolumeGrid::writeUltraliserBinaryVolume(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_BIN_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "format:1bit\n");
    fprintf(fptr,"sizes:%" PRId64 "x%" PRId64 "x%" PRId64 "\n",
            getWidth(), getHeight(), getDepth());
    fprintf(fptr, "HEADER_DONE\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Bit per Voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value = 0;
        for (uint64_t i = 0; i < 8; ++i)
        {
            if (_data->bit(voxel + i))
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

void BitVolumeGrid::writeUltraliserRawVolume(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "format:1bit\n");
    fprintf(fptr,"sizes:%" PRId64 "x%" PRId64 "x%" PRId64 "\n",
            getWidth(), getHeight(), getDepth());
    fprintf(fptr, "HEADER_DONE\n");

    LOOP_STARTS("Writing Voxels (1 Byte per voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        uint8_t value;
        if (_data->bit(voxel))
            value = FILLED_VOXEL_VALUE;
        else
            value = EMPTY_VOXEL_VALUE;

        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void BitVolumeGrid::writeUltraliserFloatVolume(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    std::string fileName = prefix + std::string(ULTRALISER_VOLUME_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    // Header
    /// NOTE: The header specifies a single bit per voxel and volume dimensions
    fprintf(fptr, "f32\n");
    fprintf(fptr,"sizes: %" PRId64 " %" PRId64 " %" PRId64 "\n",
            getWidth(), getHeight(), getDepth());

    LOOP_STARTS("Writing Voxels (1 Bit per Voxel)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        float value;
        if (_data->bit(voxel))
            value = FULLED_FLOAT_VOXEL_VALUE;
        else
            value = EMPTY_FLOAT_VOXEL_VALUE;

        fputc(value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
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
