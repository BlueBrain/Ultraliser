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

#include <data/volumes/grids/VolumeGrid.h>
#include "VoxelGrid.h"
#include <data/volumes/voxels/DMCVoxel.h>
#include <algorithms/FloodFiller.h>
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <data/volumes/grids/Projection.h>

namespace Ultraliser
{

VoxelGrid::VoxelGrid(const int64_t &width,
                     const int64_t &height,
                     const int64_t &depth,
                     const bool& preAllocateMemory)
    : VolumeGrid(width, height, depth)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

VoxelGrid::VoxelGrid(const Vec3i_64& dimensions,
                     const bool& preAllocateMemory)
    : VolumeGrid(dimensions)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

VoxelGrid::VoxelGrid(const VoxelGrid* inputGrid) : VolumeGrid(*inputGrid)
{

    // Allocate the memory to be able to copy the data
    _allocateMemory();

    // Fill the _data array
    Voxels inputData = inputGrid->getGridData();

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int64_t i = 0; i < _numberVoxels; ++i)
    {
        _data[i] = inputData[i];
    }
}

void VoxelGrid::loadBinaryVolumeData(const std::string &prefix)
{

}

void VoxelGrid::loadByteVolumeData(const std::string &prefix)
{
    // Read the volume file from the input stream
    std::string filePath = prefix + BINARY_EXTENSION;
    std::ifstream imgFileStream;
    imgFileStream.open(filePath.c_str(), std::ios::in | std::ios::binary);
    if (imgFileStream.fail())
    {
        LOG_ERROR("Could not open the volume file %s!", filePath.c_str());
    }

    // Read the file in a temporary array
    uint8_t* fileData = new uint8_t[I2UI64(_numberVoxels)];
    imgFileStream.read((char*) fileData, _numberVoxels);

    // Close the stream
    imgFileStream.close();

    // Fill the actual data array
    for (uint64_t i = 0; i < I2UI64(_numberVoxels); ++i)
    {
        _data[i].value = fileData[i];
    }

    // Release the data loaded from file
    delete [] fileData;
}

uint64_t VoxelGrid::getNumberBytes() const
{
    return I2UI64(_numberVoxels);
}

uint8_t VoxelGrid::getValue(const uint64_t &index) const
{
    return _data[index].value;
}

uint8_t VoxelGrid::getByte(uint64_t index) const
{
    return _data[index].value;
}

void VoxelGrid::addByte(const uint64_t &index, const uint8_t &byte)
{
    _data[index].value = byte;
}

void VoxelGrid::clear()
{
    for (int64_t i = 0; i < _numberVoxels; ++i)
        _data[i].value = 0;
}

void VoxelGrid::fillVoxel(const uint64_t &index)
{
    _data[index].value = 255;
}

void VoxelGrid::clearVoxel(const uint64_t &index)
{
    _data[index].value = 0;
}

bool VoxelGrid::isFilled(const uint64_t &index) const
{
    return (_data[index].value > 0);
}

bool VoxelGrid::isEmpty(const uint64_t &index) const
{
    return (_data[index].value == 0);
}

void VoxelGrid::_writeHeader(const std::string &prefix)
{
    std::string fileName = prefix + std::string(HEADER_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);
    header << getWidth() << " " << getHeight() << " " << getDepth() << std::endl;
    header.close();
}

void VoxelGrid::writeRAW(const std::string &prefix)
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
    for (int64_t voxel = 0; voxel < _numberVoxels; voxel++)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);
        image << _data[voxel].value;
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void VoxelGrid::writeBIN(const std::string &prefix)
{
//    // Starts the timer
//    TIMER_SET;

//    // Write the header file
//    _writeHeader(prefix);

//    // Write the image file
//    std::string fileName = prefix + std::string(BINARY_EXTENSION);
//    std::fstream image;
//    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

//    LOOP_STARTS("Writing Voxels (1 Bit)");
//    for (int64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
//    {
//        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

//        uint8_t value = 0;
//        for (int64_t i = 0; i < 8; ++i)
//        {
//            if (_data->bit(I2UI64(voxel + i)))
//                value |= 1 << i;
//        }
//        image << value;
//    }
//    LOOP_DONE;

//    // Statistics
//    LOG_STATS(GET_TIME_SECONDS);

//    // Close the file
//    image.close();
}

void VoxelGrid::writeNRRD(const std::string &prefix)
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
    fprintf(fptr,"sizes %" PRId64 " %" PRId64 " %" PRId64 "\n",
            getWidth(), getHeight(), getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "encoding: raw\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels");
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        fputc(_data[voxel].value, fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void VoxelGrid::andWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    Voxels inputData = static_cast< VoxelGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel].value &= inputData[voxel].value;
    }
}

void VoxelGrid::orWithAnotherGrid(VolumeGrid *anotherGrid)
{
    // Get a reference to the data
    Voxels inputData = static_cast< VoxelGrid* >(anotherGrid)->getGridData();

    // Copy the data
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        _data[voxel].value |= inputData[voxel].value;
    }
}

void VoxelGrid::_allocateMemory()
{
    // Allocate the array
    _data.resize(I2UI64(_numberVoxels));
}

void VoxelGrid::_freeMemory()
{
    _data.clear();
}

VoxelGrid::~VoxelGrid()
{
    _freeMemory();
}

}
