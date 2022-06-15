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
#include "FloatVolumeGrid.h"
#include <data/volumes/voxels/DMCVoxel.h>
#include <algorithms/FloodFiller.h>
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <data/volumes/grids/Projection.h>

namespace Ultraliser
{

template <class T>
FloatVolumeGrid<T>::FloatVolumeGrid(const uint64_t &width,
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
FloatVolumeGrid<T>::FloatVolumeGrid(const Vec3ui_64& dimensions,
                                  const bool &preAllocateMemory)
    : VolumeGrid(dimensions)
{
    // Allocate the memory
    if (preAllocateMemory)
    {
        _allocateMemory();
    }
}

template <class T>
FloatVolumeGrid<T>::FloatVolumeGrid(const FloatVolumeGrid* inputGrid) : VolumeGrid(*inputGrid)
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
void FloatVolumeGrid<T>::loadBinaryVolumeData(const std::string &prefix)
{

}

template <class T>
void FloatVolumeGrid<T>::loadByteVolumeData(const std::string &prefix)
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

template <class T>
uint64_t FloatVolumeGrid<T>::getNumberBytes() const
{
    if (typeid (T) == typeid (uint8_t))
        return _numberVoxels;
    else if (typeid (T) == typeid (uint16_t))
        return _numberVoxels * 2;
    else if (typeid (T) == typeid (uint32_t))
        return _numberVoxels * 4;
    else if (typeid (T) == typeid (uint64_t))
        return _numberVoxels * 8;
    else if (typeid (T) == typeid (float))
        return _numberVoxels * 4;
    else if (typeid (T) == typeid (double))
        return _numberVoxels * 8;
    else
        LOG_ERROR("Undefined type!");
        return 0;
}

template <class T>
uint8_t FloatVolumeGrid<T>::getValueUI8(const uint64_t &index) const
{
    LOG_ERROR("FloatVolumeGrid<T>::getValueUI8 Unimplemented!");
}

template <class T>
uint16_t FloatVolumeGrid<T>::getValueUI16(const uint64_t &index) const
{
    LOG_ERROR("FloatVolumeGrid<T>::getValueUI16 Unimplemented!");
}

template <class T>
uint32_t FloatVolumeGrid<T>::getValueUI32(const uint64_t &index) const
{
    LOG_ERROR("FloatVolumeGrid<T>::getValueUI32 Unimplemented!");
}

template <class T>
uint64_t FloatVolumeGrid<T>::getValueUI64(const uint64_t &index) const
{
    LOG_ERROR("FloatVolumeGrid<T>::getValueUI64 Unimplemented!");
}

template <class T>
float FloatVolumeGrid<T>::getValueF32(const uint64_t &index) const
{
    return static_cast< float >(_data[index]);
}

template <class T>
double FloatVolumeGrid<T>::getValueF64(const uint64_t &index) const
{
    return static_cast< double >(_data[index]);
}

template <class T>
uint8_t FloatVolumeGrid<T>::getValue(const uint64_t &index) const
{
    return _data[index];
}

template <class T>
uint8_t FloatVolumeGrid<T>::getByte(uint64_t index) const
{
    return _data[index];
}
template <class T>

void FloatVolumeGrid<T>::addByte(const uint64_t &index, const uint8_t &byte)
{
    _data[index] = byte;
}
template <class T>

void FloatVolumeGrid<T>::clear()
{
    for (int64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = 0.f;
}

template <class T>
void FloatVolumeGrid<T>::fillVoxel(const uint64_t &index)
{
    _data[index] = 1.f;
}

template <class T>
void FloatVolumeGrid<T>::clearVoxel(const uint64_t &index)
{
    _data[index] = 0.f;
}

template <class T>
bool FloatVolumeGrid<T>::isFilled(const uint64_t &index) const
{
    return _data[index] > 0.0;
}

template <class T>
bool FloatVolumeGrid<T>::isEmpty(const uint64_t &index) const
{
    return _data[index] == 0;
}

template <class T>
void FloatVolumeGrid<T>::andWithAnotherGrid(VolumeGrid *anotherGrid)
{
    LOG_ERROR("FloatVolumeGrid<T>::andWithAnotherGrid Unimplemented!");
}

template <class T>
void FloatVolumeGrid<T>::orWithAnotherGrid(VolumeGrid *anotherGrid)
{
    LOG_ERROR("FloatVolumeGrid<T>::orWithAnotherGrid Unimplemented!");
}

template <class T>
void FloatVolumeGrid<T>::_writeHeader(const std::string &prefix)
{
    std::string fileName = prefix + std::string(HEADER_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    // Volume format
    if (typeid (T) == typeid (float))
        header << "f32" << std::endl;
    else if (typeid (T) == typeid (double))
        header << "f64" << std::endl;

    // Volume dimensions
    header << getWidth() << " " << getHeight() << " " << getDepth() << std::endl;
    header.close();
}

template <class T>
void FloatVolumeGrid<T>::writeRAW(const std::string &prefix)
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

template <class T>
void FloatVolumeGrid<T>::writeNRRD(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // File name
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    if (typeid (T) == typeid (uint8_t))
        fprintf(fptr, "type: unsigned char\n");
    else if (typeid (T) == typeid (uint16_t))
        fprintf(fptr, "type: unsigned short\n");
    else if (typeid (T) == typeid (uint32_t))
        fprintf(fptr, "type: unsigned int\n");
    else if (typeid (T) == typeid (uint64_t))
        fprintf(fptr, "type: unsigned long\n");
    else if (typeid (T) == typeid (float))
        fprintf(fptr, "type: float\n");
    else if (typeid (T) == typeid (double))
        fprintf(fptr, "type: double\n");
    else
        LOG_ERROR("Undefined volume type!");

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

template <class T>
void FloatVolumeGrid<T>::writeBIN(const std::string &prefix)
{
    // Starts the timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

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

template <class T>
void FloatVolumeGrid<T>::_allocateMemory()
{
    // Allocate the array
    _data = new T[_numberVoxels];
    for (uint64_t i = 0; i < _numberVoxels; ++i)
        _data[i] = static_cast<T>(0);
}

template <class T>
void FloatVolumeGrid<T>::_freeMemory()
{
    delete [] _data;
}

template <class T>
FloatVolumeGrid<T>::~FloatVolumeGrid()
{
    _freeMemory();
}

template class FloatVolumeGrid<float>;
template class FloatVolumeGrid<double>;
}
