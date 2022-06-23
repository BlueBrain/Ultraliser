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

#pragma once

#include <math/Math.h>
#include <data/common/BitArray.h>
#include <data/images/Image.h>
#include <data/volumes/grids/VolumeGrid.h>

namespace Ultraliser
{

/**
 * @brief The FloatVolumeGrid class
 */
template <class T>
class FloatVolumeGrid : public VolumeGrid
{
public:

    /**
     * @brief Grid
     * @param width
     * @param height
     * @param depth
     */
    FloatVolumeGrid(const uint64_t &width,
                   const uint64_t &height,
                   const uint64_t &depth,
                   const bool& preAllocateMemory = true);

    /**
     * @brief FloatVolumeGrid
     * @param inputGrid
     */
    FloatVolumeGrid(const FloatVolumeGrid* inputGrid);

    ~FloatVolumeGrid();

public:
    /**
     * @brief getNumberBytes
     * @return
     */
    uint64_t getNumberBytes() const override;

    /**
     * @brief getByte
     * @param index
     * @return
     */
    uint8_t getByte(uint64_t index) const override;

    /**
     * @brief addByte
     * @param index
     * @param byte
     */
    void addByte(const uint64_t &index, const uint8_t &byte) override;

    /**
     * @brief clear
     */
    void clear() override;

    /**
     * @brief fillVoxel
     * @param index
     */
    void fillVoxel(const uint64_t &index) override;

    /**
     * @brief clearVoxel
     * @param index
     */
    void clearVoxel(const uint64_t &index) override;

    /**
     * @brief isFilled
     * @param index
     * @return
     */
    bool isFilled(const uint64_t &index) const override;

    /**
     * @brief isEmpty
     * @param index
     * @return
     */
    bool isEmpty(const uint64_t &index) const override;


    /**
     * @brief andWithAnotherGrid
     * @param anotherGrid
     */
    void andWithAnotherGrid(VolumeGrid *anotherGrid) override;

    /**
     * @brief orWithAnotherGrid
     * @param anotherGrid
     */
    void orWithAnotherGrid(VolumeGrid *anotherGrid) override;







    /**
     * @brief getGridData
     * @return
     */
    T* getGridData() const { return _data; }

    /**
     * @brief getValueUI8
     * Returns the value of a voxel specified by a given index as an 8-bit integer.
     * If the volume has 16-, 32-, 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as an 8-bit integer.
     */
    uint8_t getValueUI8(const uint64_t &index) const;

    /**
     * @brief getValueUI16
     * Returns the value of a voxel specified by a given index as a 16-bit integer.
     * If the volume has 32-, 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 16-bit integer.
     */
    uint16_t getValueUI16(const uint64_t &index) const;

    /**
     * @brief getValueUI32
     * Returns the value of a voxel specified by a given index as a 32-bit integer.
     * If the volume has 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 32-bit integer.
     */
    uint32_t getValueUI32(const uint64_t &index) const;

    /**
     * @brief getValueUI64
     * Returns the value of a voxel specified by a given index as a 64-bit integer.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 64-bit integer.
     */
    uint64_t getValueUI64(const uint64_t &index) const;

    /**
     * @brief getValueF32
     * Returns the value of a voxel specified by a given index as a single-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a single-precision float.
     */
    float getValueF32(const uint64_t &index) const;

    /**
     * @brief getValueF64
     * Returns the value of a voxel specified by a given index as a double-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a double-precision float.
     */
    double getValueF64(const uint64_t &index) const;


    /**
     * @brief writeBitVolume
     * Writes an Ultraliser-specific binary volume file (1 bit per voxel).
     * The created file contains the type of the file ('1bit'), the dimensions of the file
     * in ('x y z') format and the data of the volume grid in a binary format (1 bit per voxel).
     * @param path
     * Absolute file path.
     */
    void writeBitVolume(const std::string &prefix) const override;

    /**
     * @brief writeUnsignedVolume
     * Writes an Ultraliser-specific unsigned volume file (8-, 16-, 32-, or 64-bit file depending
     * on the type of the volume grid itself).
     * The created file contains the type of the file ('8ui, 16ui, 32ui or 64ui'), the dimensions
     * of the file in ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    void writeUnsignedVolume(const std::string &prefix) const override;

    /**
     * @brief writeFloatVolume
     * Writes an Ultraliser-specific unsigned volume file (32-, or 64-bit precision volume files
     * depending on the type of the volume grid itself).
     * The created file contains the type of the file ('32f or 64f'), the dimensions of the file in
     * ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    void writeFloatVolume(const std::string &prefix) const override;

    /**
     * @brief writeNRRDVolume
     * @param prefix
     */
    void writeNRRDVolume(const std::string &prefix) const override;

    /**
     * @brief writeRAWVolume
     * @param prefix
     */
    void writeRAWVolume(const std::string &prefix) const override;

private:

    /**
     * @brief _allocateMemory
     */
    void _allocateMemory();

    /**
     * @brief _freeMemory
     */
    void _freeMemory();

private:

    /**
     * @brief _data
     */
    T* _data;
};

/**
 * @brief VolumeGridF32
 */
typedef FloatVolumeGrid< float > VolumeGridF32;

/**
 * @brief VolumeGridF64
 */
typedef FloatVolumeGrid< double > VolumeGridF64;
}

