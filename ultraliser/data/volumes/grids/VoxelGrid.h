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

#ifndef ULTRALISER_DATA_VOLUME_VOXEL_GRID_H
#define ULTRALISER_DATA_VOLUME_VOXEL_GRID_H

#include <data/volumes/voxels/Voxel.h>
#include <data/volumes/grids/VolumeGrid.h>
#include <data/common/BitArray.h>
#include <data/images/Image.h>
#include <math/Math.h>

namespace Ultraliser
{

/**
 * @brief The VoxelGrid class
 */
template <class T>
class VoxelGrid : public VolumeGrid
{
public:

    /**
     * @brief Grid
     * @param dimensions
     */
    VoxelGrid(const Vec3ui_64 &dimensions,
              const bool& preAllocateMemory = true);

    /**
     * @brief Grid
     * @param width
     * @param height
     * @param depth
     */
    VoxelGrid(const uint64_t &width,
              const uint64_t &height,
              const uint64_t &depth,
              const bool& preAllocateMemory = true);

    VoxelGrid(const VoxelGrid* inputGrid);

    ~VoxelGrid();

public:

    /**
     * @brief loadBinaryVolumeData
     * @param prefix
     */
    void loadBinaryVolumeData(const std::string &prefix) override;

    /**
     * @brief loadUnsignedVolumeData
     * @param prefix
     */
    void loadUnsignedVolumeData(const std::string &rawvolumepath) override;

    /**
     * @brief getNumberBytes
     * @return
     */
    uint64_t getNumberBytes() const override;




//    /**
//     * @brief value
//     * @param index
//     * @return
//     */
//    uint8_t getValue(const uint64_t &index) const override;

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
     * @brief writeRAW
     * @param prefix
     */
    void writeRAW(const std::string &prefix) override;

    /**
     * @brief writeBIN
     * @param prefix
     */
    void writeBIN(const std::string &prefix) override;

    /**
     * @brief writeNRRD
     * @param prefix
     */
    void writeNRRD(const std::string &prefix) override;

    void andWithAnotherGrid(VolumeGrid *anotherGrid) override;

    void orWithAnotherGrid(VolumeGrid *anotherGrid) override;

    std::vector< Voxel< T > > getGridData() const
    {
        return _data;
    }






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


private:

    /**
     * @brief _allocateMemory
     */
    void _allocateMemory();

    /**
     * @brief _freeMemory
     */
    void _freeMemory();

    /**
     * @brief _writeHeader
     * @param prefix
     */
    void _writeHeader(const std::string &prefix) override;

private:

    /**
     * @brief _data
     */
    std::vector< Voxel< T > > _data;
};

}

#endif // ULTRALISER_DATA_VOLUME_VOXEL_GRID_H
