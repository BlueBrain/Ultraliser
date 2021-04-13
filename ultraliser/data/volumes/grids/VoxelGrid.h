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
class VoxelGrid : public VolumeGrid
{
public:

    /**
     * @brief Grid
     * @param dimensions
     */
    VoxelGrid(const Vec3i_64& dimensions,
              const bool& preAllocateMemory = true);

    /**
     * @brief Grid
     * @param width
     * @param height
     * @param depth
     */
    VoxelGrid(const int64_t &width,
              const int64_t &height,
              const int64_t &depth,
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
     * @brief loadByteVolumeData
     * @param prefix
     */
    void loadByteVolumeData(const std::string &prefix) override;

    /**
     * @brief getNumberBytes
     * @return
     */
    uint64_t getNumberBytes() const override;

    /**
     * @brief value
     * @param index
     * @return
     */
    uint8_t getValue(const uint64_t &index) const override;

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

    Voxels getGridData() const
    {
        return _data;
    }

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
    Voxels _data;
};

}

#endif // ULTRALISER_DATA_VOLUME_VOXEL_GRID_H