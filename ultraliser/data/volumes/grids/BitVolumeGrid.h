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

#pragma once

#include <math/Math.h>
#include <data/images/Image.h>
#include <data/common/BitArray.h>
#include <data/volumes/grids/VolumeGrid.h>

namespace Ultraliser
{

/**
 * @brief The BitVolumeGrid class
 */
class BitVolumeGrid : public VolumeGrid
{
public:

    /**
     * @brief BitVolumeGrid
     * @param width
     * @param height
     * @param depth
     * @param data
     */
    BitVolumeGrid(const size_t &width,
                  const size_t &height,
                  const size_t &depth,
                  BitArray *data);

    /**
     * @brief Grid
     * Constructor
     * @param width
     * The width of the volume.
     * @param height
     * The height of the volume.
     * @param depth
     * The depth of the volume.
     * @param preAllocateMemory
     * If this flag is set, the memory of the grid will be pre-allocated during the
     * construction of the object.
     */
    BitVolumeGrid(const int64_t &width,
                  const int64_t &height,
                  const int64_t &depth,
                  const bool& preAllocateMemory = true);

    /**
     * @brief BitVolumeGrid
     * Copy constructor.
     * @param inputGrid
     * An input grid to be copied.
     */
    BitVolumeGrid(const BitVolumeGrid* inputGrid);

    /// Destructor
    ~BitVolumeGrid();

public:

    /**
     * @brief getNumberBytes
     * Gets the nummber of memory bytes of the volume grid.
     * @return
     * Returns the nummber of memory bytes of the volume grid.
     */
    uint64_t getNumberBytes() const override;

    /**
     * @brief getValueUI8
     * Returns the value of a voxel specified by a given index as an 8-bit integer.
     * If the volume is a 16-, 32-, 64-bit volume, the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as an 8-bit integer.
     */
    uint8_t getValueUI8(const uint64_t &index) const override;

    /**
     * @brief getValueUI16
     * Returns the value of a voxel specified by a given index as a 16-bit integer.
     * If the volume is a 32-, 64-bit volume, the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 16-bit integer.
     */
    uint16_t getValueUI16(const uint64_t &index) const override;

    /**
     * @brief getValueUI32
     * Returns the value of a voxel specified by a given index as a 32-bit integer.
     * If the volume is a 64-bit volume, the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 32-bit integer.
     */
    uint32_t getValueUI32(const uint64_t &index) const override;

    /**
     * @brief getValueUI64
     * Returns the value of a voxel specified by a given index as a 64-bit integer.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 64-bit integer.
     */
    uint64_t getValueUI64(const uint64_t &index) const override;

    /**
     * @brief getValueF32
     * Returns the value of a voxel specified by a given index as a single-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a single-precision float.
     */
    float getValueF32(const uint64_t &index) const override;

    /**
     * @brief getValueF64
     * Returns the value of a voxel specified by a given index as a double-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a double-precision float.
     */
    double getValueF64(const uint64_t &index) const override;

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
     * @brief writeBIN
     * Writes a binary volume (1 bit per voxel) of the grid.
     * The dimensions of the volume grid will be written to a .HDR file, while the data will be
     * written to a separate .BIN file.
     * @param prefix
     * File prefix.
     */
    void writeBIN(const std::string &prefix) override;

    /**
     * @brief writeRAW
     * Writes a raw volume (1 byte per voxel) of the grid.
     * The dimensions of the volume grid will be written to a .HDR file, while the data will be
     * written to a separate .RAW file.
     * @param prefix
     * File prefix.
     */
    void writeRAW(const std::string &prefix) override;

    /**
     * @brief writeNRRD
     * Writes an NRRD file (8-bit raw) of the grid.
     * The NRRD file can be read with Paraview, where the dimensions and data are integrated
     * in the same file.
     * @param prefix
     * File prefix.
     */
    void writeNRRD(const std::string &prefix) override;

    /**
     * @brief writeUltraliserBinaryVolume
     * Writes an Ultraliser-specific binary volume file (1 bit per voxel).
     * The created file contains the type of the file ('1bit'), the dimensions of the file
     * in ('x y z') format and the data of the volume grid in a binary format (1 bit per voxel).
     * @param path
     * Absolute file path.
     */
    void writeUltraliserBinaryVolume(const std::string &prefix) override;

    /**
     * @brief writeUltraliserRawVolume
     * Writes an Ultraliser-specific unsigned volume file (8-, 16-, 32-, or 64-bit file depending
     * on the type of the volume grid itself).
     * The created file contains the type of the file ('8ui, 16ui, 32ui or 64ui'), the dimensions
     * of the file in ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    void writeUltraliserRawVolume(const std::string &prefix) override;

    /**
     * @brief writeUltraliserFloatVolume
     * Writes an Ultraliser-specific unsigned volume file (32-, or 64-bit precision volume files
     * depending on the type of the volume grid itself).
     * The created file contains the type of the file ('32f or 64f'), the dimensions of the file in
     * ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    void writeUltraliserFloatVolume(const std::string &prefix) override;

    /**
     * @brief getGridData
     * @return
     */
    BitArray* getGridData() const { return _data; }

private:

    /**
     * @brief _writeHeader
     * @param prefix
     */
    void _writeHeader(const std::string &prefix) override;

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
    BitArray* _data;
};

}
