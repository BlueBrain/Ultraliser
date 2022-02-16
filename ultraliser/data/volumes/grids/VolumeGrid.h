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

#ifndef ULTRALISER_DATA_VOLUME_VOLUME_GRID_H
#define ULTRALISER_DATA_VOLUME_VOLUME_GRID_H

#include <math/Math.h>
#include <data/common/BitArray.h>
#include <data/images/Image.h>

namespace Ultraliser
{

/**
 * @brief The Grid class
 */
class VolumeGrid
{
public:

    /**
     * @brief The TYPE enum
     */
    enum TYPE
    {
        // Each voxel is stored in a single bit
        BIT,

        // Each voxel is stored in a single byte
        BYTE,

        // Each voxel is stored in a Voxel structure
        VOXEL
    };

    enum PROJECTION
    {
        // Project the volume along the Z-axis
        XY_PROJECTION,

        // Project the volume along the Y-axis
        XZ_PROJECTION,

        // Project the volume along the X-axis
        ZY_PROJECTION
    };

    /**
     * @brief getType
     * @return
     */
    static TYPE getType(const std::string &typeString);

    /**
     * @brief getTypeString
     * @param type
     * @return
     */
    static std::string getTypeString(const VolumeGrid::TYPE &type);

public:

    /**
     * @brief VolumeGrid
     * @param dimensions
     */
    VolumeGrid(const Vec3i_64 &dimensions);

    /**
     * @brief VolumeGrid
     * Copy constructor.
     *
     * @param grid
     * Input grid.
     */
    VolumeGrid(VolumeGrid *grid);

    /**
     * @brief VolumeGrid
     * @param width
     * @param height
     * @param depth
     */
    VolumeGrid(const int64_t &width, const int64_t &height, const int64_t &depth);

    virtual ~VolumeGrid();

public:

    /**
     * @brief getDimensions
     * @return
     */
    Vec3i_64 getDimensions() const;

    /**
     * @brief getWidth
     * @return
     */
    int64_t getWidth() const;

    /**
     * @brief getHeight
     * @return
     */
    int64_t getHeight() const;

    /**
     * @brief getDepth
     * @return
     */
    int64_t getDepth() const;

    /**
     * @brief getNumberVoxels
     * @return
     */
    int64_t getNumberVoxels() const;

    /**
     * @brief loadBinaryVolumeData
     * @param prefix
     */
    virtual void loadBinaryVolumeData(const std::string &prefix) = 0;

    /**
     * @brief loadByteVolumeData
     * @param prefix
     */
    virtual void loadByteVolumeData(const std::string &prefix) = 0;

    /**
     * @brief getNumberBytes
     * @return
     */
    virtual uint64_t getNumberBytes() const = 0;

    /**
     * @brief size
     * @param dimension
     * @return
     */
    int64_t getDimension(const int32_t &i) const;

    /**
     * @brief mapToIndex
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint64_t mapToIndex(const int64_t &x, const int64_t &y, const int64_t &z, bool &outlier) const;

    /**
     * @brief clear
     */
    virtual void clear() = 0;

    /**
     * @brief fillVoxel
     * @param x
     * @param y
     * @param z
     */
    void fillVoxel(const int64_t &x, const int64_t &y, const int64_t &z);

    /**
     * @brief fillVoxel
     * @param index
     */
    virtual void fillVoxel(const uint64_t &index)= 0;

    /**
     * @brief clearVoxel
     * @param x
     * @param y
     * @param z
     */
    void clearVoxel(const int64_t &x, const int64_t &y, const int64_t &z);

    /**
     * @brief clearVoxel
     * @param index
     */
    virtual void clearVoxel(const uint64_t &index) = 0;

    /**
     * @brief isFilled
     * @param index
     * @return
     */
    virtual bool isFilled(const uint64_t &index) const = 0;

    /**
     * @brief isFilled
     * @param x
     * @param y
     * @param z
     * @return
     */
    bool isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief isEmpty
     * @param index
     * @return
     */
    virtual bool isEmpty(const uint64_t &index) const = 0;

    /**
     * @brief isEmpty
     * @param x
     * @param y
     * @param z
     * @return
     */
    bool isEmpty(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief value
     * @param index
     * @return
     */
    virtual uint8_t getValue(const uint64_t &index) const = 0;

    /**
     * @brief andWithAnotherGrid
     * @param anotherGrid
     */
    virtual void andWithAnotherGrid(VolumeGrid *anotherGrid) = 0;

    /**
     * @brief orWithAnotherGrid
     * @param anotherGrid
     */
    virtual void orWithAnotherGrid(VolumeGrid *anotherGrid) = 0;

    /**
     * @brief value
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getValue(const int64_t &x, const int64_t &y, const int64_t &z) const
    {
        bool outlier;
        uint64_t index = mapToIndex(x, y, z, outlier);
        if (outlier)
            return 0;
        else
            return getValue(index);
    }

    /**
     * @brief projectVolume
     * @param prefix
     * @param projectXYView
     * @param projectZYView
     */
    void projectVolume(const std::string &prefix,
                               const bool &projectXY = false,
                               const bool &projectXZ = false,
                               const bool &projectZY = false,
                               const bool &projectColorCoded = false);

    void writeProjection(const std::string &prefix,
                                 const PROJECTION &projection,
                                 const bool &projectColorCoded);

    /**
     * @brief writeRAW
     * @param prefix
     */
    virtual void writeRAW(const std::string &prefix) = 0;

    /**
     * @brief writeBIN
     * @param prefix
     */
    virtual void writeBIN(const std::string &prefix) = 0;

    /**
     * @brief writeNRRD
     * @param prefix
     */
    virtual void writeNRRD(const std::string &prefix) = 0;

    void floodFillSliceAlongAxis(const int64_t &x, const AXIS &axis, const uint64_t &padding = 0);

    /**
     * @brief getByte
     * @param index
     * @return
     */
    virtual uint8_t getByte(uint64_t index) const = 0;

    /**
     * @brief addByte
     * @param index
     * @param byte
     */
    virtual void addByte(const uint64_t &index, const uint8_t &byte) = 0;

    /**
     * @brief computeNumberNonZeroVoxelsPerSlice
     * @param z
     * @return
     */
    virtual uint64_t computeNumberNonZeroVoxelsPerSlice(int64_t z) const;

    /**
     * @brief computeNumberNonZeroVoxels
     * @return
     */
    virtual uint64_t computeNumberNonZeroVoxels() const;

protected:

    /**
     * @brief _writeHeader
     * Writes the header file of the data.
     * @param prefix
     */
    virtual void _writeHeader(const std::string &prefix) = 0;

protected:

    /**
     * @brief _dimensions
     * The dimensions of the grid.
     */
    Vec3i_64 _dimensions;

    /**
     * @brief _numberVoxels
     * Number of voxels in the grid.
     */
    int64_t _numberVoxels;

    /**
     * @brief _numberBytes
     * Number of bytes in the grid.
     */
    int64_t _numberBytes;

    /**
     * @brief _projectionTime
     */
    double _projectionTime;
};

}

#endif // ULTRALISER_DATA_VOLUME_VOLUME_GRID_H
