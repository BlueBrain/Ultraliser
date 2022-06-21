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
        // Each voxel is stored in a single bit (class 1)
        BIT,

        // Each voxel is stored in a single byte (8 bits) as uint8_t (class 2)
        UI8,

        // Each voxel is stored in a single word (16 bits) as uint16_t (class 2)
        UI16,

        // Each voxel is stored in a doubleword (32 bits) as uint32_t (class 2)
        UI32,

        // Each voxel is stored in a quad-word (64 bits) as uint64_t (class 2)
        UI64,

        // Each voxel is stored as a single-precision floating-point or float value (class 3)
        F32,

        // Each voxel is stored as a double-precision floating-point or double value (class 3)
        F64,

        // Undefined type
        UNDEFINED
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
     * @brief getVolumeTypeFromHdrFile
     * Parses the header file of the volume and returns its type
     * @param filePrefix
     * The prefix of the input file.
     * @return
     * Volume type.
     */
    static TYPE getVolumeTypeFromHdrFile(const std::string& filePrefix);

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
    VolumeGrid(const Vec3ui_64 &dimensions);

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
    VolumeGrid(const uint64_t &width, const uint64_t &height, const uint64_t &depth);

    virtual ~VolumeGrid();

public:




    /**
     * @brief getDimensions
     * Gets the dimensions of the volume.
     * @return
     * The dimensions of the volume in a Vec3 vector format.
     */
    Vec3ui_64 getDimensions() const;

    /**
     * @brief getDimension
     * Gets the dimension of the volume along a given axis.
     * @param i
     * The index of the axis.
     * @return
     * The dimension of the volume along a given axis.
     */
    uint64_t getDimension(const uint64_t &i) const;

    /**
     * @brief mapToIndex
     * Maps a given voxel 3D index into a 1D index.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The one-dimensional index of the data from a three-dimensional voxels.
     */
    uint64_t mapToIndex(const uint64_t &x, const uint64_t &y, const uint64_t &z,
                        bool &outlier) const;

    /**
     * @brief getWidth
     * Gets the width of the volume.
     * @return
     * The width of the volume.
     */
    uint64_t getWidth() const;

    /**
     * @brief getHeight
     * Gets the height of the volume.
     * @return
     * The height of the volume.
     */
    uint64_t getHeight() const;

    /**
     * @brief getDepth
     * Gets the depth of the volume.
     * @return
     * The depth of the volume.
     */
    uint64_t getDepth() const;

    /**
     * @brief getNumberVoxels
     * Gets the number of voxels in the volume.
     * @return
     * The number of voxels in the volume.
     */
    uint64_t getNumberVoxels() const;

    /**
     * @brief getNumberBytes
     * Gets the number of bytes used to store the volume data.
     * @return
     * The number of bytes used to store the volume data.
     */
    virtual uint64_t getNumberBytes() const = 0;

    /**
     * @brief loadBinaryVolumeData
     * @param prefix
     */
    virtual void loadBinaryVolumeData(const std::string &prefix) = 0;

    /**
     * @brief loadUnsignedVolumeData
     * @param prefix
     */
    virtual void loadUnsignedVolumeData(const std::string &rawvolumepath) = 0;


    virtual void readUVOLBData(const std::string &prefix) = 0;

    virtual void readUVOLData(const std::string &prefix) = 0;







    /**
     * @brief clear
     * Clears the voxel data of the volume.
     */
    virtual void clear() = 0;

    /**
     * @brief fillVoxel
     * Fills a given voxel in the volume specified by its three-dimensional indices.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
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
    // virtual uint8_t getValue(const uint64_t &index) const = 0;

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
     * @brief getValueUI8
     * Returns the value of a voxel specified by a given index as an 8-bit integer.
     * If the volume has 16-, 32-, 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as an 8-bit integer.
     */
    virtual uint8_t getValueUI8(const uint64_t &index) const = 0;

    /**
     * @brief getValueUI16
     * Returns the value of a voxel specified by a given index as a 16-bit integer.
     * If the volume has 32-, 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 16-bit integer.
     */
    virtual uint16_t getValueUI16(const uint64_t &index) const = 0;

    /**
     * @brief getValueUI32
     * Returns the value of a voxel specified by a given index as a 32-bit integer.
     * If the volume has 64-bit volume the return value is zero.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 32-bit integer.
     */
    virtual uint32_t getValueUI32(const uint64_t &index) const = 0;

    /**
     * @brief getValueUI64
     * Returns the value of a voxel specified by a given index as a 64-bit integer.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a 64-bit integer.
     */
    virtual uint64_t getValueUI64(const uint64_t &index) const = 0;

    /**
     * @brief getValueF32
     * Returns the value of a voxel specified by a given index as a single-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a single-precision float.
     */
    virtual float getValueF32(const uint64_t &index) const = 0;

    /**
     * @brief getValueF64
     * Returns the value of a voxel specified by a given index as a double-precision float.
     * @param index
     * The one-dimensional index of the voxel.
     * @return
     * The value of the voxel as a double-precision float.
     */
    virtual double getValueF64(const uint64_t &index) const = 0;

    /**
     * @brief getValueUI8
     * Returns the value of a voxel specified by a given index as an 8-bit integer.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as an 8-bit integer.
     */
    uint8_t getValueUI8(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI16
     * Returns the value of a voxel specified by a given index as a 16-bit integer.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as a 16-bit integer.
     */
    uint16_t getValueUI16(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI32
     * Returns the value of a voxel specified by a given index as a 32-bit integer.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as a 32-bit integer.
     */
    uint32_t getValueUI32(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI64
     * Returns the value of a voxel specified by a given index as a 64-bit integer.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as a 64-bit integer.
     */
    uint64_t getValueUI64(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI16
     * Returns the value of a voxel specified by a given index as a single-precision float.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as a single-precision float.
     */
    float getValueF32(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueF64
     * Returns the value of a voxel specified by a given index as a double-precision float.
     * @param x
     * X-axis index.
     * @param y
     * Y-axis index.
     * @param z
     * Z-axis index.
     * @return
     * The value of the voxel as a double-precision float.
     */
    double getValueF64(const int64_t &x, const int64_t &y, const int64_t &z) const;

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

    /**
     * @brief writeProjection
     * @param prefix
     * @param projection
     * @param projectColorCoded
     */
    void writeProjection(const std::string &prefix,
                         const PROJECTION &projection,
                         const bool &projectColorCoded);

    /**
     * @brief writeBIN
     * Writes a binary volume (1 bit per voxel) of the grid.
     * The dimensions of the volume grid will be written to a .HDR file, while the data will be
     * written to a separate .BIN file.
     * @param prefix
     * File prefix.
     */
    virtual void writeBIN(const std::string &prefix) = 0;

    /**
     * @brief writeRAW
     * Writes a raw volume (1 byte per voxel) of the grid.
     * The dimensions of the volume grid will be written to a .HDR file, while the data will be
     * written to a separate .RAW file.
     * @param prefix
     * File prefix.
     */
    virtual void writeRAW(const std::string &prefix) = 0;

    /**
     * @brief writeNRRD
     * Writes an NRRD file (8-bit raw) of the grid.
     * The NRRD file can be read with Paraview, where the dimensions and data are integrated
     * in the same file.
     * @param prefix
     * File prefix.
     */
    virtual void writeNRRD(const std::string &prefix) = 0;

    /**
     * @brief writeUltraliserBinaryVolume
     * Writes an Ultraliser-specific binary volume file (1 bit per voxel).
     * The created file contains the type of the file ('1bit'), the dimensions of the file
     * in ('x y z') format and the data of the volume grid in a binary format (1 bit per voxel).
     * @param prefix
     * File prefix.
     */
    virtual void writeUltraliserBinaryVolume(const std::string &path) = 0;

    /**
     * @brief writeUltraliserRawVolume
     * Writes an Ultraliser-specific unsigned volume file (8-, 16-, 32-, or 64-bit file depending
     * on the type of the volume grid itself).
     * The created file contains the type of the file ('8ui, 16ui, 32ui or 64ui'), the dimensions
     * of the file in ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    virtual void writeUltraliserRawVolume(const std::string &path) = 0;

    /**
     * @brief writeUltraliserFloatVolume
     * Writes an Ultraliser-specific unsigned volume file (32-, or 64-bit precision volume files
     * depending on the type of the volume grid itself).
     * The created file contains the type of the file ('32f or 64f'), the dimensions of the file in
     * ('x y z') format and the data of the volume grid.
     * @param prefix
     * File prefix.
     */
    virtual void writeUltraliserFloatVolume(const std::string &path) = 0;






    /**
     * @brief floodFillSliceAlongAxis
     * @param x
     * @param axis
     * @param padding
     */
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
     * Vector containing the dimensions of the grid along each axis.
     */
    Vec3ui_64 _dimensions;

    /**
     * @brief _numberVoxels
     * Total number of voxels in the grid.
     */
    uint64_t _numberVoxels;

    /**
     * @brief _numberBytes
     * Total number of bytes in the grid.
     */
    uint64_t _numberBytes;

    /**
     * @brief _projectionTime
     * The time it takes to project the volume grid into a 2D image.
     */
    double _projectionTime;
};

}

#endif // ULTRALISER_DATA_VOLUME_VOLUME_GRID_H
