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

#ifndef ULTRALISER_DATA_VOLUME_H
#define ULTRALISER_DATA_VOLUME_H

#include <common/Common.h>
#include <data/common/GridIndex.h>
#include <data/volumes/grids/BitVolumeGrid.h>
#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/morphologies/Morphologies.h>
#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{

/**
 * @brief Te Volume class
 * Builds a voxel grid for a mesh where 0 stands for empty cell and 1 stands
 * for a filledd cell.
 */
class Volume
{

public:
    enum SOLID_VOXELIZATION_AXIS
    {
        X = 0,
        Y = 1,
        Z = 2,
        XYZ = 3
    };

public:

    /**
     * @brief Volume
     * Constructor
     *
     * @param pMin
     *
     * @param pMax
     *
     * @param baseResolution
     * The base resolution of the volume.
     * By default, this resolution is set to 512.
     *
     * @param voxelPadding
     * Number of voxels used to pad the volume.
     * By default, the volume has zero-padding.
     */
    Volume(const Vector3f& pMin,
           const Vector3f& pMax,
           const uint64_t &baseResolution = 512,
           const float &voxelPadding = 0.0,
           const VolumeGrid::TYPE& gridType = VolumeGrid::TYPE::BIT);

    /**
     * @brief Volume
     * Constructor
     * @param width
     * @param height
     * @param depth
     */
    Volume(const int64_t width,
           const int64_t height,
           const int64_t depth,
           const Vector3f pMin,
           const Vector3f pMax,
           const VolumeGrid::TYPE& gridType = VolumeGrid::TYPE::BIT);

    /**
     * @brief Volume
     * Constructor
     * @param prefix
     * @param type
     */
    Volume(const std::string &prefix,
           const VolumeGrid::TYPE& gridType = VolumeGrid::TYPE::BIT);
    ~Volume();

    /**
     * @brief getByte
     * @param index
     * @return
     */
    uint8_t getByte(const uint64_t index) const;

    /**
     * @brief getValue
     * @param index
     * @return
     */
    uint8_t getValue(const uint64_t index) const;

    /**
     * @brief getValue
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getValue(const int64_t &x,
                     const int64_t &y,
                     const int64_t &z) const;

    /**
     * @brief fillVoxel
     * @param x
     * @param y
     * @param z
     */
    void fillVoxel(const int64_t &x,
                   const int64_t &y,
                   const int64_t &z);

    /**
     * @brief getByte
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getByte(const int64_t &x,
                    const int64_t &y,
                    const int64_t &z) const;

    /**
     * @brief addByte
     * @param index
     * @param byte
     */
    void addByte(const uint64_t &index, const uint8_t byte);

    /**
     * @brief clear
     */
    void clear(void);

    /**
     * @brief surfaceVoxelization
     * Performs a surface voxelization on the input mesh.
     * @param mesh
     */
    void surfaceVoxelization(Mesh* mesh,
                             const bool& verbose = false,
                             const bool parallel = false);

    /**
     * @brief surfaceVoxelization
     * Surface voxelization of advanced mesh.
     * @param mesh
     */
    void surfaceVoxelization(AdvancedMesh* mesh);

    /**
     * @brief surfaceVoxelization
     * Performs a surface voxelization on a list of input meshes in parallel.
     * @param meshFiles
     */
    void surfaceVoxelization(const std::string &inputDirectory,
                             const std::vector< std::string>& meshFiles);

    /**
     * @brief surfaceVoxelizeVasculatureMorphology
     * @param morphology
     */
    void surfaceVoxelizeVasculatureMorphologyParallel(VasculatureMorphology* vasculatureMorphology);

    /**
     * @brief solidVoxelization
     * lood fill the exterior to figure out the interior.
     * @param algorithm
     */
    void solidVoxelization(const SOLID_VOXELIZATION_AXIS& axis=X);

    /**
     * @brief exportToMesh
     * @param prefix
     */
    void exportToMesh(const std::string &prefix,
                      const bool &formatOBJ = false,
                      const bool &formatPLY = false,
                      const bool &formatOFF = false,
                      const bool &formatSTL = false) const;

    /**
     * @brief writeStackXY
     * @param prefix
     */
    void writeStackXY(const std::string &outputDirectory,
                        const std::string &prefix) const;

    /**
     * @brief writeStackXZ
     * @param outputDirectory
     * @param prefix
     */
    void writeStackXZ(const std::string &outputDirectory,
                        const std::string &prefix) const;

    /**
     * @brief writeStackZY
     * @param outputDirectory
     * @param prefix
     */
    void writeStackZY(const std::string &outputDirectory,
                        const std::string &prefix) const;

    /**
     * @brief writeVolumes
     * @param prefix
     * @param binaryFormat
     * @param rawFormat
     */
    void writeVolumes(const std::string &prefix,
                      const bool& binaryFormat = false,
                      const bool& rawFormat = false,
                      const bool& nrrdFormat = false) const;

    /**
     * @brief writeStacks
     * @param outputDirectory
     * @param prefix
     * @param xy
     * @param zy
     */
    void writeStacks(const std::string &outputDirectory,
                     const std::string &prefix,
                     const bool& xy = false,
                     const bool& xz = false,
                     const bool& zy = false) const;

    /**
     * @brief project
     * @param prefix
     * @param xy
     * @param zy
     */
    void project(const std::string prefix,
                 const bool xy = false,
                 const bool xz = false,
                 const bool zy = false,
                 const bool &projectColorCoded = false) const;

    /**
     * @brief getWidth
     * @return
     */
    int64_t getWidth(void) const;

    /**
     * @brief getHeight
     * @return
     */
    int64_t getHeight(void) const;

    /**
     * @brief getDepth
     * @return
     */
    int64_t getDepth(void) const;

    /**
     * @brief getNumberVoxels
     * @return
     */
    uint64_t getNumberVoxels(void) const;

    /**
     * @brief computeNumberNonZeroVoxels
     * @return
     */
    uint64_t computeNumberNonZeroVoxels(void) const;

    /**
     * @brief computeVolume
     * @return
     */
    float computeVolume() const;

    /**
     * @brief getNumberBytes
     * @return
     */
    uint64_t getNumberBytes(void) const;

    /**
     * @brief getSurfaceVoxelizationTime
     * @return
     */
    double getSurfaceVoxelizationTime(void) const;

    /**
     * @brief getSolidVoxelizationTime
     * @return
     */
    double getSolidVoxelizationTime(void) const;

    /**
     * @brief getAppendingVolumePassTime
     * @return
     */
    double getAppendingVolumePassTime(void) const;

    /**
     * @brief printBenchmarks
     * @param mesh
     */
    //void printBenchmarks(OriginalMesh* mesh) const;

    /**
     * @brief addVolume
     * Adds volume with the same dimensions to the current one from a file.
     * @param volumePrefix
     */
    void addVolume(const std::string &volumePrefix);

    /**
     * @brief addVolumePass
     * Adds volume with the same dimensions to the current one.
     * @param volume
     */
    void addVolumePass(const Volume* volume);

    /**
     * @brief mapToIndex
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint64_t mapToIndex(const int64_t &x,
                        const int64_t &y,
                        const int64_t &z) const;

    /**
     * @brief fill
     * @param index
     */
    void fill(const u_int64_t& index);

    /**
     * @brief clear
     * @param index
     */
    void clear(const u_int64_t& index);

    /**
     * @brief fill
     * @param x
     * @param y
     * @param z
     */
    void fill(const int64_t &x, const int64_t &y, const int64_t &z);

    /**
     * @brief clear
     * @param x
     * @param y
     * @param z
     */
    void clear(const int64_t &x, const int64_t &y, const int64_t &z);

    /**
     * @brief isFilled
     * @param index
     * @return
     */
    bool isFilled(const u_int64_t& index) const;

    /**
     * @brief isFilled
     * @param x
     * @param y
     * @param z
     * @return
     */
    bool isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getFormatString
     * @return
     */
    std::string getFormatString() const;

    /**
     * @brief printVolumeStats
     * @param reference
     * @param prefix
     */
    void printStats(const std::string &reference,
                          const std::string* prefix = nullptr) const;

public:

    /**
     * @brief getLargestDimension
     * @param dimensions
     * @return
     */
    static int32_t getLargestDimension(const Vector3f& dimensions);

    /**
     * @brief constructFromTiffMask
     * Construct a volume from a tiff mask.
     * @return
     * A pointer to the reconstructed volume from the mask.
     */
    static Volume* constructFromTiffMask(const std::string &maskDirectory,
                                         const int64_t &maskWidth,
                                         const int64_t &maskHeight,
                                         const VolumeGrid::TYPE &gridType);

    /**
     * @brief constructIsoValueVolume
     * Constructs a binary volume (1 bit per voxel) from a byte volume (1 byte
     * per voxel) for a given iso value that is defined by the user.
     * @param volume
     * A given byte volume.
     * @param isoValue
     * An iso-value between 0 and 255.
     * @param padding
     * Number of voxels that are used to zero-pad the volume, by default
     * 32 in each dimension.
     * @return
     * A binary volume (1 bit per voxel) corresponding to the given iso value.
     */
    static Volume* constructIsoValueVolume(const Volume* volume,
                                           const uint8_t &isoValue,
                                           const int64_t &padding = 32);

    /**
     * @brief constructFullRangeVolume
     * Constructs a binary volume (1 bit per voxel) from a byte volume (1 byte
     * per voxel). If the voxel is filled, then it sets the value of the bit
     * to 1, otherwise it is zero.
     * @param volume
     * A given byte volume.
     * @param padding
     * Number of voxels that are used to zero-pad the volume, by default
     * 32 in each dimension.
     * @return
     * A binary volume (1 bit per voxel) corresponding to the given iso value.
     */
    static Volume* constructFullRangeVolume(const Volume* volume,
                                            const int64_t &padding = 32);

    /**
     * @brief createHistogram
     * Creates the histogram of the given volume.
     * @param volume
     * A given byte volume.
     * @return
     * Histogram array.
     */
    static std::vector<uint64_t> createHistogram(const Volume* volume);

private:

    /**
     * @brief _loadHeaderData
     * @param prefix
     */
    void _loadHeaderData(const std::string &prefix);

    /**
     * @brief _loadByteVolumeData
     * @param prefix
     */
    void _loadByteVolumeData(const std::string &prefix);

    /**
     * @brief _loadBinaryVolumeData
     * @param prefix
     */
    void _loadBinaryVolumeData(const std::string &prefix);

    /**
     * @brief _allocateGrid
     */
    void _allocateGrid();

    /**
     * @brief _createGrid
     * Creates the volume grid and initializes it to zeros.
     */
    void _createGrid(void);

    /**
     * @brief _clampIndex
     * @param idx
     * @param dimension
     * @return
     */
    uint64_t _clampIndex(uint64_t idx, uint64_t dimension);

    /**
     * @brief _getBoundingBox
     * @param mesh
     * @param i
     * @param tMin
     * @param tMax
     */
    void _getBoundingBox(Mesh* mesh ,
                         uint64_t i, int64_t *tMin, int64_t *tMax);

    /**
     * @brief _getTriangleBoundingBox
     * @param triangle
     * @param tMin
     * @param tMax
     */
    void _getTriangleBoundingBox(AdvancedTriangle triangle,
                                 int64_t *tMin, int64_t *tMax);

    /**
     * @brief _rasterize
     * @param mesh
     * @param grid
     */
    void _rasterize(Mesh* mesh, VolumeGrid* grid, const bool& verbose = false);

    /**
     * @brief _rasterizeParallel
     * @param mesh
     * @param grid
     * @param verbose
     */
    void _rasterizeParallel(Mesh* mesh, VolumeGrid* grid);

    /**
     * @brief floodFill2D
     * Performs a 2D flood-filling operation.
     */
    void _floodFill2D(const SOLID_VOXELIZATION_AXIS& axis);

    void _floodFillAlongAxis(VolumeGrid* grid, const SOLID_VOXELIZATION_AXIS &axis);

    /**
     * @brief _floodFillAlongXYZ
     * @param grid
     */
    void _floodFillAlongXYZ(VolumeGrid *grid);

    /**
     * @brief floodFill3D
     * Performs a 3D flood-filling operation.
     */
    void _floodFill3D(void);

    /**
     * @brief _rasterize
     * @param mesh
     * @param grid
     */
    void _rasterize(AdvancedMesh *mesh, VolumeGrid* grid);

    /**
     * @brief vec2grid
     * Computes grid index given coordinate
     * @param v
     * @param grid
     */
    void _vec2grid(const Vector3f &v, GridIndex& grid);

    /**
     * @brief _testTriangleCubeIntersection
     * @param mesh
     * @param tIdx
     * @param voxel
     * @return
     */
    bool _testTriangleCubeIntersection(Mesh* mesh,
                                       uint64_t tIdx,
                                       const GridIndex& voxel);

    /**
     * @brief _testTriangleGridIntersection
     * @param triangle
     * @param voxel
     * @return
     */
    bool _testTriangleGridIntersection(AdvancedTriangle triangle,
                                       const GridIndex& voxel);

    /**
     * @brief _triangleCubeSign
     * @param mesh
     * @param tIdx
     * @param _boundingBox
     * @return 2 if cube center is on the positive half space (outside)
     * of triangle.
     * @return 1 of tje cube center is on the negative side, but the
     * projection is outside of the triangle.
     */
    int _triangleCubeSign(Mesh* mesh, int tIdx,
                          const GridIndex& boundingBox);

private:

    /**
     * @brief grid
     */
    VolumeGrid* _grid;

    /**
     * @brief _gridType
     */
    const VolumeGrid::TYPE _gridType;

    /**
     * @brief _pMin
     */
    Vector3f _pMin;

    /**
     * @brief _pMax
     */
    Vector3f _pMax;

    /**
     * @brief voxelPadding
     * Additional layer of voxels around the bounding box of the object.
     */
    float _voxelPadding;

    /**
     * @brief _gridDimensions
     * The dimensions of the volume grid.
     */
    Vec3i_64 _gridDimensions;

    /**
     * @brief _baseResolution
     * Base resolution of the volume.
     */
    uint64_t _baseResolution;

    /**
     * @brief _largestDimensionIdx
     * The index of the largest dimension
     */
    int _largestDimensionIdx;

    /**
     * @brief _voxelSize
     * The size of the voxel.
     */
    float _voxelSize;

    /**
     * @brief _volumeOrigin
     * The origin of the mesh model.
     */
    Vector3f _meshOrigin;

    /**
     * @brief _surfaceVoxelizationTime
     */
    double _surfaceVoxelizationTime;

    /**
     * @brief _solidVoxelizationTime
     */
    double _solidVoxelizationTime;

    /**
     * @brief _solidVoxelizationTime
     */
    double _addingVolumePassTime;

public:

    /**
     * @brief getSolidvoxelizationAxis
     * @param argumentString
     * @return
     */
    static SOLID_VOXELIZATION_AXIS getSolidvoxelizationAxis(const std::string &argumentString);
};

}

#endif // ULTRALISER_DATA_VOLUME_H
