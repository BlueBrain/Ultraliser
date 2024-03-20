/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

#include <common/Common.h>
#include <data/common/GridIndex.h>
#include <data/volumes/grids/BitVolumeGrid.h>
#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/morphologies/Morphologies.h>
#include <data/meshes/simple/Mesh.h>
#include <data/volumes/utilities/VolumeType.hh>
#include <data/volumes/utilities/VolumeData.hh>
#include <algorithms/skeletonization/thinning/Thinning6Iterations.h>
#include <data/volumes/voxels/CandiateVoxel.h>
#include <data/fields/Field.h>

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

    /**
     * @brief The SOLID_VOXELIZATION_AXIS enum
     */
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
     * Copy constructor
     * @param volume Input volume to be copied.
     */
    Volume(const Volume* volume);

    /**
     * @brief Volume
     * @param filePath
     */
    Volume(const std::string &filePath);

    /**
     * @brief Volume
     * Constructor
     *
     * @param pMin
     *
     * @param pMax

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
           const size_t &baseResolution = 512,
           const float &expansionRatio = 0.0,
           const VOLUME_TYPE& gridType = VOLUME_TYPE::BIT,
           const bool verbose = VERBOSE);

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
           const Vector3f pMin = Vector3f::ZERO,
           const Vector3f pMax = Vector3f::ZERO,
           const VOLUME_TYPE& gridType = VOLUME_TYPE::BIT,
           const float expansionRatio = 0.0);

    ~Volume();

    /**
     * @brief getByte
     * @param index
     * @return
     */
    uint8_t getByte(const size_t index) const;

    /**
     * @brief getConfirmedValue
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getConfirmedValue(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI8
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getValueUI8(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI16
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint16_t getValueUI16(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI32
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint32_t getValueUI32(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueUI64
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint64_t getValueUI64(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueF32
     * @param x
     * @param y
     * @param z
     * @return
     */
    float getValueF32(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief getValueF64
     * @param x
     * @param y
     * @param z
     * @return
     */
    double getValueF64(const int64_t &x, const int64_t &y, const int64_t &z) const;


    /**
     * @brief fillVoxel
     * @param x
     * @param y
     * @param z
     */
    void fillVoxel(const int64_t &x, const int64_t &y, const int64_t &z);

    /**
     * @brief getByte
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getByte(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief addByte
     * @param index
     * @param byte
     */
    void addByte(const size_t &index, const uint8_t byte);

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
                             const bool parallel = false,
                             const float&sideRatio = 1.0);

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
     * @brief surfaceVoxelizeNeuronMorphology
     * @param morphology
     */
    void surfaceVoxelizeNeuronMorphology(NeuronMorphology* neuronMorphology,
                                         const std::string &packingAlgorithm,
                                         const bool& verbose = SILENT);

    /**
     * @brief surfaceVoxelizeSpineMorphology
     * @param spineMorphology
     * @param packingAlgorithm
     */
    void surfaceVoxelizeSpineMorphology(SpineMorphology* spineMorphology,
                                        const std::string &packingAlgorithm,
                                        const bool verbose = SILENT);

    /**
     * @brief surfaceVoxelizeAstrocyteMorphology
     * @param morphology
     */
    void surfaceVoxelizeAstrocyteMorphology(
            const AstrocyteMorphology *astrocyteMorphology, float threshold = 0.75,
            const std::string &packingAlgorithm = POLYLINE_PACKING);

    /**
     * @brief surfaceVoxelizeVasculatureMorphology
     * Create the volumetric shell of the vasculature morphology.
     * @param vasculatureMorphology
     * The input vascular morphology.
     * @param packingAlgorithm
     * The used packing algorithm.
     */
    void surfaceVoxelizeVasculatureMorphology(
            VasculatureMorphology* vasculatureMorphology,
            const std::string& packingAlgorithm = POLYLINE_PACKING);

    /**
     * @brief solidVoxelization
     * @param axis
     * @param verbose
     */
    void solidVoxelization(const SOLID_VOXELIZATION_AXIS& axis=X, const bool &verbose = VERBOSE);

    /**
     * @brief solidVoxelizationROI
     * @param axis
     * The voxelization axis.
     * @param x1
     * The lower bound along the X-axis.
     * @param x2
     * The upper bound along the X-axis.
     * @param y1
     * The lower bound along the Y-axis.
     * @param y2
     * The upper bound along the Y-axis.
     * @param z1
     * The lower bound along the Z-axis.
     * @param z2
     * The upper bound along the Z-axis.
     * @param verbose
     */
    void solidVoxelizationROI(const SOLID_VOXELIZATION_AXIS& axis,
                                      const size_t& x1, const size_t x2,
                                      const size_t& y1, const size_t y2,
                                      const size_t& z1, const size_t z2,
                                      const bool& verbose = true);

    /**
     * @brief exportToMesh
     * Exports the volume into a mesh file in several file formats.
     * @param prefix
     * @param formatOBJ
     * @param formatPLY
     * @param formatOFF
     * @param formatSTL
     */
    void exportToMesh(const std::string &prefix,
                      const bool &formatOBJ = false,
                      const bool &formatPLY = false,
                      const bool &formatOFF = false,
                      const bool &formatSTL = false) const;

    /**
     * @brief exportVolumeGridMesh
     * Export the volumegrid into a mesh file in several file formats.
     * @param prefix
     * @param formatOBJ
     * @param formatPLY
     * @param formatOFF
     * @param formatSTL
     */
    void exportVolumeGridToMesh(const std::string &prefix,
                                const bool &formatOBJ = false,
                                const bool &formatPLY = false,
                                const bool &formatOFF = false,
                                const bool &formatSTL = false) const;
    /**
     * @brief exportBoundingBoxMesh
     * Exoprts the bounding box of the mesh in several file formats.
     * @param prefix
     * @param formatOBJ
     * @param formatPLY
     * @param formatOFF
     * @param formatSTL
     */
    void exportBoundingBoxMesh(const std::string &prefix,
                               const bool &formatOBJ = false,
                               const bool &formatPLY = false,
                               const bool &formatOFF = false,
                               const bool &formatSTL = false) const;

    /**
     * @brief writeStackXY
     * @param prefix
     */
    void writeStackXY(const std::string &outputDirectory, const std::string &prefix) const;

    /**
     * @brief writeStackXZ
     * @param outputDirectory
     * @param prefix
     */
    void writeStackXZ(const std::string &outputDirectory, const std::string &prefix) const;

    /**
     * @brief writeStackZY
     * @param outputDirectory
     * @param prefix
     */
    void writeStackZY(const std::string &outputDirectory, const std::string &prefix) const;

    /**
     * @brief writeVolumes
     * Write the volume to the file system using different file formats.
     * @param prefix
     * File prefix.
     * @param binaryFormat
     * The volume will be stored in .HDR/.BIN files, where each voxel will be stored in a
     * single bit.
     * @param rawFormat
     * The volume will be stored in .HDR/.IMG files, where each voxel will be stored
     * @param nrrdFormat
     * @param ultraliserFormat
     */
    void writeVolumes(const std::string &prefix,
                      const bool& bitFormat = false,
                      const bool& unsignedFormat = false,
                      const bool& floatFormat = false,
                      const bool& nrrdFormat = false,
                      const bool &rawFormat = false) const;

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
    size_t getNumberVoxels(void) const;

    /**
     * @brief computeNumberNonZeroVoxels
     * @return
     */
    size_t computeNumberNonZeroVoxels(void) const;

    /**
     * @brief computeVolume
     * @return
     */
    float computeVolume() const;

    /**
     * @brief getNumberBytes
     * @return
     */
    size_t getNumberBytes(void) const;

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
     * @brief mapTo1DIndexWithoutBoundCheck
     * @param x
     * @param y
     * @param z
     * @return
     */
    size_t mapTo1DIndexWithoutBoundCheck(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief mapToIndex
     * @param x
     * @param y
     * @param z
     * @return
     */
    size_t mapToIndex(const int64_t &x, const int64_t &y, const int64_t &z, bool &outlier) const;

    /**
     * @brief searchForBorderVoxels
     * @param verbose
     * @return
     */
    std::vector< std::vector< Vec3ui_64> > searchForBorderVoxels(const bool verbose = VERBOSE) ;

    /**
     * @brief searchForCandidateVoxels
     * @return
     */
    std::vector< CandidateVoxels > searchForCandidateVoxels() ;

    /**
     * @brief searchForCandidateVoxelsOne
     * @return
     */
    CandidateVoxels searchForCandidateVoxelsOne();


    /**
     * @brief deleteCandidateVoxels
     * @param thinning
     * @return
     */
    size_t deleteCandidateVoxels(std::unique_ptr< Thinning6Iterations > &thinning,
                                 const bool& displayProgress = true);

    size_t deleteCandidateVoxelsParallel(std::unique_ptr< Thinning6Iterations > &thinning);

    size_t deleteBorderVoxelsUsingThinningVoxels(std::unique_ptr< Thinning6Iterations > &thinning,
            ThinningVoxelsUI16List &thinningVoxels);

    /**
     * @brief searchForDeletableVoxels
     * @param perSliceBorderVoxels
     * @param thinning
     * @param direction
     * @return
     */
    std::vector< Vec3ui_64 > searchForDeletableVoxels(std::vector< std::vector< Vec3ui_64> > &perSliceBorderVoxels,
            std::unique_ptr<Thinning6Iterations> &thinning,
            int direction) const;

    /**
     * @brief confirmDeletableVoxels
     * @param candidateVoxels
     * @param thinning
     * @param direction
     */
    void confirmDeletableVoxels(CandidateVoxels& candidateVoxels,
                                std::unique_ptr< Thinning6Iterations > &thinning,
                                int direction) const;

    /**
     * @brief mapToXYZ
     * @param index
     * @param x
     * @param y
     * @param z
     */
    void mapToXYZ(const size_t index, size_t &x, size_t &y, size_t &z) const;

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
     * @brief isFilledWithoutBoundCheck
     * @param x
     * @param y
     * @param z
     * @return
     */
    bool isFilledWithoutBoundCheck(const int64_t &x, const int64_t &y, const int64_t &z) const;

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

    /**
     * @brief verifyBorderVoxels
     * @param mesh
     * @param sideRatio
     * @param verbose
     * @return
     */
    std::vector< Vec3ui_64 > verifyBorderVoxels(Mesh* mesh,
                                                const float& sideRatio, const bool& verbose);

    /**
     * @brief getBrick
     * @param x1
     * @param x2
     * @param y1
     * @param y2
     * @param z1
     * @param z2
     * @return
     */
    Volume* getBrick(const size_t& x1, const size_t& x2,
                     const size_t& y1, const size_t& y2,
                     const size_t& z1, const size_t& z2,
                     const size_t emptyShellThickness = 1);

    void addBrickToVolume(const Volume* brick,
                                  const size_t& xMin, const size_t& xMax,
                                  const size_t& yMin, const size_t& yMax,
                                  const size_t& zMin, const size_t& zMax,
                                  const size_t emptyShellThickness);

    void andBrickToVolume(const Volume* brick,
                          const size_t& xMin, const size_t& xMax,
                          const size_t& yMin, const size_t& yMax,
                          const size_t& zMin, const size_t& zMax,
                          const size_t emptyShellThickness);



    /**
     * @brief extractBrickFromVolume
     * Extracts a subset of the volume represented by a volume brick, which has the same data
     * structure of the volume.
     * @param xVolumeStart
     * The starting x-axis index of the volume, where the brick will be copied.
     * @param xVolumeEnd
     * The ending x-axis index of the volume, where the brick will be copied.
     * @param yVolumeStart
     * The starting y-axis index of the volume, where the brick will be copied.
     * @param yVolumeEnd
     * The ending y-axis index of the volume, where the brick will be copied.
     * @param zVolumeStart
     * The starting z-axis index of the volume, where the brick will be copied.
     * @param zVolumeEnd
     * The ending z-axis index of the volume, where the brick will be copied.
     * @return The extracted volume brick.
     */
    Volume* extractBrickFromVolume(const size_t& xVolumeStart, const size_t& xVolumeEnd,
                                   const size_t& yVolumeStart, const size_t& yVolumeEnd,
                                   const size_t& zVolumeStart, const size_t& zVolumeEnd,
                                   const bool& displayProgress = true);



    Volume* extractBoundedBrickFromVolume(const size_t& xVolumeStart, const size_t& xVolumeEnd,
                                          const size_t& yVolumeStart, const size_t& yVolumeEnd,
                                          const size_t& zVolumeStart, const size_t& zVolumeEnd,
                                          const size_t& numberBoundaryVoxels = 2,
                                          const bool& displayProgress = true) const;

    /**
     * @brief insertBrickToVolume
     * Copies (or overwrites, inserts) to content of a given brick to the volume. The spatial
     * extent of which the brick is going to be copied in the volume is defined by the given
     * volume indices.
     * @param brick
     * A given brick to be copied to the volume.
     * @param xVolumeStart
     * The starting x-axis index of the volume, where the brick will be copied.
     * @param xVolumeEnd
     * The ending x-axis index of the volume, where the brick will be copied.
     * @param yVolumeStart
     * The starting y-axis index of the volume, where the brick will be copied.
     * @param yVolumeEnd
     * The ending y-axis index of the volume, where the brick will be copied.
     * @param zVolumeStart
     * The starting z-axis index of the volume, where the brick will be copied.
     * @param zVolumeEnd
     * The ending z-axis index of the volume, where the brick will be copied.
     * @return True or false depending on the success of the operation.
     */
    bool insertBrickToVolume(const Volume* brick,
                             const size_t& xVolumeStart, const size_t& xVolumeEnd,
                             const size_t& yVolumeStart, const size_t& yVolumeEnd,
                             const size_t& zVolumeStart, const size_t& zVolumeEnd,
                             const bool &displayProgress = true);

    /**
     * @brief insertBoundedBrickToVolume
     * @param brick
     * @param xVolumeStart
     * @param xVolumeEnd
     * @param yVolumeStart
     * @param yVolumeEnd
     * @param zVolumeStart
     * @param zVolumeEnd
     * @param numberBoundaryVoxels
     * @return
     */
    bool insertBoundedBrickToVolume(const Volume* brick,
                                    const size_t& xVolumeStart, const size_t& xVolumeEnd,
                                    const size_t& yVolumeStart, const size_t& yVolumeEnd,
                                    const size_t& zVolumeStart, const size_t& zVolumeEnd,
                                    const size_t& numberBoundaryVoxels,
                                    const bool &displayProgress = true);

    bool insertOverlappingBoundedBrickToVolume(const Volume* brick,
                                               const size_t& xVolumeStart, const size_t& xVolumeEnd,
                                               const size_t& yVolumeStart, const size_t& yVolumeEnd,
                                               const size_t& zVolumeStart, const size_t& zVolumeEnd,
                                               const size_t& xOverlappingVoxels,
                                               const size_t &yOverlappingVoxels,
                                               const size_t &zOverlappingVoxels,
                                               const size_t& numberBoundaryVoxels,
                                               const bool &displayProgress = true);

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
                                         const VOLUME_TYPE &gridType);

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
                                           const size_t &isoValue);


    /**
     * @brief constructVolumeWithMinimumIsoValue
     * @param volume
     * @param minIsoValue
     * @return
     */
    static Volume* constructVolumeWithMinimumIsoValue(const Volume* volume,
                                                      const size_t& minIsoValue);

    /**
     * @brief constructVolumeWithMaximumIsoValue
     * @param volume
     * @param maxIsoValue
     * @return
     */
    static Volume* constructVolumeWithMaximumIsoValue(const Volume* volume,
                                                      const size_t& maxIsoValue);

    /**
     * @brief constructVolumeWithIsoRange
     * @param volume
     * @param minIsoValue
     * @param maxIsoValue
     * @return
     */
    static Volume* constructVolumeWithIsoRange(const Volume* volume,
                                               const size_t& minIsoValue,
                                               const size_t& maxIsoValue);



    /**
     * @brief constructNonZeroVolume
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
    static Volume* constructNonZeroVolume(const Volume* volume);


    /**
     * @brief constructIsoValuesVolume
     * @param volume
     * @param isoValues
     * @return
     */
    static Volume* constructIsoValuesVolume(const Volume* volume,
                                            const std::vector< size_t > &isoValues);

    /**
     * @brief createHistogram
     * Creates the histogram of the given volume.
     * @param volume
     * A given byte volume.
     * @return
     * Histogram array.
     */
    static std::vector<size_t> createHistogram(const Volume* volume, const VOLUME_TYPE &type);

    /**
     * @brief getVoxelBoundingBox
     * Gets the bounding box of a specific voxel.
     * @param x
     * X-index in the volume.
     * @param y
     * Y-index in the volume.
     * @param z
     * Z-index in the volume.
     * @param pMin
     * Returned pMin of the voxel.
     * @param pMax
     * Returned pMax of the voxel.
     */
    void getVoxelBoundingBox(const int64_t& x, const int64_t& y, const int64_t& z,
                             Vector3f& pMin, Vector3f& pMax) const;

    VolumeGrid* getGrid() const { return _grid; }

    VOLUME_TYPE getVolumeType() const { return _gridType; }

    /**
     * @brief getVolumeBoundingBox
     * Gets the boundin box of the volume.
     * @param pMin
     * Returned pMin of the volume.
     * @param pMax
     * Returned pMax of the volume.
     */
    void getVolumeBoundingBox(Vector3f& pMin, Vector3f& pMax) const;

    /**
     * @brief getVoxelSize
     * Gets the size of the voxel.
     * @return
     * Returns the size of the voxel.
     */
    float getVoxelSize() const { return _voxelSize; }

    /**
     * @brief getPMin
     * @return
     */
    Vector3f getPMin() const { return _pMin; }

    /**
     * @brief getPMax
     * @return
     */
    Vector3f getPMax() const { return _pMax; }

    /**
     * @brief getCenter
     * @return
     */
    Vector3f getCenter() const { return _center; }

    /**
     * @brief getScale
     * @return
     */
    Vector3f getScale() const { return _scale; }

    /**
     * @brief getType
     * @return
     */
    VOLUME_TYPE getType() const { return _gridType; }

    /**
     * @brief isBorderVoxel
     * @param x
     * @param y
     * @param z
     * @return
     */
    bool isBorderVoxel(const int64_t& x, const int64_t& y,const int64_t& z) const;

    /**
     * @brief surfaceVoxelizeSections
     * @param section
     * @param verbose
     */
    void surfaceVoxelizeSections(const Sections &section, const bool& verbose = VERBOSE);



    /**
     * @brief surfaceVoxelizationRegion
     * @param mesh
     * @param pMinRegion
     * @param pMaxRegion
     * @param verbose
     */
    void surfaceVoxelizationRegion(Mesh* mesh,
                                   const Vector3f& pMinRegion,
                                   const Vector3f& pMaxRegion,
                                   const bool& verbose = false);


    float getExpansionRatio() const { return _expansionRatio; }

    Bounds3D_ui64 getROIBounds(const Vector3f& pMin, const Vector3f& pMax);

    Bounds3D_ui64 getActiveRegionBounds();


    void projectXY(const std::string& prefix, const bool &projectColorCoded = false);
    void projectYZ(const std::string& prefix, const bool &projectColorCoded = false) const;
    void projectXZ(const std::string& prefix, const bool &projectColorCoded = false) const;

    ThinningVoxelsUI16List& getThinningVoxelsList(const bool &rebuild = false,
                                                  const bool verbose = VERBOSE)
    {
        return _grid->getThinningVoxelsList(rebuild, verbose);
    }

    size_t getBaseResolution() const
    {
        return _baseResolution;
    }

private:

    void _buildOccupancyRanges();

    void _buildVolumeOccupancy();

    /**
     * @brief _allocateGrid
     * @param width
     * @param height
     * @param depth
     */
    void _allocateGrid(const size_t& width, const size_t& height, const size_t& depth);

    /**
     * @brief _createGrid
     * Creates the volume grid and initializes it to zeros.
     */
    void _createGrid(const bool verbose = VERBOSE);

    /**
     * @brief _createGrid
     * @param headerFilePath
     */
    void _createGrid(const std::string& headerFilePath);

    /**
     * @brief _createGrid
     * @param volumeData
     */
    void _createGrid(const UltraliserVolumeData* volumeData);

#ifdef ULTRALISER_USE_NRRD
    /**
     * @brief _createGrid
     * @param volumeData
     */
    void _createGrid(const NRRDVolumeData* volumeData);
#endif

    /**
     * @brief _clampIndex
     * @param idx
     * @param dimension
     * @return
     */
    size_t _clampIndex(size_t idx, size_t dimension);

    /**
     * @brief _getBoundingBox
     * @param mesh
     * @param i
     * @param tMin
     * @param tMax
     */
    void _getBoundingBox(Mesh* mesh, size_t i, int64_t *tMin, int64_t *tMax);

    /**
     * @brief _getBoundingBox
     * @param sample
     * @param tMin
     * @param tMax
     */
    void _getBoundingBox(Sample* sample, int64_t *tMin, int64_t *tMax);

    /**
     * @brief _getBoundingBox
     * @param sample0
     * @param sample1
     * @param tMin
     * @param tMax
     */
    void _getBoundingBox(Sample* sample0, Sample* sample1, int64_t *tMin, int64_t *tMax);

    /**
     * @brief _getTriangleBoundingBox
     * @param triangle
     * @param tMin
     * @param tMax
     */
    void _getTriangleBoundingBox(AdvancedTriangle triangle, int64_t *tMin, int64_t *tMax);

    /**
     * @brief _rasterize
     * @param mesh
     * @param grid
     */
    void _rasterize(Mesh* mesh, VolumeGrid* grid,
                    const float &sideRatio = 1.0, const bool& verbose = false);

    /**
     * @brief _rasterizeRegion
     * @param mesh
     * @param grid
     * @param pMinRegion
     * @param pMaxRegion
     * @param verbose
     */
    void _rasterizeRegion(Mesh* mesh, VolumeGrid* grid,
                          const Vector3f& pMinRegion,
                          const Vector3f& pMaxRegion,
                          const bool& verbose = false);

    /**
     * @brief _rasterize
     * @param sample
     * @param grid
     */
    void _rasterize(Sample* sample, VolumeGrid* grid);

    /**
     * @brief _rasterize
     * @param sample0
     * @param sample1
     * @param grid
     * @param stepAlpha
     */
    void _rasterize(Sample* sample0, Sample* sample1, VolumeGrid* grid, float stepAlpha = 0.1f);

    /**
     * @brief _rasterizeParallel
     * @param mesh
     * @param grid
     * @param verbose
     */
    void _rasterizeParallel(Mesh* mesh, VolumeGrid* grid, const float &sideRatio = 1.0);

    /**
     * @brief _floodFill2D
     * Performs a 2D flood-filling operation.
     * @param axis
     * @param verbose
     */
    void _floodFill2D(const SOLID_VOXELIZATION_AXIS& axis, const bool &verbose = true);

    /**
     * @brief _floodFill2DROI
     * @param axis
     * @param x1
     * @param x2
     * @param y1
     * @param y2
     * @param z1
     * @param z2
     * @param verbose
     */
    void _floodFill2DROI(const SOLID_VOXELIZATION_AXIS& axis,
                         const size_t& x1, const size_t x2,
                         const size_t& y1, const size_t y2,
                         const size_t& z1, const size_t z2,
                         const bool &verbose = true);

    /**
     * @brief _floodFillAlongAxis
     * @param grid
     * @param axis
     * @param verbose
     */
    void _floodFillAlongAxis(VolumeGrid* grid, const SOLID_VOXELIZATION_AXIS &axis,
                             const bool &verbose = true);

    void _floodFillAlongAxisROI(VolumeGrid* grid, const SOLID_VOXELIZATION_AXIS &axis,
                                const size_t& x1, const size_t x2,
                                const size_t& y1, const size_t y2,
                                const size_t& z1, const size_t z2,
                                const bool &verbose);

    /**
     * @brief _floodFillAlongXYZ
     * @param grid
     * @param verbose
     */
    void _floodFillAlongXYZ(VolumeGrid *grid, const bool &verbose = true);

    /**
     * @brief _floodFillAlongXYZROI
     * @param grid
     * @param x1
     * @param x2
     * @param y1
     * @param y2
     * @param z1
     * @param z2
     * @param verbose
     */
    void _floodFillAlongXYZROI(VolumeGrid *grid,
                               const size_t& x1, const size_t x2,
                               const size_t& y1, const size_t y2,
                               const size_t& z1, const size_t z2,
                               const bool &verbose = true);

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
                                       size_t tIdx,
                                       const GridIndex& voxel, const float& sideRatio=1.0);

    /**
     * @brief _testSampleCubeIntersection
     * @param sample
     * @param voxel
     * @return
     */
    bool _testSampleCubeIntersection(Sample* sample, const GridIndex& voxel);

    /**
     * @brief _testTriangleGridIntersection
     * @param triangle
     * @param voxel
     * @return
     */
    bool _testTriangleGridIntersection(AdvancedTriangle triangle, const GridIndex& voxel);

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
    int _triangleCubeSign(Mesh* mesh, int tIdx, const GridIndex& boundingBox);

private:

    /**
     * @brief grid
     */
    VolumeGrid* _grid;

    /**
     * @brief _gridType
     */
    VOLUME_TYPE _gridType;

    /**
     * @brief _pMin
     */
    Vector3f _pMin;

    /**
     * @brief _pMax
     */
    Vector3f _pMax;

    /**
     * @brief _center
     */
    Vector3f _center;

    /**
     * @brief _scale
     * The diagonal, i.e. the difference between _pMax and _pMin.
     */
    Vector3f _scale;

    /**
     * @brief _expansionRatio
     * Additional layer of voxels around the bounding box of the object.
     */
    float _expansionRatio;

    /**
     * @brief rawFileName
     */
    std::string _rawFileName;

    /**
     * @brief _baseResolution
     * Base resolution of the volume.
     */
    size_t _baseResolution;

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
    // Vector3f _meshOrigin;

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
