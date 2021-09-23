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

#ifndef ULTRALISER_SYSTEM_OPTIONS_HH
#define ULTRALISER_SYSTEM_OPTIONS_HH

#include <common/Common.h>
#include <arguments/Argument.h>
#include <data/volumes/Volumes.h>

namespace Ultraliser
{

/**
 * @brief The Options class
 */
class AppOptions
{

public:

    /**
     * @brief AppOptions
     */
    AppOptions() { /* EMPTY CONSTRUCTOR */};

public:

    /**
     * @brief verifyInputMeshArgument
     * Verifies the input mesh argument.
     */
    void verifyInputMeshArgument();

    /**
     * @brief verifyInputMeshesDirectoryArgument
     * Verifies the directory that contains a list of meshes.
     */
    void verifyInputMeshesDirectoryArgument();

    /**
     * @brief verifyInputMorphologyArgument
     * Verifies the input morphology argument.
     */
    void verifyInputMorphologyArgument();

    /**
     * @brief verifyInputMaskDirectoryArgument
     * Verifies the presence of the directory that contains the mask.
     */
    void verifyInputMaskDirectoryArgument();

    /**
     * @brief verifyOutputDirectoryArgument
     * Verifies the output arguments.
     */
    void verifyOutputDirectoryArgument();

    /**
     * @brief verifyBoudsFileArgument
     * Verifies the bounds file.
     */
    void verifyBoudsFileArgument();

    /**
     * @brief verifyMeshPrefixArgument
     * Verifies the prefix arguments for a given mesh.
     */
    void verifyMeshPrefixArgument();

    /**
     * @brief verifyMorphologyPrefixArgument
     * Verifies the prefix arguments for a given morphology.
     */
    void verifyMorphologyPrefixArgument();

    /**
     * @brief verifyMeshesPrefixArgument
     * Verifies the prefix arguments for a given directory with meshes.
     */
    void verifyMeshesPrefixArgument();

    /**
     * @brief verifyMaskPrefixArgument
     * Verifies the prefix arguments for a given mask.
     */
    void verifyMaskPrefixArgument();

    /**
     * @brief verifyMaskDimensionsArguments
     * Verifies the dimensions of the given mask NOT to be zero.
     */
    void verifyMaskDimensionsArguments();

    /**
     * @brief verifyMeshExportArguments
     * Verifies the mesh export options.
     */
    void verifyMeshExportArguments();

    /**
     * @brief verifyIsoSurfaceExtractionArgument
     * Verifies the isosurface extraction technique.
     */
    void verifyIsoSurfaceExtractionArgument();

    /**
     * @brief verifyVolumeExportArguments
     * Verifies the volume export options.
     */
    void verifyVolumeExportArguments();

    /**
     * @brief verifyMeshExportLogic
     * Runs a simple logic to ensure that there is at least ONE output mesh will be exported.
     */
    void verifyMeshExportLogic();

    /**
     * @brief createRespectiveDirectories
     * Create the directory tree where the artifacts will be genrated.
     */
    void createRespectiveDirectories();

    /**
     * @brief initializeContext
     * Initialization of the Ultralisation context.
     */
    void initializeContext();

public:

    /**
     * @brief serialExecution
     * Execute the workflows in a single thread for validation.
     */
    bool serialExecution = false;

    /**
     * @brief inputMeshPath
     * An input mesh file path.
     */
    std::string inputMeshPath;

    /**
     * @brief inputVolumePath
     * An input volume file path.
     */
    std::string inputVolumePath;

    /**
     * @brief inputMeshesDirectory
     * A directory that contains a lists of meshes that will be all loaded in Ultraliser at once.
     */
    std::string inputMeshesDirectory;

    /**
     * @brief inputMaskDirectory
     * The directory that contains a list of .tiff files that correspond to a segmented mask.
     */
    std::string inputMaskDirectory;

    /**
     * @brief inputMorphologyPath
     * An input morphology file, whether for neurons, astrocytes or vasculature.
     */
    std::string inputMorphologyPath;

    /**
     * @brief outputDirectory
     * The directory where the resulting data or artifacts will be created including meshes,
     * volumes, projections, stacks, etc.
     */
    std::string outputDirectory;

    /**
     * @brief projectionPrefix
     */
    std::string projectionPrefix;

    /**
     * @brief meshPrefix
     */
    std::string meshPrefix;

    /**
     * @brief morphologyPrefix
     */
    std::string morphologyPrefix;

    /**
     * @brief volumePrefix
     */
    std::string volumePrefix;

    /**
     * @brief statisticsPrefix
     */
    std::string statisticsPrefix;

    /**
     * @brief distributionsPrefix
     */
    std::string distributionsPrefix;

    /**
     * @brief maskWidth
     */
    uint64_t maskWidth;

    /**
     * @brief maskHeight
     */
    uint64_t maskHeight;

    /**
     * @brief boundsFile
     * Use a bounds file to only voxelize part of the mesh or even a greater space. If the file
     * is not given, the bounding box of the input mesh will be used in addition to a little delta
     * to avoid intersection.
     */
    std::string boundsFile;

    /**
     * @brief isoValue
     * The iso value where the volume will get segmented, default 127.
     */
    uint64_t isoValue;

    /**
     * @brief isosurfaceTechnique
     * The technique that is used to extract the iso surface from the mesh.
     * Either dmc (Dual Marching Cubes) or mc (Marching Cubes)
     */
    std::string isosurfaceTechnique;

    /**
     * @brief fullRangeIsoValue
     * If the voxel contains any value, then use it.
     */
    bool fullRangeIsoValue;

    /**
     * @brief volumeResolution
     * The base resolution of the volume that corresponds to the largest dimension.
     */
    uint64_t volumeResolution;

    /**
     * @brief autoResolution
     * Sets the resolution of the volume based on mesh dimensions.
     */
    bool autoResolution = false;

    /**
     * @brief voxelsPerMicron
     * Number of voxels per micron in case of auto resolution.
     */
    uint64_t voxelsPerMicron;

    /**
     * @brief edgeGap
     */
    float edgeGap;

    /**
     * @brief zeroPaddingVoxels
     * The number of zero-padding voxels that will be appended to the volume to avoid any clipping
     * artifacts, default 0.
     */
    uint64_t zeroPaddingVoxels;

    /**
     * @brief solid
     * Use solid voxelization to fill the volume.
     */
    bool useSolidVoxelization = false;

    /**
     * @brief voxelizationAxis
     * The axis where the solid voxelization operation will be performed.
     */
    Ultraliser::Volume::SOLID_VOXELIZATION_AXIS voxelizationAxis;

    /**
     * @brief volumeType
     * Use a specific volume for the voxelization process.
     */
    std::string volumeType;

    /**
     * @brief projectXY
     * If this flag is set, the XY projection of the volume will be saved to a PNG image.
     * This flag is set to validate the output volume.
     */
    bool projectXY = false;

    /**
     * @brief projectXZ
     * If this flag is set, the XY projection of the volume will be saved to a PNG image.
     * This flag is set to validate the output volume.
     */
    bool projectXZ = false;

    /**
     * @brief projectZY
     * If this flag is set, the ZY projection of the volume will be saved to a PNG image.
     * This flag is set to validate the output volume.
     */
    bool projectZY = false;

    /**
     * @brief projectColorCoded
     * If this flag is set, a series of color-coded projections with different color maps will
     * be generated.
     */
    bool projectColorCoded = false;

    /**
     * @brief exportStackXY
     * Create an image stack along the XY plane and export it.
     */
    bool exportStackXY = false;

    /**
     * @brief exportStackXZ
     * Create an image stack along the XZ plane and export it.
     */
    bool exportStackXZ = false;

    /**
     * @brief exportStackZY
     * Create an image stack along the ZY plane and export it.
     */
    bool exportStackZY = false;

    /**
     * @brief createBinaryVolume
     * If this flag is set, a binary volume will be created. This volume has
     * 1 bit per voxel.
     */
    bool exportBitVolume = false;

    /**
     * @brief createByteVolume
     * If this flag is set, a default raw volume will be created. This volume
     * has 1 byte per voxel.
     */
    bool exportByteVolume = false;

    /**
     * @brief exportNRRDVolume
     * If this flag is set, the volume will be written to an NRRD file that is
     * compatible with VTK.
     */
    bool exportNRRDVolume = false;

    /**
     * @brief exportVolumeMesh
     * Export a mesh that represents the volume where each voxel will be represented by a cube.
     */
    bool exportVolumeMesh = false;

    /**
     * @brief exportVolumeBoundingBoxMesh
     * Export a mesh that represents the bounding box of the volume.
     */
    bool exportVolumeBoundingBoxMesh = false;

    /**
     * @brief exportVolumeGridMesh
     * Export a mesh that represents the volumetric grid used to voxelize the mesh.
     */
    bool exportVolumeGridMesh = false;

    /**
     * @brief ignoreLaplacianSmoothing
     * Ignore Laplacian smoothing to smooth the reconstructed mesh.
     */
    bool ignoreLaplacianSmoothing = false;

    /**
     * @brief laplacianIterations
     * Number of iterations of the Laplacian smoothing filter.
     */
    int64_t laplacianIterations;

    /**
     * @brief optimizeMeshHomogenous
     * Optimize the reconstructed mesh using the default optimization strategy.
     */
    bool optimizeMeshHomogenous;

    /**
     * @brief optimizeMeshAdaptively
     * Optimize the mesh using the adaptive optimization strategy.
     */
    bool optimizeMeshAdaptively = false;

    /**
     * @brief voxelizeMesh
     * If this flag is set in some application, the input mesh will be voxelized.
     */
    bool voxelizeMesh = false;

    /**
     * @brief simpleMeshJoin
     * When merging a group of meshes together, if this flag is set, the meshes will be just
     * appended together without any volume reconstruction operation and therefore they will
     * not be optimized.
     */
    bool simpleMeshJoin = false;

    /**
     * @brief optimizationIterations
     * Number of iterations of optimizing the mesh surface.
     * By default it is set to 1, but more accurate numbers depend on the
     * given mesh with trial and error.
     */
    float optimizationIterations;

    /**
     * @brief smoothingIterations
     * Number of iterations required to smooth the optimized mesh, by default 10.
     */
    int64_t smoothingIterations;

    /**
     * @brief flatFactor
     * A factor that is used for the coarseFlat function.
     * Default value is 0.1.
     */
    float flatFactor;

    /**
     * @brief smoothingFactor
     */
    float smoothingFactor;

    /**
     * @brief denseFactor
     * A factor that is used for the coarseDense function.
     * Default value is 5.0.
     */
    float denseFactor;

    /**
     * @brief exportOBJ
     * Export any reconstructed mesh to .OBJ file.
     */
    bool exportOBJ = false;

    /**
     * @brief exportPLY
     * Export any reconstructed mesh to .PLY file.
     */
    bool exportPLY = false;

    /**
     * @brief exportOFF
     * Export any reconstructed mesh to .OFF file.
     */
    bool exportOFF = false;

    /**
     * @brief exportSTL
     * Export any reconstructed mesh to .STL file.
     */
    bool exportSTL = false;

    /**
     * @brief preservePartitions
     * Keeps all the mesh partitions in the optimized mesh.
     */
    bool preservePartitions;

    /**
     * @brief prefix
     * Just a prefix that will be used to label the output files. If this
     * is not given by the user, the name of the mesh file will be used.
     */
    std::string prefix;

    /**
     * @brief writeStatictics
     * Write the statictics.
     */
    bool writeStatistics = false;

    /**
     * @brief writeDistributions
     * Write distributions of morphologies and meshes.
     */
    bool writeDistributions = false;

    /**
     * @brief writeHistogram
     * Write the histogram of the volume into a text file.
     */
    bool writeHistogram = false;

    /**
     * @brief outputPrefix
     * Simply, the [OUTPUT_DIRECTORY]/[PREFIX]. This variable is just added to make the code simpler.
     */
    std::string outputPrefix;

    /**
     * @brief writeMarchingCubeMesh
     * By default the reconstructed mesh with the marching cubes algorithm is not written. But for
     * debugging purposes, we might need to write it.
     */
    bool writeMarchingCubeMesh = false;

    /**
     * @brief ignoreSelfIntersections
     * Ignore if the mesh has self intersections, and do NOT repair them.
     */
    bool ignoreSelfIntersections = false;

    /**
     * @brief ensureWatertight
     */
    bool watertight = false;

    /**
     * @brief ignoreOptimizedNonWatertightMesh
     * Ignores writing the optimized mesh that is not watertight.
     */
    bool ignoreOptimizedNonWatertightMesh = false;

    /**
     * @brief xScaleFactor
     * A scale factor along the X-axis.
     */
    float xScaleFactor = 1.f;

    /**
     * @brief yScaleFactor
     * A scale factor along the Y-axis.
     */
    float yScaleFactor = 1.f;

    /**
     * @brief zScaleFactor
     * A scale factor along the Z-axis.
     */
    float zScaleFactor = 1.f;
};

}

#endif // ULTRALISER_ARGUMENTS_ARGUMENTS_ULTRALISER_OPTIONS_HH

