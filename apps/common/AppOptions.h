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

#include <common/Common.h>
#include <arguments/Argument.h>
#include <data/volumes/Volume.h>
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
     * @brief verifyProcessingArguments
     * Verifies the processing arguments.
     */
    void verifyProcessingArguments();

    /**
     * @brief verifyInputMeshArgument
     * Verifies the input mesh argument.
     */
    void verifyInputMeshArgument();

    /**
     * @brief verifyTargetMeshArgument
     * Verifies the target mesh argument.
     */
    void verifyTargetMeshArgument();

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
     * @brief verifyVolumePrefixArgument
     * Verifies the prefix arguments for a given volume.
     */
    void verifyVolumePrefixArgument();

    /**
     * @brief verifyMaskDimensionsArguments
     * Verifies the dimensions of the given mask NOT to be zero.
     */
    void verifyMaskDimensionsArguments();

    /**
     * @brief verifyMorphologyExtractionArguments
     */
    void verifyMorphologyExtractionArguments();

    /**
     * @brief verifyPackingAlgorithmArgument
     * Verifies the option given to the packing algorithm
     */
    void verifyPackingAlgorithmArgument();

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
     * @brief verifyIsoOptionArgument
     */
    void verifyIsoOptionArgument();

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
     * @brief threads
     * The number of threads used to process the parallel chunks in the code.
     */
    size_t threads;

    /**
     * @brief inputMeshPath
     * An input mesh file path.
     */
    std::string inputMeshPath;

    /**
     * @brief targetMeshPath
     * The path to the target mesh where an input mesh will be mapped to.
     */
    std::string targetMeshPath;

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
     * @brief morphologyPrefix
     */
    std::string morphologyPrefix;

    /**
     * @brief meshPrefix
     */
    std::string meshPrefix;

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
    size_t maskWidth;

    /**
     * @brief maskHeight
     */
    size_t maskHeight;

    /**
     * @brief bboxWidth
     */
    float bboxWidth;

    /**
     * @brief bboxHeight
     */
    float bboxHeight;

    /**
     * @brief bboxDepth
     */
    float bboxDepth;

    /**
     * @brief bboxCenterX
     */
    float bboxCenterX;

    /**
     * @brief bboxCenterY
     */
    float bboxCenterY;

    /**
     * @brief bboxCenterZ
     */
    float bboxCenterZ;

    /**
     * @brief boundsFile
     * Use a bounds file to only voxelize part of the mesh or even a greater space. If the file
     * is not given, the bounding box of the input mesh will be used in addition to a little delta
     * to avoid intersection.
     */
    std::string boundsFile;

    /**
     * @brief isoOption
     * isoOption can be one of the following:
     * 1. isovalue: a single iso value is used to segment the volume (--isovalue).
     * 2. min: all the voxels with values over a given iso value will be selected (--min-isovalue)
     * 3. max: all the voxels with values below a given iso value will be selected (--max-isovalue)
     * 4. isovalues: a list of iso values are used to segment the volume (--iso-values-file)
     * 5. fullrange: all the non-zero voxels of the volumes will be used for the meshing (full-range)
     */
    std::string isoOption;

    /**
     * @brief isoValue
     * The isovalue where the volume will get segmented, default 127.
     */
    size_t isoValue;

    /**
     * @brief isovaluesFile
     * A file containing a list of isovalues used to construct the segmented volume.
     */
    std::string isovaluesFile;

    /**
     * @brief minIsoValue
     * Select all the values that are greater than or equal this isovalue.
     */
    size_t minIsoValue;

    /**
     * @brief maxIsoValue
     * Select all the values that are lower than this isovalue.
     */
    size_t maxIsoValue;

    /**
     * @brief nonZeroVoxels
     * If the voxel contains any value except it, then use it.
     */
    bool nonZeroVoxels;

    /**
     * @brief isosurfaceTechnique
     * The technique that is used to extract the iso surface from the mesh.
     * Either dmc (Dual Marching Cubes) or mc (Marching Cubes)
     */
    std::string isosurfaceTechnique;

    /**
     * @brief volumeResolution
     * The base resolution of the volume that corresponds to the largest dimension.
     */
    size_t volumeResolution;

    /**
     * @brief scaledResolution
     * Sets the resolution of the volume based on mesh dimensions.
     */
    bool scaledResolution = false;

    /**
     * @brief voxelsPerMicron
     * Number of voxels per micron in case of auto resolution.
     */
    float voxelsPerMicron;

    /**
     * @brief edgeGap
     */
    float edgeGap;

    /**
     * @brief zeroPaddingVoxels
     * The number of zero-padding voxels that will be appended to the volume to avoid any clipping
     * artifacts, default 0.
     */
    size_t zeroPaddingVoxels;

    /**
     * @brief solid
     * Use solid voxelization to fill the volume.
     */
    bool useSolidVoxelization = false;

    /**
     * @brief packingAlgorithm
     * The packing algorithm used to create the proxy mesh.
     */
    std::string packingAlgorithm;

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
     * @brief exportMorphology
     * Export the re-sampled morphologies in *Morpho2Mesh applications.
     */
    bool exportMorphology= false;

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
     * @brief exportUltraBitVolume
     */
    bool exportUltraBitVolume = false;

    /**
     * @brief createByteVolume
     * If this flag is set, a default raw volume will be created. This volume
     * has 1 byte per voxel.
     */
    bool exportRawVolume = false;

    /**
     * @brief exportUnsignedVolume
     */
    bool exportUnsignedVolume = false;

    /**
     * @brief exportFloatVolume
     */
    bool exportFloatVolume = false;

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
     * @brief minDihedralAngle
     * The minimum value of the dihedral angle, by default it should be 0.1 for the mesh to
     * be watertight.
     */
    float minDihedralAngle;

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
     * @brief autoParameters
     * Use the best parameters computed from the statistical analysis of the data.
     */
    bool autoParameters;
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
     * @brief ignoreMarchingCubesMesh
     * If this flag is set, the mesh reconstructed with the marching cubes algorithm will not
     * be written to disk.
     */
    bool ignoreMarchingCubesMesh = false;

    /**
     * @brief ignoreLaplacianMesh
     * If this flag is set, the mesh resulting from the application of the Laplacian operator will
     * be ignored and will not be written to disk.
     */
    bool ignoreLaplacianMesh = false;

    /**
     * @brief ignoreOptimizedMesh
     * If this flag is set, the optimized mesh will not be written to disk.
     */
    bool ignoreOptimizedMesh = false;

    /**
     * @brief ignoreWatertightMesh
     * If this flag is set, the watertight mesh will not be written to disk.
     */
    bool ignoreWatertightMesh = false;

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

    /**
     * @brief minSampleRadius
     * The radius of the smallest sample in the morphology.
     */
    float minSampleRadius = 0.05;

    /** 
     * @brief axonBranchOrder
     * Maximum branch order applies in the neuron morphology axon reconstruction 
     */
    size_t axonBranchOrder = INT_MAX;

    /** 
     * @brief basalBranchOrder
     * Maximum branch order applies in the neuron morphology basal dendrites reconstruction 
     */
    size_t basalBranchOrder = INT_MAX;

    /** 
     * @brief apicalBranchOrder
     * Maximum branch order applies in the neuron morphology apical dendrites reconstruction 
     */
    size_t apicalBranchOrder = INT_MAX;

    /**
     * @brief exportAstrocyteAtOrigin
     * Exports the astrocyte mesh at the origin, where the soma is located.
     */
    bool exportAstrocyteAtOrigin = false;
};

}
