#ifndef ULTRALISER_SYSTEM_OPTIONS_HH
#define ULTRALISER_SYSTEM_OPTIONS_HH
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

#include <common/Common.h>
#include <data/volumes/Volumes.h>

namespace Ultraliser
{

/**
 * @brief The Options struct
 */
struct Options
{
    /**
     * @brief inputMesh
     * An input mesh file.
     */
    std::string inputMesh;

    /**
     * @brief inputMorphology
     * An input morphology file, whether for neurons, astrocytes or vasculature.
     */
    std::string inputMorphology;

    /**
     * @brief outputDirectory
     * The directory where the resulting data or artifacts will be created including meshes,
     * volumes, projections, stacks, etc.
     */
    std::string outputDirectory;

    std::string projectionPrefix;
    std::string meshPrefix;
    std::string volumePrefix;
    std::string statisticsPrefix;
    std::string distributionsPrefix;


    /**
     * @brief boundsFile
     * Use a bounds file to only voxelize part of the mesh or even a greater space. If the file
     * is not given, the bounding box of the input mesh will be used in addition to a little delta
     * to avoid intersection.
     */
    std::string boundsFile;

    /**
     * @brief volumeResolution
     * The base resolution of the volume that corresponds to the largest dimension.
     */
    uint64_t volumeResolution;

    /**
     * @brief autoResolution
     * Sets the resolution of the volume based on mesh dimensions.
     */
    bool autoResolution;

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
     * @brief solid
     * Use solid voxelization to fill the volume.
     */
    bool useSolidVoxelization;

    /**
     * @brief solid
     * Fill the interior of the volume using solid voxelization.
     */
    Ultraliser::Volume::SOLID_VOXELIZATION_AXIS VoxelizationAxis;

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
    bool projectXY;

    /**
     * @brief projectXZ
     * If this flag is set, the XY projection of the volume will be saved to a PNG image.
     * This flag is set to validate the output volume.
     */
    bool projectXZ;

    /**
     * @brief projectZY
     * If this flag is set, the ZY projection of the volume will be saved to a PNG image.
     * This flag is set to validate the output volume.
     */
    bool projectZY;

    /**
     * @brief projectColorCoded
     * If this flag is set, a series of color-coded projections with different color maps will
     * be generated.
     */
    bool projectColorCoded;

    /**
     * @brief stackXY
     * Create an image stack along the XY plane.
     */
    bool stackXY;

    /**
     * @brief stackXZ
     * Create an image stack along the XZ plane.
     */
    bool stackXZ;

    /**
     * @brief stackZY
     * Create an image stack along the ZY plane.
     */
    bool stackZY;

    /**
     * @brief createBinaryVolume
     * If this flag is set, a binary volume will be created. This volume has
     * 1 bit per voxel.
     */
    bool writeBitVolume;

    /**
     * @brief createByteVolume
     * If this flag is set, a default raw volume will be created. This volume
     * has 1 byte per voxel.
     */
    bool writeByteVolume;

    /**
     * @brief writeNRRDVolume
     * If this flag is set, the volume will be written to an NRRD file that is
     * compatible with VTK.
     */
    bool writeNRRDVolume;

    /**
     * @brief exportVolumeMesh
     * Export a mesh that represents the volume where each voxel will be represented by a cube.
     */
    bool exportVolumeMesh;

    /**
     * @brief useLaplacian
     * Use Laplacian smoothing to clear the grid artifacts.
     */
    bool useLaplacian;

    /**
     * @brief laplacianIterations
     * Number of iterations of the Laplacian smoothing filter.
     */
    int64_t laplacianIterations;

    /**
     * @brief optimizeMesh
     * Optimize the reconstructed mesh using the default optimization strategy.
     */
    bool optimizeMesh;

    /**
     * @brief optimizeMeshAdaptively
     * Optimize the mesh using the adaptive optimization strategy.
     */
    bool optimizeMeshAdaptively;

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
     * @brief denseFactor
     * A factor that is used for the coarseDense function.
     * Default value is 5.0.
     */
    float denseFactor;

    /**
     * @brief exportOBJ
     * Export any reconstructed mesh to .OBJ file.
     */
    bool exportOBJ;

    /**
     * @brief exportPLY
     * Export any reconstructed mesh to .PLY file.
     */
    bool exportPLY;

    /**
     * @brief exportOFF
     * Export any reconstructed mesh to .OFF file.
     */
    bool exportOFF;

    /**
     * @brief exportSTL
     * Export any reconstructed mesh to .STL file.
     */
    bool exportSTL;

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
    bool writeStatistics;

    float smoothingFactor;



    /**
     * @brief writeDistributions
     * Write distributions of morphologies and meshes.
     */
    bool writeDistributions;

    /**
     * @brief outputPrefix
     * Simply, the [OUTPUT_DIRECTORY]/[PREFIX]. This variable is just added to make the code simpler.
     */
    std::string outputPrefix;

    /**
     * @brief ignoreDMCMesh
     * If this flag is set, the resulting mesh from the DMC stage will be
     * ignored and not exported to disk. Note that if this flag is set and
     * the ignoreSelfIntersection flag is set as well, there will no be
     * any exported mesh from this process.
     */
    bool ignoreDMCMesh;

    /**
     * @brief ignoreSelfIntersections
     * Ignore if the mesh has self intersections, and do NOT repair them.
     */
    bool ignoreSelfIntersections;

    /**
     * @brief ignoreOptimizedNonWatertightMesh
     * Ignores writing the optimized mesh that is not watertight.
     */
    bool ignoreOptimizedNonWatertightMesh;
};

}

#endif // ULTRALISER_ARGUMENTS_ARGUMENTS_ULTRALISER_OPTIONS_HH

