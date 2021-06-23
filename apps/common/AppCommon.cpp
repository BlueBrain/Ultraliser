/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 **************************************************************************************************/

#include <Ultraliser.h>
#include <AppOptions.h>

namespace Ultraliser
{

Volume* createVolumeGrid(Mesh *mesh, const AppOptions* options)
{
    // Get relaxed bounding box to build the volume
    Ultraliser::Vector3f pMinInput, pMaxInput;
    mesh->computeBoundingBox(pMinInput, pMaxInput);
    Vector3f meshBoundingBox = pMaxInput - pMinInput;

    // Get the largest dimension
    const float largestDimension = meshBoundingBox.getLargestDimension();

    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_SUCCESS("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    return new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                      Ultraliser::VolumeGrid::getType(options->volumeType));
}

void computeBoundingBoxForMeshes(const std::string& boundsFile,
                                 const std::string& inputMeshesDirectory,
                                 std::vector< std::string > meshFiles,
                                 Vector3f& pMax, Vector3f& pMin)
{
    if (boundsFile == EMPTY)
    {
        LOG_STATUS("Computing Bounding Box");
        TIMER_SET;

        // Vectors containing all the pMin and pMax of all the objects
        std::vector< Ultraliser::Vector3f > pMinVector, pMaxVector;

        // Resize to perform in parallel
        pMinVector.resize(meshFiles.size());
        pMaxVector.resize(meshFiles.size());

        LOOP_STARTS("Loading Meshes");
        size_t loadedMeshCount = 0;
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for (uint64_t iMesh = 0; iMesh < meshFiles.size(); ++iMesh)
        {
            // Create and load the mesh from the file
            std::string meshName= meshFiles[iMesh];
            std::string meshFile = meshName;
            if (inputMeshesDirectory != EMPTY )
                meshFile = inputMeshesDirectory + "/" + meshName;

            if (Ultraliser::File::exists(meshFile))
            {
                // Load the mesh
                auto mesh = new Ultraliser::Mesh(meshFile, false);

                // Compute its bounding box
                Ultraliser::Vector3f pMinMesh, pMaxMesh;
                mesh->computeBoundingBox(pMinMesh, pMaxMesh);
                pMinVector[iMesh] = pMinMesh;
                pMaxVector[iMesh] = pMaxMesh;
                mesh->~Mesh();
            }
            else
            {
                LOG_WARNING("Ignoring Mesh: [ %s ]", meshFile.c_str());
            }

            PROGRESS_UPDATE;
            LOOP_PROGRESS(PROGRESS,  meshFiles.size());
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);

        if (loadedMeshCount == 0 )
            LOG_ERROR("No Loaded Meshes");
        else
            LOG_DETAIL("Loaded Meshes: [%zu/%zu]", loadedMeshCount, meshFiles.size());


        // Compute the bounding box of the group
        pMax.x() = std::numeric_limits< float >::min();
        pMax.y() = std::numeric_limits< float >::min();
        pMax.z() = std::numeric_limits< float >::min();

        pMin.x() = std::numeric_limits< float >::max();
        pMin.z() = std::numeric_limits< float >::max();
        pMin.y() = std::numeric_limits< float >::max();

        LOOP_STARTS("Computing Bounding Box");
        TIMER_RESET;
        for(uint64_t iMesh = 0; iMesh < meshFiles.size(); ++iMesh)
        {
            LOOP_PROGRESS(iMesh,  meshFiles.size());

            Ultraliser::Vector3f pMinObject = pMinVector[ iMesh ];
            Ultraliser::Vector3f pMaxObject = pMaxVector[ iMesh ];

            if (pMinObject.x() < pMin.x()) pMin.x() = pMinObject.x();
            if (pMinObject.y() < pMin.y()) pMin.y() = pMinObject.y();
            if (pMinObject.z() < pMin.z()) pMin.z() = pMinObject.z();

            if (pMaxObject.x() > pMax.x()) pMax.x() = pMaxObject.x();
            if (pMaxObject.y() > pMax.y()) pMax.y() = pMaxObject.y();
            if (pMaxObject.z() > pMax.z()) pMax.z() = pMaxObject.z();
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOG_STATUS_IMPORTANT("Loading Bounding Box from [ %s ]", boundsFile.c_str());

        // Verify the bounding box file
        if (File::exists(boundsFile))
            File::parseBoundsFile(boundsFile, pMin, pMax);
        else
            LOG_ERROR("No Bounding Box File is Provided !");
    }
}

void applyLaplacianOperator(Mesh *mesh, const AppOptions* options)
{
    // Apply the Laplacian filter
    mesh->smoothLaplacian(options->laplacianIterations);

    // Export the mesh
    mesh->exportMesh(options->meshPrefix + LAPLACIAN_SUFFIX,
                     options->exportOBJ, options->exportPLY,
                     options->exportOFF, options->exportSTL);

    // Print the mesh statistcs
    if (options->writeStatistics)
        mesh->printStats(LAPLACIAN_STRING, &options->statisticsPrefix);

    // Print the mesh distributions
    if (options->writeDistributions)
        mesh->printStats(LAPLACIAN_STRING, &options->distributionsPrefix);
}

void createWatertightMesh(const Mesh* mesh, const AppOptions* options)
{
    // Create an advanced mesh to process the manifold mesh and make it watertight if it has any
    // self intersections
    std::unique_ptr< AdvancedMesh> watertightMesh = std::make_unique< AdvancedMesh >(
                mesh->getVertices(), mesh->getNumberVertices(),
                mesh->getTriangles(), mesh->getNumberTriangles());

    if (options->preservePartitions)
    {
        // Split the mesh into partitions
        std::vector < AdvancedMesh* > partitions = watertightMesh->splitPartitions();

        // Ensure watertightness for the rest of the partitions
        for (auto mesh : partitions)
        {
            try {
                mesh->ensureWatertightness();
            }  catch (...) {
                LOG_WARNING("Soma partition is invalid");
            }
        }

        // Ensures that the mesh is truly two-manifold with no self intersections
        watertightMesh->ensureWatertightness();

        // Merge back after checking the watertightness
        watertightMesh->appendMeshes(partitions);

        // Free
        for (auto mesh : partitions)
            delete mesh;
    }
    else
    {
        // Ensures that the mesh is truly two-advanced with no self intersections
        watertightMesh->ensureWatertightness();
    }

    // Print the mesh statistcs
    if (options->writeStatistics)
        watertightMesh->printStats(WATERTIGHT_STRING, &options->statisticsPrefix);

    // Print the mesh distributions
    if (options->writeDistributions)
        watertightMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

    // Export the repaired mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        watertightMesh->exportMesh(options->meshPrefix + WATERTIGHT_SUFFIX,
                                   options->exportOBJ, options->exportPLY,
                                   options->exportOFF, options->exportSTL);
}

void generateOptimizedMesh(Mesh *mesh, const AppOptions* options)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        mesh->optimizeAdaptively(options->optimizationIterations, options->smoothingIterations,
                                    options->flatFactor, options->denseFactor);

        mesh->smooth();
        mesh->smoothNormals();
    }
    else
    {
        // Default optimization
        if (options->optimizeMeshHomogenous)
            mesh->optimize(options->optimizationIterations,
                              options->smoothingIterations,
                              options->denseFactor);
    }

    if (!options->ignoreOptimizedNonWatertightMesh)
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            mesh->printStats(OPTIMIZED_STRING, &options->statisticsPrefix);

        // Print the mesh statistcs
        if (options->writeDistributions)
            mesh->writeDistributions(OPTIMIZED_STRING, &options->distributionsPrefix);

        // Export the mesh
        if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
            mesh->exportMesh(options->meshPrefix + OPTIMIZED_SUFFIX,
                                options->exportOBJ, options->exportPLY,
                                options->exportOFF, options->exportSTL);
    }

    // Fix self-intersections if any to create the watertight mesh
    if (!options->ignoreSelfIntersections)
        createWatertightMesh(mesh, options);
}

Volume* reconstructVolumeFromMesh(Mesh* inputMesh, const AppOptions* options,
                                  const bool& releaseInputMesh)
{
    auto volume = createVolumeGrid(inputMesh, options);

    // Surface voxelization
    volume->surfaceVoxelization(inputMesh, true, true);

    // Free the input mesh, if asked for
    if (releaseInputMesh)
        delete inputMesh;

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    return volume;
}

Mesh* reconstructMeshFromVolume(Volume* volume, const AppOptions* options)
{
    // Generate the reconstructed mesh from any of the marching cubes algorithms
    if (options->isosurfaceTechnique == DMC_STRING)
        return DualMarchingCubes::generateMeshFromVolume(volume, options->serialExecution);
    else
        return MarchingCubes::generateMeshFromVolume(volume, options->serialExecution);
}

void generateMarchingCubesMeshArtifacts(const Mesh *mesh, const AppOptions* options)
{   
    // Write the statistics of the reconstructed mesh from the marhcing cubes algorithm
    if (options->writeStatistics)
        mesh->printStats(MC_STRING, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        mesh->writeDistributions(MC_STRING, &options->distributionsPrefix);

    // Export the MC mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        mesh->exportMesh(options->meshPrefix + MC_SUFFIX,
                         options->exportOBJ, options->exportPLY,
                         options->exportOFF, options->exportSTL);
}


Mesh* loadInputMesh(const AppOptions* options)
{
    // Load the mesh and construct the mesh object, and generate its artifacts if needed
    auto mesh = new Mesh(options->inputMeshPath);

    // Write the statistics of the input mesh
    if (options->writeStatistics)
        mesh->printStats(INPUT_STRING, &options->statisticsPrefix);

    // Write the statistics of the input mesh
    if (options->writeDistributions)
        mesh->writeDistributions(INPUT_STRING, &options->statisticsPrefix);

    return mesh;
}

void generateReconstructedMeshArtifacts(Mesh* mesh, const AppOptions* options)
{
    // MC mesh output
    if (options->writeMarchingCubeMesh)
        generateMarchingCubesMeshArtifacts(mesh, options);

    // Laplacian smoorhing
    if (!options->ignoreLaplacianSmoothing || options->laplacianIterations > 0)
        applyLaplacianOperator(mesh, options);

    // Optimize the mesh and create a watertight mesh
    if (options->optimizeMeshHomogenous || options->optimizeMeshAdaptively)
        generateOptimizedMesh(mesh, options);
}

void generateVolumeArtifacts(const Volume* volume, const AppOptions* options)
{
    // Projecting the volume to validate its content
    if (options->projectXY || options->projectXZ || options->projectZY)
        volume->project(options->projectionPrefix,
                        options->projectXY, options->projectXZ, options->projectZY,
                        options->projectColorCoded);

    // Write the volume
    if (options->exportBitVolume || options->exportByteVolume || options->exportNRRDVolume)
        volume->writeVolumes(options->volumePrefix,
                             options->exportBitVolume,
                             options->exportByteVolume,
                             options->exportNRRDVolume);

    // Write the stacks
    if (options->exportStackXY || options->exportStackXZ || options->exportStackZY)
        volume->writeStacks(options->outputDirectory + "/" + STACKS_SIRECTORY, options->prefix,
                            options->exportStackXY, options->exportStackXZ, options->exportStackZY);

    // Export volume mesh
    if (options->exportVolumeMesh)
        volume->exportToMesh(options->meshPrefix,
                             options->exportOBJ, options->exportPLY,
                             options->exportOFF, options->exportSTL);

    // Export volume bounding box mesh
    if (options->exportVolumeBoundingBoxMesh)
        volume->exportBoundingBoxMesh(options->meshPrefix,
                                      options->exportOBJ, options->exportPLY,
                                      options->exportOFF, options->exportSTL);

    // Export volume grid mesh
    if (options->exportVolumeGridMesh)
        volume->exportVolumeGridToMesh(options->meshPrefix,
                                       options->exportOBJ, options->exportPLY,
                                       options->exportOFF, options->exportSTL);


    // Print the volume statistics
    if (options->writeStatistics)
        volume->printStats(options->prefix, &options->statisticsPrefix);
}

}
