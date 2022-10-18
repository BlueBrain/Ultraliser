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

    size_t resolution;
    if (options->scaledResolution)
        resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);
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
                                 Vector3f& pMax, Vector3f& pMin,
                                 const float& xScale,const float& yScale,const float& zScale)
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
        for (size_t iMesh = 0; iMesh < meshFiles.size(); ++iMesh)
        {
            // Create and load the mesh from the file
            std::string meshName= meshFiles[iMesh];
            std::string meshFile = meshName;
            if (inputMeshesDirectory != EMPTY )
                meshFile = inputMeshesDirectory + "/" + meshName;

            if (Ultraliser::File::exists(meshFile))
            {
                // One more mesh is readable
                loadedMeshCount++;

                // Load the mesh
                auto mesh = std::make_unique< Ultraliser::Mesh >(meshFile, false);

                // Scale the mesh
                mesh->scale(xScale, yScale, zScale);

                // Compute its bounding box
                Ultraliser::Vector3f pMinMesh, pMaxMesh;
                mesh->computeBoundingBox(pMinMesh, pMaxMesh);
                pMinVector[iMesh] = pMinMesh;
                pMaxVector[iMesh] = pMaxMesh;
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
        pMax.x() = std::numeric_limits< float >::lowest();
        pMax.y() = std::numeric_limits< float >::lowest();
        pMax.z() = std::numeric_limits< float >::lowest();

        pMin.x() = std::numeric_limits< float >::max();
        pMin.z() = std::numeric_limits< float >::max();
        pMin.y() = std::numeric_limits< float >::max();

        LOOP_STARTS("Computing Bounding Box");
        TIMER_RESET;
        for(size_t iMesh = 0; iMesh < meshFiles.size(); ++iMesh)
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

void ensureWatertightness(Mesh* mesh, const AppOptions* options)
{
    std::unique_ptr< AdvancedMesh > watertightMesh =  std::make_unique< AdvancedMesh >(
                mesh->getVertices(), mesh->getNumberVertices(),
                mesh->getTriangles(), mesh->getNumberTriangles());

    // Release the data of the mesh to keep some free space for processing
    mesh->relaseData();

    if (options->preservePartitions)
    {
        // Split the mesh into partitions
        std::vector < AdvancedMesh* > partitions = watertightMesh->splitPartitions();

        // Ensure watertightness for the rest of the partitions
        for (auto partition : partitions)
        {
            try {
                partition->ensureWatertightness();
            } catch (...) {
                LOG_WARNING("Some partition is invalid");
            }
        }

        // Ensures that the mesh is truly two-manifold with no self intersections
        watertightMesh->ensureWatertightness();

        // Merge back after checking the watertightness
        watertightMesh->appendMeshes(partitions);

        // Free
        for (auto partition : partitions)
            delete partition;
    }
    else
    {
        // Ensures that the mesh is truly two-advanced with no self intersections
        watertightMesh->ensureWatertightness();
    }

    // Update the simple mesh data again
    watertightMesh->toSimpleMesh(mesh);
}

void applySmoothingOperator(Mesh *mesh, const AppOptions* options)
{   
    // Apply the smoothing filter
    mesh->smoothSurface(options->laplacianIterations);

    // TODO: We don't realy need it
    // ensureWatertightness(mesh, options);

    if (!options->ignoreLaplacianMesh)
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            mesh->printStats(LAPLACIAN_STRING, &options->statisticsPrefix);

        // Print the mesh distributions
        if (options->writeDistributions)
            mesh->printStats(LAPLACIAN_STRING, &options->distributionsPrefix);

        // Export the laplacian mesh
        mesh->exportMesh(options->meshPrefix + LAPLACIAN_SUFFIX,
                         options->exportOBJ, options->exportPLY,
                         options->exportOFF, options->exportSTL);
    }
}

Mesh* removeUnwantedPartitions(Mesh* mesh, const AppOptions* options)
{
    // Create an advanced mesh to process the manifold mesh and make it watertight if it has any
    // self intersections
    std::unique_ptr< AdvancedMesh> watertightMesh = std::make_unique< AdvancedMesh >(
                mesh->getVertices(), mesh->getNumberVertices(),
                mesh->getTriangles(), mesh->getNumberTriangles());

    // Deallocate the mesh
    mesh->~Mesh();

    // Split the mesh and return the partition with largest geometry
    AdvancedMesh* partition = watertightMesh->split();

            // Ensures that the mesh is truly two-advanced with no self intersections
    partition->ensureWatertightness();

    // Return a simple mesh with a single partition for optimization
    return partition->toSimpleMesh();
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
                LOG_WARNING("Some partition in the mesh is invalid");
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

        if (!watertightMesh->checkMinDihedralAngle(options->minDihedralAngle))
        {
            watertightMesh->removeTrianglesWithDihedralAngles(options->minDihedralAngle);
        }

        watertightMesh->ensureWatertightness();
    }

    // Verification on the watertightness again
    // watertightMesh->ensureWatertightness();

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

void generateOptimizedMeshWithROI(Mesh *mesh, const AppOptions* options, const ROIs& regions)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        mesh->optimizeAdapttivelyWithROI(options->optimizationIterations,
                                         options->smoothingIterations,
                                         options->flatFactor,
                                         options->denseFactor,
                                         regions);
    }
    else
    {
        mesh->optimizeWithROIs(options->optimizationIterations,
                               options->smoothingIterations,
                               options->denseFactor,
                               regions);
    }

    if (!options->ignoreOptimizedMesh)
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
    if (!options->ignoreWatertightMesh)
        createWatertightMesh(mesh, options);
}

void optimizeMesh(Mesh *mesh, const AppOptions* options)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        mesh->optimizeAdaptively(options->optimizationIterations,
                                 options->smoothingIterations,
                                 options->flatFactor,
                                 options->denseFactor);
    }
    else
    {
        // Default optimization
        if (options->optimizeMeshHomogenous)
            mesh->optimize(options->optimizationIterations,
                           options->smoothingIterations,
                           options->denseFactor);
    }
}

void generateOptimizedMesh(Mesh *mesh, const AppOptions* options)
{
    // Optimize the mesh
    optimizeMesh(mesh, options);

    if (!options->ignoreOptimizedMesh)
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

}

void generateOptimizedMeshWithROI(Mesh *mesh, const AppOptions* options, ROIs regions)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        mesh->optimizeAdapttivelyWithROI(options->optimizationIterations,
                                         options->smoothingIterations,
                                         options->flatFactor,
                                         options->denseFactor, regions);
    }
    else
    {
        // Default optimization
        if (options->optimizeMeshHomogenous)
            mesh->optimizeWithROIs(options->optimizationIterations,
                                   options->smoothingIterations,
                                   options->denseFactor, regions);
    }

    if (!options->ignoreOptimizedMesh)
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

}




Volume* reconstructVolumeFromMesh(Mesh* inputMesh, const AppOptions* options,
                                  const bool& releaseInputMesh)
{
    // Scale the mesh
    inputMesh->scale(options->xScaleFactor,
                     options->yScaleFactor,
                     options->zScaleFactor);

    // Create the volume from the mesh
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
        return DualMarchingCubes::generateMeshFromVolume(volume);
    else
        return MarchingCubes::generateMeshFromVolume(volume);
}

AdvancedMesh* reconstructAdvancedMeshFromVolume(Volume* volume, const AppOptions* options)
{
    if (options->isosurfaceTechnique == DMC_STRING)
        return DualMarchingCubes::generateAdvancedMeshFromVolume(volume);
    else
        return MarchingCubes::generateAdvancedMeshFromVolume(volume);
}

void optimizeMeshWithPartitions(AdvancedMesh* mesh, const AppOptions* options)
{
    // Simple mesh partitions
    std::vector< Mesh* > simplePartitions;

    // Split the mesh into partitions
    std::vector < Ultraliser::AdvancedMesh* > partitions = mesh->splitPartitions(false);

    // Handle the principal partition
    try
    {
        // Convert the advanced mesh to the mesh
        Mesh* simpleMesh = mesh->toSimpleMesh();

        // Release the advanced mesh
        delete mesh;

        // Laplacian smoorhing on a per-partition-basis
        if (options->laplacianIterations > 0)
            simpleMesh->smoothSurface(options->laplacianIterations);

        // Optimize the simple mesh
        optimizeMesh(simpleMesh, options);

        // Create an advanced mesh again for watertightness check
        mesh = new AdvancedMesh(simpleMesh->getVertices(),
                                simpleMesh->getNumberVertices(),
                                simpleMesh->getTriangles(),
                                simpleMesh->getNumberTriangles());

        // Ensure its watertightness
        mesh->ensureWatertightness();
    }
    catch (...)
    {
        LOG_ERROR("There was an issue processing the main partition!");
    }

    // Process the rest of the partitions
    for (size_t i = 0; i < partitions.size(); ++i)
    {
        LOG_SUCCESS("Partition [%d / %d]", i, partitions.size());

        AdvancedMesh* advancedMesh = partitions[i];

        // Convert the advanced mesh to the mesh
        Mesh* simpleMesh = advancedMesh->toSimpleMesh();

        // Release the advanced mesh
        delete advancedMesh;

        // Laplacian smoorhing on a per-partition-basis
        if (options->laplacianIterations > 0)
            simpleMesh->smoothSurface(options->laplacianIterations);

        // Optimize the simple mesh
        try {
            optimizeMesh(simpleMesh, options);
        }  catch (...) {
            LOG_WARNING("Cannot optimize partition [ %d ] of the mesh! Ignoring it", i);
        }

        // Create an advanced mesh again for watertightness check
        advancedMesh = new AdvancedMesh(simpleMesh->getVertices(),
                                        simpleMesh->getNumberVertices(),
                                        simpleMesh->getTriangles(),
                                        simpleMesh->getNumberTriangles());

        // Ensure its watertightness
        advancedMesh ->ensureWatertightness();

        // Put the mesh back to the partitions
        partitions[i] = advancedMesh;
    }

    // Merge back after checking the watertightness
    mesh->appendMeshes(partitions);

    // Free all the partitions
    for (auto partition : partitions)
        delete partition;
}

void generateMarchingCubesMeshArtifacts(const Mesh *mesh, const AppOptions* options)
{   
    // Set the string based on the iso-surface reconstruction algorithm MC or DMC
    std::string artifactString;
    std::string artifactSuffix;
    if (options->isosurfaceTechnique == DMC_STRING)
    {
        artifactString = DMC_STRING; artifactSuffix = DMC_SUFFIX;
    }
    else
    {
        artifactString = MC_STRING; artifactSuffix = MC_SUFFIX;
    }

    // Write the statistics of the reconstructed mesh from the marhcing cubes algorithm
    if (options->writeStatistics)
        mesh->printStats(artifactString, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        mesh->writeDistributions(artifactString, &options->distributionsPrefix);

    // Export the MC mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        mesh->exportMesh(options->meshPrefix + artifactSuffix,
                         options->exportOBJ, options->exportPLY,
                         options->exportOFF, options->exportSTL);
}

Mesh* loadInputMesh(const AppOptions* options)
{
    // Load the mesh and construct the mesh object, and generate its artifacts if needed
    auto mesh = new Mesh(options->inputMeshPath);

    // If the dimensions of the input mesh are larger than

    // Write the statistics of the input mesh
    if (options->writeStatistics)
        mesh->printStats(INPUT_STRING, &options->statisticsPrefix);

    // Write the statistics of the input mesh
    if (options->writeDistributions)
        mesh->writeDistributions(INPUT_STRING, &options->distributionsPrefix);

    return mesh;
}

void generateReconstructedMeshArtifacts(Mesh* mesh, const AppOptions* options)
{
    // Write the mesh reconstructed from the marching cubes algorithm
    if (!options->ignoreMarchingCubesMesh)
        generateMarchingCubesMeshArtifacts(mesh, options);

    // Apply laplacian smoothing
    if (options->laplacianIterations > 0)
        applySmoothingOperator(mesh, options);

    // Create an optimized version of the mesh
    if (options->optimizeMeshHomogenous || options->optimizeMeshAdaptively)
        generateOptimizedMesh(mesh, options);

    // Before creating the watertight mesh, improve the topology even if the mesh was not optimized
    mesh->improveTopology(options->smoothingIterations);

    // Create the final watertight mesh    
    createWatertightMesh(mesh, options);
}

void generateVolumeArtifacts(const Volume* volume, const AppOptions* options)
{
    // Projecting the volume to validate its content
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        volume->project(options->projectionPrefix,
                        options->projectXY, options->projectXZ, options->projectZY,
                        options->projectColorCoded);
    }

    // Write the volume
    if (options->exportBitVolume || options->exportUnsignedVolume || options->exportFloatVolume ||
        options->exportNRRDVolume || options->exportRawVolume )
    {
        volume->writeVolumes(options->volumePrefix,
                             options->exportBitVolume,
                             options->exportUnsignedVolume,
                             options->exportFloatVolume,
                             options->exportNRRDVolume,
                             options->exportRawVolume);
    }

    // Write the stacks
    if (options->exportStackXY || options->exportStackXZ || options->exportStackZY)
    {
        volume->writeStacks(options->outputDirectory + "/" + STACKS_SIRECTORY, options->prefix,
                            options->exportStackXY, options->exportStackXZ, options->exportStackZY);
    }

    // Export volume mesh
    if (options->exportVolumeMesh)
    {
        volume->exportToMesh(options->meshPrefix,
                             options->exportOBJ, options->exportPLY,
                             options->exportOFF, options->exportSTL);
    }

    // Export volume bounding box mesh
    if (options->exportVolumeBoundingBoxMesh)
    {
        volume->exportBoundingBoxMesh(options->meshPrefix,
                                      options->exportOBJ, options->exportPLY,
                                      options->exportOFF, options->exportSTL);
    }

    // Export volume grid mesh
    if (options->exportVolumeGridMesh)
    {
        volume->exportVolumeGridToMesh(options->meshPrefix,
                                       options->exportOBJ, options->exportPLY,
                                       options->exportOFF, options->exportSTL);
    }

    // Print the volume statistics
    if (options->writeStatistics)
    {
        volume->printStats(VOLUME_STRING, &options->statisticsPrefix);
    }
}

}
