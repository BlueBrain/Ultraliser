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

namespace Ultraliser
{

void applyLaplacianOperator(Mesh *mesh, const Options* options)
{
    // Apply the Laplacian filter
    mesh->applyLaplacianSmooth(options->laplacianIterations, 0.2, 0.1);

    // Export the mesh
    mesh->exportMesh(options->meshPrefix + LAPLACIAN_SUFFIX,
                     options->exportOBJ, options->exportPLY,
                     options->exportOFF, options->exportSTL);

    // Print the mesh statistcs
    if (options->writeStatistics)
        mesh->printStats(LAPLACIAN_STRING, &options->statisticsPrefix);
}

void createWatertightMesh(const Mesh* mesh, const Options* options)
{
    // Create an advanced mesh to process the manifold mesh and make it watertight if it has any
    // self intersections
    std::unique_ptr< Ultraliser::AdvancedMesh> watertightMesh =
        std::make_unique< Ultraliser::AdvancedMesh >(
                mesh->getVertices(), mesh->getNumberVertices(),
                mesh->getTriangles(), mesh->getNumberTriangles());

    if (options->preservePartitions)
    {
        // Split the mesh into partitions
        std::vector < Ultraliser::AdvancedMesh* > partitions = watertightMesh->splitPartitions();

        // Ensure watertightness for the rest of the partitions
        for (auto mesh : partitions)
            mesh->ensureWatertightness();

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

void optimizeMesh(Mesh *dmcMesh, const Options* options)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        dmcMesh->optimizeAdaptively(options->optimizationIterations, options->smoothingIterations,
                                    options->flatFactor, options->denseFactor);

        dmcMesh->smooth();
        dmcMesh->smoothNormals();
    }
    else
    {
        // Default optimization
        if (options->optimizeMeshHomogenous)
            dmcMesh->optimize(options->optimizationIterations,
                              options->smoothingIterations,
                              options->denseFactor);
    }

    if (!options->ignoreOptimizedNonWatertightMesh)
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            dmcMesh->printStats(OPTIMIZED_STRING, &options->statisticsPrefix);

        // Print the mesh statistcs
        if (options->writeDistributions)
            dmcMesh->writeDistributions(OPTIMIZED_STRING, &options->distributionsPrefix);

        // Export the mesh
        if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
            dmcMesh->exportMesh(options->meshPrefix + OPTIMIZED_SUFFIX,
                                options->exportOBJ, options->exportPLY,
                                options->exportOFF, options->exportSTL);

    }

    // Fix self-intersections if any to create the watertight mesh
    if (!options->ignoreSelfIntersections)
        createWatertightMesh(dmcMesh, options);
}

void generateDMCMeshArtifacts(const Mesh *mesh, const Options* options)
{   
    // Write the statistics of the DMC mesh
    if (options->writeStatistics)
        mesh->printStats(DMC_STRING, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        mesh->writeDistributions(DMC_STRING, &options->distributionsPrefix);

    // Export the DMC mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        mesh->exportMesh(options->meshPrefix + DMC_SUFFIX,
                         options->exportOBJ, options->exportPLY,
                         options->exportOFF, options->exportSTL);
}

void generateVolumeArtifacts(const Volume* volume, const Options* options)
{
    // Projecting the volume to validate its content
    if (options->projectXY || options->projectXZ || options->projectZY)
        volume->project(options->projectionPrefix,
                        options->projectXY, options->projectXZ, options->projectZY,
                        options->projectColorCoded);

    // Write the volume
    if (options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume)
        volume->writeVolumes(options->volumePrefix,
                             options->writeBitVolume,
                             options->writeByteVolume,
                             options->writeNRRDVolume);

    // Write the stacks
    if (options->stackXY || options->stackXZ || options->stackZY)
        volume->writeStacks(options->outputDirectory + "/" + STACKS_SIRECTORY, options->prefix,
                            options->stackXY, options->stackXZ, options->stackZY);

    // Export volume mesh
    if (options->exportVolumeMesh)
        volume->exportToMesh(options->meshPrefix,
                             options->exportOBJ, options->exportPLY,
                             options->exportOFF, options->exportSTL);

    // Print the volume statistics
    if (options->writeStatistics)
        volume->printStats(options->prefix, &options->statisticsPrefix);
}

}
