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
#include <AppCommon.h>
#include <AppArguments.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique < AppArguments >(
        argc, argv, COPYRIGHT
        "This application takes an input mesh that is relatively valid, i.e with no crappy "
        "geometries, but might have holes, non-manifold edges and vertices or have "
        "self-intersections. The resulting mesh will be watertight, and could be optimized too."
        "The application is better than using the voxelization-based remeshing "
        "in ultraMesh2Mesh in terms of time and space complexity, but it is not guranteed to work "
        "if the mesh is severely fragmented, for example if a triangle soupe is given.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addMeshVoxelizationArgument();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyBoudsFileArgument();
    options->verifyMeshPrefixArgument();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc, argv);

    // Load the mesh
    std::unique_ptr<Ultraliser::AdvancedMesh> inputMesh =
            std::make_unique<Ultraliser::AdvancedMesh>(options->inputMeshPath);

    // Write the statistics of the original mesh
    if (options->writeStatistics)
        inputMesh->printStats(INPUT_STRING, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        inputMesh->writeDistributions(INPUT_STRING, &options->distributionsPrefix);

    if (options->preservePartitions)
    {
        // Split the mesh into partitions
        std::vector < Ultraliser::AdvancedMesh* > partitions = inputMesh->splitPartitions();

        // Ensure watertightness for the rest of the partitions
        for (auto mesh : partitions)
            mesh->ensureWatertightness();

        // Ensures that the mesh is truly two-manifold with no self intersections
        inputMesh->ensureWatertightness();

        // Merge back after checking the watertightness
        inputMesh->appendMeshes(partitions);

        // Free
        for (auto mesh : partitions)
            delete mesh;
    }
    else
    {
        // Ensures that the mesh is truly two-advanced with no self intersections
        inputMesh->ensureWatertightness();
    }

    // Voxelize the mesh
    if (options->voxelizeMesh)
    {
        // Get relaxed bounding box to build the volume
        Ultraliser::AdvancedPoint pMinInput, pMaxInput;
        inputMesh->getBoundingBox(pMinInput, pMaxInput);

        // Extend the bounding box a little bit to avoid edge issues
        Ultraliser::AdvancedPoint inputBB = pMaxInput - pMinInput;

        // Get the largest dimension
        float largestDimension = inputBB.x;
        if (inputBB.y > largestDimension)
            largestDimension = inputBB.y;
        if (inputBB.z > largestDimension)
            largestDimension = inputBB.z;

        size_t resolution;
        if (options->scaledResolution)
            resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);
        else
            resolution = options->volumeResolution;
        LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

        Ultraliser::Volume* volume = new Ultraliser::Volume(
                    Ultraliser::Vector3f(pMinInput.x, pMinInput.y, pMinInput.z),
                    Ultraliser::Vector3f(pMaxInput.x, pMaxInput.y, pMaxInput.z),
                    resolution, options->edgeGap,
                    Ultraliser::VolumeGrid::getType(options->volumeType));

        // Surface voxelization
        volume->surfaceVoxelization(inputMesh.get());

        // Enable solid voxelization
        if (options->useSolidVoxelization)
            volume->solidVoxelization();

        // Generate the volume artifacts based on the given options
        generateVolumeArtifacts(volume, options);

        // Destructor
        delete volume;
    }

    Ultraliser::Vertex* vertexArray;
    Ultraliser::Triangle* triangleArray;
    size_t numberVertices;
    size_t numberTriangles;
    inputMesh->getVerticesAndTrianglesArray(vertexArray, triangleArray,
                                            numberVertices, numberTriangles);

    if (options->optimizationIterations > 0)
    {
        // Construct the optimization mesh
        Ultraliser::Mesh* optimizationMesh = new Ultraliser::Mesh(numberVertices, numberTriangles);
        for (size_t i = 0; i < numberVertices; ++i)
            optimizationMesh->_vertices[i] = vertexArray[i];
        for (size_t i = 0; i < numberTriangles; ++i)
            optimizationMesh->_triangles[i] = triangleArray[i];

        optimizationMesh->optimizeAdaptively(options->optimizationIterations,
                                             options->smoothingIterations,
                                             options->flatFactor,
                                             options->denseFactor);

        // Fix the mesh if it has any self intersections
        std::unique_ptr<Ultraliser::AdvancedMesh> watertightMesh =
                std::make_unique<Ultraliser::AdvancedMesh>
                (optimizationMesh->getVertices(),
                 optimizationMesh->getNumberVertices(),
                 optimizationMesh->getTriangles(),
                 optimizationMesh->getNumberTriangles());

        // Free the input mesh
        delete optimizationMesh;

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

        // Distributions
        if (options->writeDistributions)
            watertightMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->meshPrefix + MANIFOLD_SUFFIX;
        watertightMesh->exportMesh(filePrefix,
                                   options->exportOBJ,
                                   options->exportPLY,
                                   options->exportOFF,
                                   options->exportSTL);
    }
    else
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            inputMesh->printStats(WATERTIGHT_STRING, &options->statisticsPrefix);

        // Distributions
        if (options->writeDistributions)
            inputMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->outputPrefix + MANIFOLD_SUFFIX;
        inputMesh->exportMesh(filePrefix,
                                 options->exportOBJ,
                                 options->exportPLY,
                                 options->exportOFF,
                                 options->exportSTL);
    }
}

}

int main(int argc , const char** argv)
{
    TIMER_SET;

    Ultraliser::run(argc, argv);

    LOG_STATUS_IMPORTANT("Ultralization Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    ULTRALISER_DONE;
}
