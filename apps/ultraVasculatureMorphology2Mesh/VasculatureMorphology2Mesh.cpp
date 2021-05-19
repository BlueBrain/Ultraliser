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
#include <AppCommon.h>
#include <AppArguments.h>

namespace Ultraliser
{

Options* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tool reconstructs a watertight polygonal mesh from an input vasculature "
              "morphology. The generated mesh can be also optimized to reduce the number of "
              "triangles while preserving the volume. "
              "The output mesh is guaranteed in all cases to be two-advanced with no self-intersecting "
              "faces unless the --ignore-self-intersections flag is enabled.");

    args->addInputMorphologyArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();

    // Get all the options
    Options* options = args->getOptions();

    LOG_TITLE("Creating Context");

    /// Validate the arguments
    if (!File::exists(options->inputMorphology))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMorphology.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the input will be computed on the fly");
        options->boundsFile = EMPTY;
    }
    else
    {
        LOG_WARNING("The bounding box of the volume will be loaded from [ %s ]",
                    options->boundsFile.c_str());
    }

    // Exporting formats, at least one of them must be there
    if (!(options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-ply, --export-off, --export-stl]");
    }

    // If no prefix is given, use the file name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = File::getName(options->inputMorphology);
    }

    // Initialize context
    initializeContext(options);

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Read the file into a morphology structure
    auto vasculatureMorphology = readVascularMorphology(options->inputMorphology);

    if (options->writeStatistics)
        vasculatureMorphology->printStats(options->prefix, &options->statisticsPrefix);

    if (options->writeDistributions)
        vasculatureMorphology->printDistributions(&options->distributionsPrefix);

    // Get relaxed bounding box to build the volume
    Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
    vasculatureMorphology->getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    // Calculate the volume resolution based on the largest dimension in the morphology
    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                VolumeGrid::getType(options->volumeType));

    // Voxelize morphology
    volume->surfaceVoxelizeVasculatureMorphologyParallel(vasculatureMorphology);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Generate the DMC mesh, scaled and translated to the original location
    auto reconstructedMesh = DualMarchingCubes::generateMeshFromVolume(volume);
    reconstructedMesh->scaleAndTranslate(inputCenter, inputBB);

    // Free the volume, we do not need it anymore
    delete volume;

    // Generate the artifacts of the optimized mesh
    if (!options->ignoreDMCMesh)
        generateDMCMeshArtifacts(reconstructedMesh, options);

    // Laplacian smoorhing
    if (options->useLaplacian)
        applyLaplacianOperator(reconstructedMesh, options);

    // Optimize the mesh
    if (options->optimizeMeshHomogenous || options->optimizeMeshAdaptively)
        generateOptimizedMesh(reconstructedMesh, options);

    // Free
    delete reconstructedMesh;
    delete options;
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
