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
    std::unique_ptr< AppArguments > args = std::make_unique< AppArguments >(
        argc, argv, COPYRIGHT
        "This application reconstructs a watertight polygonal mesh from an input vasculature "
        "morphology. The generated mesh can be also optimized to reduce the number of "
        "triangles while preserving the volume. "
        "The output mesh is guaranteed in all cases to be two-manifold with no self-intersecting "
        "faces unless the --ignore-self-intersections flag is enabled.");

    args->addInputMorphologyArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addPackingAlgorithmArguments();
    args->addProcessingArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMorphologyArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMeshExportArguments();
    options->verifyMorphologyPrefixArgument();
    options->verifyPackingAlgorithmArgument();
    options->verifyProcessingArguments();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Read the file into a morphology structure
    auto vasculatureMorphology = readVascularMorphology(options->inputMorphologyPath);

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
    size_t resolution;
    if (options->scaledResolution)
        resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                VolumeGrid::getType(options->volumeType));

    // Voxelize morphology
    volume->surfaceVoxelizeVasculatureMorphology(vasculatureMorphology, options->packingAlgorithm);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Generate the reconstructde mesh, scaled and translated to the original location
    auto advancedMesh = reconstructMeshFromVolume(volume, options);

    // Free the volume, we do not need it anymore
    delete volume;

    // Generate the mesh artifacts
    generateReconstructedMeshArtifacts(advancedMesh, options);

    // Compute the ROIs
    // auto regions = vasculatureMorphology->collectRegionsWithThinStructures(0.4);

    // Free
    delete vasculatureMorphology;
    delete advancedMesh;
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
