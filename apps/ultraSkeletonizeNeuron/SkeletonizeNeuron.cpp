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

#include <Ultraliser.h>
#include <AppCommon.h>
#include <AppArguments.h>
#include <algorithms/skeletonization/NeuronSkeletonizer.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< AppArguments > args = std::make_unique< AppArguments >(
        argc, argv, COPYRIGHT
        "This application reconstructs a high quality neuronal morphology skeleton (directed "
        "acyclic graph) from an input mesh model of the neuron. If this mesh contains spines, the "
        "application is capable of segmenting those spines and reconstruct high quality skeletons "
        "of the individual spines and providing some information on their types and locations. "
        "The application requires a neuronal mesh to create a valid morphology skeleton. "
        "The scale of the input mesh must be microns.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addNeuronalMorphologyExportArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addProcessingArguments();
    args->addSkeletonizationAccelerationArgument();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMeshPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();
    options->verifyNeuronalMorphologyExportArguments();
    options->verifyProcessingArguments();

    // Initialize context, once everything is in place and all the options are verified
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Load the input mesh
    auto inputMesh = loadInputMesh(options);

    // Create the volume from the mesh
    auto solidVolume = createVolumeGrid(inputMesh, options);

    // Adaptive and conservative Voxelization
    solidVolume->surfaceVoxelization(inputMesh, false, false, 1.0);
    solidVolume->solidVoxelization(options->voxelizationAxis);
    solidVolume->surfaceVoxelization(inputMesh, false, false, 0.5);

    // Project the volume created from the input mesh
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        solidVolume->project(options->projectionPrefix,
                             options->projectXY, options->projectXZ, options->projectZY);
    }

    // Create a skeletonization object
    NeuronSkeletonizer* skeletonizer = new NeuronSkeletonizer(
                solidVolume, options->useAccelerationStructures,
                options->debugSkeletonization, options->morphologyPrefix);

    // Initialize the skeltonizer
    skeletonizer->initialize();

    // Skeletonize the volume to obtain the centerlines
    skeletonizer->skeletonizeVolumeToCenterLines();

    // Project the center lines of the skeleton
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        const std::string prefix = options->projectionPrefix + SKELETON_SUFFIX;
        solidVolume->project(prefix, options->projectXY, options->projectXZ, options->projectZY);
    }

    // Construct the neuron graph from the volume
    skeletonizer->constructGraph();

    // Segment the different components of the graph
    skeletonizer->segmentComponents();

    // Export the SWC file of the neuron
    if (options->exportSWC)
    {
        skeletonizer->exportSWCFile(options->morphologyPrefix, options->resampleSkeleton);
    }

    // Export the somatic mesh
    skeletonizer->exportSomaMesh(options->meshPrefix,
                                 options->exportOBJ, options->exportPLY,
                                 options->exportOFF, options->exportSTL);
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
