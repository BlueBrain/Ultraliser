/***************************************************************************************************
 * Copyright (c) 2016 - 2024
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
#include <algorithms/skeletonization/SpineSkeletonizer.h>


// Defines
#define NEURON_SMOOTHING_ITERATIONS 10

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

Mesh* remeshNeuron(Mesh* inputNeuronMesh, AppOptions* options, const bool verbose = VERBOSE)
{
    // Compute the bounding box of the input neuron mesh
    Vector3f pMinInput, pMaxInput;
    inputNeuronMesh->computeBoundingBox(pMinInput, pMaxInput);
    const auto& meshBoundingBox = pMaxInput - pMinInput;

    // Compute the resolution of the volume
    const auto largestDimension = meshBoundingBox.getLargestDimension();
    size_t resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);

    // Construct the volume from the input mesh
    auto volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                             VolumeGrid::getType(options->volumeType), verbose);

    // Apply surface and solid voxelization to the input neuron mesh
    volume->surfaceVoxelization(inputNeuronMesh, false, false, 1.0);
    volume->solidVoxelization(options->voxelizationAxis);

    // Remove the border voxels that span less than half the voxel
    auto bordeVoxels = volume->searchForBorderVoxels();
    for (size_t i = 0; i < bordeVoxels.size(); ++i)
    {
        for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
        {
            auto voxel = bordeVoxels[i][j]; volume->clear(voxel.x(), voxel.y(), voxel.z());
        }
        bordeVoxels[i].clear();
    }
    bordeVoxels.clear();
    volume->surfaceVoxelization(inputNeuronMesh, false, false, 0.5);

    // Construct the mesh using the DMC technique
    auto reconstructedNeuronMesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Smooth the resulting surface mesh
    reconstructedNeuronMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, verbose);

    // Return a pointer to the resulting neuron
    return reconstructedNeuronMesh;
}

Volume* createNeuronVolume(Mesh* neuronMesh, AppOptions* options, const bool verbose = VERBOSE)
{
    // Create the volume from the mesh
    auto neuronVolume = createVolumeGrid(neuronMesh, options);

    // Adaptive and conservative Voxelization
    neuronVolume->surfaceVoxelization(neuronMesh, false, false, 1.0);
    neuronVolume->solidVoxelization(options->voxelizationAxis);

    // Remove the border voxels that span less than half the voxel
    // TODO: VERIFY neuronVolume->surfaceVoxelization(neuronMesh, false, false, 0.5);
    auto bordeVoxels = neuronVolume->searchForBorderVoxels();
    for (size_t i = 0; i < bordeVoxels.size(); ++i)
    {
        for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
        {
            auto voxel = bordeVoxels[i][j];
            neuronVolume->clear(voxel.x(), voxel.y(), voxel.z());
        }
        bordeVoxels[i].clear();
    }
    bordeVoxels.clear();
    neuronVolume->surfaceVoxelization(neuronMesh, false, false, 0.5);

    // Return the volume
    return neuronVolume;
}

Mesh* reconstructNeuronMeshFromVolume(Volume* neuronVolume,
                                      AppOptions* options, const bool verbose = VERBOSE)
{
    // Construct the mesh using the DMC technique
    auto reconstructedNeuronMesh = DualMarchingCubes::generateMeshFromVolume(neuronVolume);

    // Smooth the resulting surface mesh
    reconstructedNeuronMesh->smoothSurface(NEURON_SMOOTHING_ITERATIONS, verbose);

    // Return a pointer to the resulting neuron
    return reconstructedNeuronMesh;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Load the input mesh of the neuron
    auto inputMesh = loadInputMesh(options);

    // Construct the neuron volume
    auto neuronVolume = createNeuronVolume(inputMesh, options, VERBOSE);

    // TODO: Make this as an option
    auto remeshedNeuron = reconstructNeuronMeshFromVolume(neuronVolume, options, VERBOSE);

    // Export optimized neuron mesh
    if (options->exportOptimizedNeuronMesh)
    {
        // Extract the mesh from the volume again
        auto reconstructedMesh = reconstructMeshFromVolume(neuronVolume, options);

        // Generate the artifacts of the mesh
        generateReconstructedMeshArtifacts(reconstructedMesh, options);
    }

    // Project the volume created from the input mesh
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        neuronVolume->project(options->projectionPrefix,
                             options->projectXY, options->projectXZ, options->projectZY);
    }

    // Create a skeletonization object
    NeuronSkeletonizer* skeletonizer = new NeuronSkeletonizer(
                neuronVolume, options->removeSpines, options->useAccelerationStructures,
                options->debugSkeletonization, options->morphologyPrefix);

    // Initialize the skeltonizer
    skeletonizer->initialize();

    // Skeletonize the volume to obtain the centerlines
    skeletonizer->skeletonizeVolumeToCenterLines();

    // Project the center lines of the skeleton
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        const std::string prefix = options->projectionPrefix + SKELETON_SUFFIX;
        neuronVolume->project(prefix, options->projectXY, options->projectXZ, options->projectZY);
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

    {
        auto spineMeshes = skeletonizer->reconstructSpineMeshes(inputMesh, 20, 0.5);

        for (size_t i = 0; i < spineMeshes.size(); ++i)
        {
            auto spineMesh = spineMeshes[i];
            std::stringstream stream;
            stream << options->morphologyPrefix << "_spine_" << i;
            spineMesh->exportMesh(stream.str(), true, false, false, false, SILENT);


            // Construct the spine volume from the mesh
            // Get relaxed bounding box to build the volume
            Vector3f pMinInput, pMaxInput;
            spineMesh->computeBoundingBox(pMinInput, pMaxInput);
            const auto& meshBoundingBox = pMaxInput - pMinInput;

            // Get the largest dimension
            float largestDimension = meshBoundingBox.getLargestDimension();
            auto voxelsPerMicron = 100;
            size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

            // Construct the volume
            Volume* volume = new Volume(pMinInput, pMaxInput, resolution, 0.1,
                                        VOLUME_TYPE::BIT, SILENT);

            // Rasterize the neuron mesh within the bounding box
            volume->surfaceVoxelization(spineMesh);
            volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

            std::stringstream prefixStream;
            prefixStream << options->morphologyPrefix << "_spine" << i;
            volume->projectXY(prefixStream.str());

            std::unique_ptr< SpineSkeletonizer > spineSkeletonizer =
                 std::make_unique< SpineSkeletonizer >(volume, true, false, prefixStream.str());
            spineSkeletonizer->run(SILENT);
        }

        // skeletonizer->_exportSpineExtents(options->morphologyPrefix);
    }

    // Export the somatic proxy mesh
    if (options->exportProxySomaMesh)
    {
        skeletonizer->exportSomaProxyMesh(options->meshPrefix,
            options->exportOBJ, options->exportPLY, options->exportOFF, options->exportSTL);
    }

    // Export the somatic mesh
    if (options->exportSomaMesh)
    {
        skeletonizer->exportSomaMesh(options->meshPrefix,
            options->exportOBJ, options->exportPLY, options->exportOFF, options->exportSTL);
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
