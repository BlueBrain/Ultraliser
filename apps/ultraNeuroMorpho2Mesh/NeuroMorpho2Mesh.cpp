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

#include <AppArguments.h>
#include <AppCommon.h>
#include <Ultraliser.h>

namespace Ultraliser
{
AppOptions* parseArguments(const int& argc, const char** argv)
{
    // Arguments
    std::unique_ptr<AppArguments> args = std::make_unique<AppArguments>(
        argc, argv,
        "This tool reconstructs a watertight polygonal mesh from an input neuron morphology. "
        "The generated mesh can be also optimized to reduce the number of triangles while "
        "preserving the volume. "
        "The output mesh is guaranteed in all cases to be two-manifold with no self-intersections"
        "unless the --ignore-self-intersections flag is enabled.");

    // Adding the arguments
    args->addInputMorphologyArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addNeuronMorphologyBranchOrderArguments();
    args->addMorphologyAdjustmentParameters();
    args->addMorphologyExportArguments();
    args->addPackingAlgorithmArguments();
    args->addProcessingArguments();
    args->addAutoParametersArguments();


    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMorphologyArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMorphologyPrefixArgument();
    options->verifyPackingAlgorithmArgument();
    options->verifyProcessingArguments();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc, const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Read the file into a morphology structure
    auto neuronMorphology = readNeuronMorphology(options->inputMorphologyPath);

    // Write the statistical anaylsis results of the original morphology
    writeMorphologyAnalysisResults(neuronMorphology, options, "original");

    // Update the minimum sample radius to a given value
    neuronMorphology->verifyMinimumSampleRadius(options->minSampleRadius);

    // Resample the morphology
    neuronMorphology->resampleSectionsAdaptively();

    // Write the statistical anaylsis results of the resampled morphology
    // writeMorphologyAnalysisResults(neuronMorphology, options, "resampled");

    // Trim the neuron morphology following the branch order options
    if (options->axonBranchOrder < INT_MAX |
        options->basalBranchOrder < INT_MAX |
        options->apicalBranchOrder < INT_MAX)
    {
        neuronMorphology->trim(
            options->axonBranchOrder, options->basalBranchOrder, options->apicalBranchOrder);
    }

    // Export the morphology skeleton after the processing
    if (options->exportMorphology)
    {
        neuronMorphology->exportToSWC(options->morphologyPrefix);
    }

    // Get relaxed bounding box to build the volume
    Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
    neuronMorphology->getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    // Calculate the grid resolution to ensure correct voxelization and optimized re-meshing
    size_t resolution;
    float minimumSampleRadius;

    // After the analysis and morphology updates, you can safely now determine the resolution
    if (options->autoParameters)
    {
        minimumSampleRadius = neuronMorphology->computeMinimumSampleRadius();

        // To capture the full details of the neuron, use the minimum sample radius
        const auto voxelSize = minimumSampleRadius * 2.0;

        // Construct the voxels per micron
        const auto voxelsPerMicron = 1.f / voxelSize;

        // Now, you can compute the resolution
        resolution = voxelsPerMicron * largestDimension;
    }
    else
    {
        minimumSampleRadius = options->minSampleRadius;

        // If using the scaled resolution, then detect how many voxels per microns
        if (options->scaledResolution)
            resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);

        // Otherwise, use the volume resolution specified by the user
        else
            resolution = options->volumeResolution;
    }

    // Collect the ROIs where the radii are small
    ROIs regions;
    regions = neuronMorphology->collectRegionsWithThinStructures(minimumSampleRadius * 4.0);

    // Collect the ROIs where the radii are small
    // auto regionsSub = neuronMorphology->collectRegionsWithThinStructures(minimumSampleRadius * 2.0);

    LOG_WARNING("Volume Resolution [%d], Largest Dimension [%f], Min. Sample radius [%f]",
                resolution, largestDimension, minimumSampleRadius);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                VolumeGrid::getType(options->volumeType));

    // Voxelize morphology
    volume->surfaceVoxelizeNeuronMorphology(neuronMorphology, options->packingAlgorithm);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Generate the reconstructde mesh, scaled and translated to the original
    // location
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Free the volume, we do not need it anymore
    delete volume;

    // Refine the ROIs
    // mesh->refineROIs(regionsSub, 1);

    // Write the mesh reconstructed from the marching cubes algorithm
    if (!options->ignoreMarchingCubesMesh)
        generateMarchingCubesMeshArtifacts(mesh, options);

    // Apply laplacian smoothing
    if (options->laplacianIterations > 0)
        applySmoothingOperator(mesh, options);

    // Create an optimized version of the mesh
    if (options->optimizeMeshHomogenous || options->optimizeMeshAdaptively)
        generateOptimizedMeshWithROI(mesh, options, regions);

    // Free
    delete mesh;
    delete options;
}

}  // namespace Ultraliser

int main(int argc, const char** argv)
{
    TIMER_SET;

    Ultraliser::run(argc, argv);

    LOG_STATUS_IMPORTANT("Ultralization Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    ULTRALISER_DONE;
}
