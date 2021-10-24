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

AppOptions* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tool reconstructs a watertight mesh from a given volume."
              "The volume is given in .img/.hdr format."
              "The reconstructed mesh can be optimized to create a mesh with "
              "nicer topology and less tessellation.");


    args->addInputVolumeArguments();
    args->addInputVolumeParametersArguments();
    args->addVolumeProjectionArguments();
    args->addSolidVoxelizationArguments();
    args->addOutputArguments();
    args->addStacksArguments();
    args->addVolumeExportArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    /// Validate the arguments
    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    // Exporting formats, at least one of them must be there
    if (!(options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj-mesh, --export-ply-mesh, "
                  "                 --export-off-mesh, --export-stl-mesh]");
    }

    // If no prefix is given, use the file name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::File::getName(options->inputVolumePath);
    }

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Construct a volume from the file
    Ultraliser::Volume* loadedVolume = new Ultraliser::Volume(
                options->inputVolumePath, Ultraliser::VolumeGrid::TYPE::BYTE);

    std::stringstream prefix;
    if (options->fullRangeIsoValue)
        prefix << options->outputPrefix;
    else
        prefix << options->outputPrefix << "-" << options->isoValue;

    if (options->writeHistogram)
    {
        // Create the histogram
        std::vector<uint64_t> histogram = Ultraliser::Volume::createHistogram(loadedVolume);

        // Write the histogram to a file
        const std::string path = prefix.str() + std::string(".histogram");
        File::writeIntegerDistributionToFile(path, histogram);
    }

    // Construct a volume that will be used for the mesh reconstruction
    Ultraliser::Volume* volume;
    if (options->fullRangeIsoValue)
    {
        // Construct a bit volume with a specific iso value
        volume = Ultraliser::Volume::constructFullRangeVolume(
                    loadedVolume, options->zeroPaddingVoxels);
    }
    else
    {
        // Construct a bit volume with a specific iso value
        volume = Ultraliser::Volume::constructIsoValueVolume(
                    loadedVolume, I2UI8(options->isoValue), options->zeroPaddingVoxels);
    }

    // Free the loaded volume
    delete loadedVolume;

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Extract the mesh from the volume again
    auto reconstructedMesh = reconstructMeshFromVolume(volume, options);

    // If a scale factor is given, not 1.0, scale the mesh, otherwise avoid the expensive operation
    if (!(isEqual(options->xScaleFactor, 1.f) &&
          isEqual(options->xScaleFactor, 1.f) &&
          isEqual(options->xScaleFactor, 1.f)))
    {
        // Scale the mesh
        reconstructedMesh->scale(options->xScaleFactor,
                                 options->yScaleFactor,
                                 options->zScaleFactor);
    }

    // Free the voulme
    delete volume;

    // Generate the mesh artifacts
    generateMarchingCubesMeshArtifacts(reconstructedMesh, options);

    // Generate the reconstructed mesh artifacts
    generateReconstructedMeshArtifacts(reconstructedMesh, options);

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



