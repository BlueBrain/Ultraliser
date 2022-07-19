/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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
              "ultraVolume2Mesh reconstructs a watertight mesh from a given volume."
              "The volume can be a .NRRD file, a RAW file in .HDR/.IMG or an Ulrraliser-specific "
              "volume in .UVOL or .UVOLB. The reconstructed mesh is optimized to create clean "
              "topology and less tessellation.");

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

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyIsoSurfaceExtractionArgument();
    options->verifyIsoOptionArgument();
    options->verifyVolumePrefixArgument();

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
    Volume* loadedVolume = new Ultraliser::Volume(options->inputVolumePath);

    // Compute the projection of the loaded volume (used for verification)
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        loadedVolume->project(options->projectionPrefix + "-in",
                              options->projectXY, options->projectXZ, options->projectZY,
                              options->projectColorCoded);
    }

    // Compute the projection of the histogram
    if (options->writeHistogram)
    {
        // Create the histogram
        std::vector< size_t > histogram = Volume::createHistogram(loadedVolume,
                                                                  loadedVolume->getType());

        // Write the histogram to a file
        const std::string path = options->outputPrefix + HISTOGRAM_EXTENSION;
        File::writeIntegerDistributionToFile(path, histogram);
    }

    // Construct the iso-volume that will be used for the mesh reconstruction
    Ultraliser::Volume* volume;
    if (options->isoOption == ISOVALUE_STRING)
    {
        // Construct a volume with a specific iso value
        volume = Volume::constructIsoValueVolume(loadedVolume, options->isoValue);
    }
    else if (options->isoOption == ISOVALUES_STRING)
    {
        // Parse the iso-values from the file into a list
        const std::vector< size_t > isoValues = File::parseIsovaluesFile(options->isovaluesFile);

        // Construct a volume with a list of values
        volume = Volume::constructIsoValuesVolume(loadedVolume, isoValues);
    }
    else if (options->isoOption == MIN_ISOVALUE_STRING)
    {
        // Construct a volume with a minimum value
        volume = Volume::constructVolumeWithMinimumIsoValue(loadedVolume, options->minIsoValue);
    }
    else if (options->isoOption == MAX_ISOVALUE_STRING)
    {
        // Construct a volume with a maximum value
        volume = Volume::constructVolumeWithMaximumIsoValue(loadedVolume, options->maxIsoValue);
    }
    else if (options->isoOption == ISOVALUE_RANGE_STRING)
    {
        // Construct a volume with a given range
        volume = Volume::constructVolumeWithIsoRange(loadedVolume,
                                                     options->minIsoValue, options->maxIsoValue);
    }
    else if (options->isoOption == NON_ZERO_STRING)
    {
        // Construct a volume containing all the non-zero voxels of the loaded volume
        volume = Volume::constructNonZeroVolume(loadedVolume);
    }
    else
    {
        LOG_ERROR("The selected isooption [%s] is NOT valid.", options->isoOption.c_str());
    }

    Vector3f scale;
    scale.x() = loadedVolume->getScale().x(); //loadedVolume->getWidth();
    scale.y() = loadedVolume->getScale().y(); // loadedVolume->getHeight();
    scale.z() = loadedVolume->getScale().z(); // loadedVolume->getDepth();

    Vector3f center = loadedVolume->getCenter();

    // Free the loaded volume
    delete loadedVolume;

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Extract the mesh from the volume again
    auto reconstructedMesh = reconstructMeshFromVolume(volume, options);

    /*
    reconstructedMesh->scale(scale.x(), scale.y(), scale.z());
    reconstructedMesh->translate(center);
*/
    // If a scale factor is given, not 1.0, scale the mesh, otherwise avoid the expensive operation
    if (!(isEqual(options->xScaleFactor, 1.f) &&
          isEqual(options->xScaleFactor, 1.f) &&
          isEqual(options->xScaleFactor, 1.f)))
    {
        reconstructedMesh->scale(options->xScaleFactor, options->yScaleFactor, options->zScaleFactor);
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
