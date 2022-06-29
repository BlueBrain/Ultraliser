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
#include <nrrdloader/NRRDLoader.h>

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

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyIsoSurfaceExtractionArgument();
    options->verifyIsoOptionArgument();


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
    Volume* loadedVolume = new Ultraliser::Volume(options->inputVolumePath);

    loadedVolume->project(options->projectionPrefix + "p",
                    options->projectXY, options->projectXZ, options->projectZY,
                    options->projectColorCoded);

    std::stringstream prefix;
    if (options->fullRangeIsoValue)
        prefix << options->outputPrefix;
    else
        prefix << options->outputPrefix << "-" << options->isoValue;

//    if (options->writeHistogram)
//    {
//        // Create the histogram
//        std::vector<uint64_t> histogram = Volume::createHistogram(loadedVolume,
//                                                                  loadedVolume->getType());

//        // Write the histogram to a file
//        const std::string path = prefix.str() + std::string(".histogram");
//        File::writeIntegerDistributionToFile(path, histogram);
//    }

    // Construct a volume that will be used for the mesh reconstruction
    Ultraliser::Volume* volume;
    if (options->isoOption == ISOVALUE_STRING)
    {
        // Construct a bit volume with a specific iso value
        volume = Volume::constructIsoValueVolume(loadedVolume, options->isoValue);
    }
    else if (options->isoOption == ISOVALUES_STRING)
    {
        // Parse the isp-values file into a list
        const std::vector< uint64_t > isoValues = File::parseIsovaluesFile(options->isovaluesFile);

        volume = Volume::constructIsoValuesVolume(loadedVolume, isoValues);
    }
    else if (options->isoOption == MIN_ISOVALUE_STRING)
    {
        volume = Volume::constructVolumeWithMinimumIsoValue(loadedVolume, options->minIsoValue);
    }
    else if (options->isoOption == MAX_ISOVALUE_STRING)
    {
        volume = Volume::constructVolumeWithMaximumIsoValue(loadedVolume, options->maxIsoValue);
    }
    else if (options->isoOption == FULL_RANGE_STRING)
    {
        // Construct a bit volume with a specific iso value
        volume = Volume::constructFullRangeVolume(loadedVolume);
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
