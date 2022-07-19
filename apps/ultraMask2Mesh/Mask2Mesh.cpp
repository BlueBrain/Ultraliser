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
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tool reconstructs a watertight mesh from a .tiff mask "
              "extracted from an EM stack. The .tiff mask is given as a series "
              "of .tif images in a single directory. The reconstructed mesh "
              "can be optimized to create a mesh with nicer topology and less "
              "tessellation.");

    args->addInputMaskDirectoryArguments();
    args->addMaskArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMaskDirectoryArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyMaskDimensionsArguments();
    options->verifyMaskPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Construct a volume from the mask
    Ultraliser::Volume* volume =
            Ultraliser::Volume::constructFromTiffMask(
                options->inputMaskDirectory, options->maskWidth, options->maskHeight,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);


    // Extract the mesh from the volume again
    auto reconstructedMesh = reconstructMeshFromVolume(volume, options);

    // Free the volume, it is not needed any further
    delete volume;

    // If a scale factor is given, not 1.0, scale the mesh
    if ((options->xScaleFactor > 1.f || options->yScaleFactor > 1.f || options->zScaleFactor > 1.f)
     || (options->xScaleFactor < 1.f || options->yScaleFactor < 1.f || options->zScaleFactor < 1.f))
    {
        // Scale the mesh
        reconstructedMesh->scale(options->xScaleFactor, options->yScaleFactor, options->zScaleFactor);
    }

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

