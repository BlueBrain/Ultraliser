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

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< AppArguments > args = std::make_unique < AppArguments >(
        argc, argv, COPYRIGHT
        "This application reconstructs a watertight polygonal mesh from an input "
        "non-watertight mesh. The generated mesh can be also optimized to "
        "reduce the number of triangles while preserving the volume. "
        "The output mesh is guaranteed in all cases to be two-manifold.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addProcessingArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMeshExportArguments();
    options->verifyMeshPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();
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

    // Scale the mesh
    if (options->xScaleFactor > 0.0 || options->yScaleFactor > 0.0 || options->zScaleFactor > 0.0 )
    {
        inputMesh->scale(options->xScaleFactor, options->yScaleFactor, options->zScaleFactor);
    }

    // Creates the volume grid to start processing the mesh
    /// NOTE: The input mesh is release within this operation
    auto volume = reconstructVolumeFromMesh(inputMesh, options);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Extract the mesh from the volume again
    auto reconstructedMesh = reconstructMeshFromVolume(volume, options);

    // Free the volume, it is not needed any further
    delete volume;

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
