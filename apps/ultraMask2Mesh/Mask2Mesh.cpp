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
    Options* options = args->getOptions();

    LOG_TITLE("Creating Context");

    /// Validate the arguments
    if (!Ultraliser::Directory::exists(options->inputMaskDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist! ", options->inputMaskDirectory.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    if (options->maskWidth == 0 || options->maskHeight == 0)
    {
        LOG_ERROR("Mask dimensions cannot be zero: [%d x %d]",
                  options->maskWidth, options->maskHeight);
    }

    // Exporting formats, at least one of them must be there
    if (!(options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-ply, --export-off, --export-stl]");
    }

    if (!options->writeMarchingCubeMesh && options->ignoreLaplacianSmoothing
            && options->ignoreSelfIntersections)
    {
        LOG_ERROR("No meshes will be created since you ignored the meshes "
                  "resulting from the marching cubes stage and also did not use the "
                  "optimization flag to produce an optimized mesh. Enable the "
                  "optimization flag --optimize-mesh to create an optimized "
                  "mesh or remove the --ignore-self-intersections flag.");
    }

    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::Directory::getName(options->inputMaskDirectory);
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

    // Construct a volume from the mask
    Ultraliser::Volume* volume =
            Ultraliser::Volume::constructFromTiffMask(
                options->inputMaskDirectory, options->maskWidth, options->maskHeight,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);


    // Generate the mesh using the DMC algorithm and adjust its scale
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Free the volume, it is not needed any further
    delete volume;

    // If a scale factor is given, not 1.0, scale the mesh
    if (!(Ultraliser::isEqual(options->xScaleFactor, 1.f) &&
          Ultraliser::isEqual(options->xScaleFactor, 1.f) &&
          Ultraliser::isEqual(options->xScaleFactor, 1.f)))
    {
        // Scale the mesh
        mesh->scale(options->xScaleFactor, options->yScaleFactor, options->zScaleFactor);
    }

    // Generate the mesh artifacts
    generateMeshArtifacts(mesh, options);

    // Free
    delete mesh;
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

