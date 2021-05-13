/*******************************************************************************
 * Copyright (c) 2016 - 2019
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 ******************************************************************************/

#include <Ultraliser.h>
#include "Args.h"

/**
 * @brief defineOptions
 * Define the command line arguments options here, parse them and then return a
 * struct with all the options to make it easy to use it in the code.
 * @param parser
 * Command line arguments parser.
 */
Options* parseArguments(Args* args)
{
    Ultraliser::Argument maskDirectory(
                "--mask-directory",
                ARGUMENT_TYPE::STRING,
                "The full path to directory that contains the mask.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&maskDirectory);

    Ultraliser::Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the meshes or data will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Ultraliser::Argument maskWidth(
                "--mask-width",
                ARGUMENT_TYPE::INTEGER,
                "The width of the mask.",
                ARGUMENT_PRESENCE::MANDATORY,
                "0");
    args->addArgument(&maskWidth);

    Ultraliser::Argument maskHeight(
                "--mask-height",
                ARGUMENT_TYPE::INTEGER,
                "The height of the mask.",
                ARGUMENT_PRESENCE::MANDATORY,
                "0");
    args->addArgument(&maskHeight);

    Ultraliser::Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "Just a prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the mesh file "
                "will be used.");
    args->addArgument(&prefix);

    Ultraliser::Argument projectXY(
                "--project-xy",
                ARGUMENT_TYPE::BOOL,
                "Project an XY projection of the volume.");
    args->addArgument(&projectXY);

    Ultraliser::Argument projectZY(
                "--project-zy",
                ARGUMENT_TYPE::BOOL,
                "Project a ZY projection of the volume.");
    args->addArgument(&projectZY);

    Ultraliser::Argument stackXY(
                "--stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Ctrate an image stack along the XY direction.");
    args->addArgument(&stackXY);

    Ultraliser::Argument stackZY(
                "--stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Ctrate an image stack along the ZY direction.");
    args->addArgument(&stackZY);

    Ultraliser::Argument writeBitVolume(
                "--write-bit-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a bit volume, where each voxel is stored in "
                "a single bit.");
    args->addArgument(&writeBitVolume);

    Ultraliser::Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in "
                "a single byte.");
    args->addArgument(&writeByteVolume);

    Ultraliser::Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume to perform the voxelization: "
                "[bit, byte, voxel]. By default, it is a bit volume.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    args->addArgument(&volumeType);

    Ultraliser::Argument solid(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Create a solid volume where the interior is filled.");
    args->addArgument(&solid);

    Ultraliser::Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the meshes to an OBJ file.");
    args->addArgument(&exportOBJ);

    Ultraliser::Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the meshes to an OFF file.");
    args->addArgument(&exportOFF);

    Ultraliser::Argument optimizeMesh(
                "--optimize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Optimize the created mesh.");
    args->addArgument(&optimizeMesh);

    Ultraliser::Argument smoothingIterations(
                "--smooth-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh, "
                "default 10.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "10");
    args->addArgument(&smoothingIterations);

    Ultraliser::Argument smoothingFactor(
                "--smooth-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor used to remove unnecessary geometry from the "
                "reconstructed mesh, by default 10.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "10");
    args->addArgument(&smoothingFactor);

    Ultraliser::Argument xScale(
                "--x-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the X-axis, , default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    args->addArgument(&xScale);

    Ultraliser::Argument yScale(
                "--y-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Y-axis, , default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    args->addArgument(&yScale);

    Ultraliser::Argument zScale(
                "--z-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Z-axis, default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    args->addArgument(&zScale);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    /// Get all the options
    options->maskDirectory = args->getStringValue(&maskDirectory);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->maskWidth = args->getIntegrValue(&maskWidth);
    options->maskHeight = args->getIntegrValue(&maskHeight);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->solid = args->getBoolValue(&solid);
    options->volumeType = args->getStringValue(&volumeType);
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectZY = args->getBoolValue(&projectZY);
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackZY = args->getBoolValue(&stackZY);
    options->optimizeMesh = args->getBoolValue(&optimizeMesh);
    options->smoothingFactor = args->getFloatValue(&smoothingFactor);
    options->smoothingIterations =  args->getIntegrValue(&smoothingIterations);
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->xScale = args->getFloatValue(&xScale);
    options->yScale = args->getFloatValue(&yScale);
    options->zScale = args->getFloatValue(&zScale);
    options->prefix = args->getStringValue(&prefix);


    /// Validate the arguments
    if (!Ultraliser::Directory::exists(options->maskDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist! ",
                  options->maskDirectory.c_str());
    }

    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!",
                  options->outputDirectory.c_str());
    }

    if (options->maskWidth == 0 || options->maskHeight == 0)
    {
        LOG_ERROR("Mask dimensions cannot be zero: [%d x %d]",
                  options->maskWidth, options->maskHeight);
    }

    if (!(options->exportOBJ || options->exportOFF))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-off]");
    }

    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::Directory::getName(options->maskDirectory);
        LOG_INFO("%s", options->prefix.c_str());
    }

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;

    LOG_TITLE("Ultralizing");

    // Return the executable options
    return options;
}

int main(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool reconstructs a watertight mesh from a .tiff mask "
              "extracted from an EM stack. The .tiff mask is given as a series "
              "of .tif images in a single directory. The reconstructed mesh "
              "can be optimized to create a mesh with nicer topology and less "
              "tessellation.");

    // Parse the arguments and get the values
    Options* options = parseArguments(&args);

    // Construct a volume from the mask
    Ultraliser::Volume* maskVolume =
            Ultraliser::Volume::constructFromTiffMask(
                options->maskDirectory, options->maskWidth, options->maskHeight,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Project the mask
    maskVolume->project(options->outputPrefix,
                        options->projectXY, options->projectZY);

    // Print the volume statistics
    maskVolume->printStats(options->prefix, &options->outputPrefix);

    // Write the volume
    maskVolume->writeVolumes(options->outputPrefix,
                             options->writeBitVolume,
                             options->writeByteVolume);

    // Write the stacks
    maskVolume->writeStacks(options->outputDirectory, options->prefix,
                            options->stackXY, options->stackZY);

    // Reconstruct a mesh from the mask with DMC
    Ultraliser::DualMarchingCubes* dmc =
            new Ultraliser::DualMarchingCubes(maskVolume);
    Ultraliser::Mesh* mesh = dmc->generateMesh();

    // Free the voulme
    maskVolume->~Volume();

    // If a scale factor is given, not 1.0, scale the mesh
    if (!(Ultraliser::isEqual(options->xScale, 1.f) &&
          Ultraliser::isEqual(options->xScale, 1.f) &&
          Ultraliser::isEqual(options->xScale, 1.f)))
    {
        // Scale the mesh
        mesh->scale(options->xScale, options->yScale, options->zScale);
    }

    // Print statistics of the reconstructed mesh
    mesh->printStats("dmc", &options->outputPrefix);

    // Export the original mesh
    mesh->exportMesh(options->outputPrefix,
                        options->exportOBJ,
                        options->exportOFF, false);

    // Optimize the mesh if required
    if (options->optimizeMesh)
    {
        // Optimize
        mesh->optimize(options->smoothingIterations,
                       options->smoothingIterations,
                       options->smoothingFactor);

        // Print statistics of the reconstructed mesh
        mesh->printStats("optimized", &options->outputPrefix);

        // Export the original mesh
        mesh->exportMesh(options->outputPrefix,
                         options->exportOBJ,
                         options->exportOFF,
                         false);
    }

    // Free the mesh
    // dmcMesh->~();

    ULTRALISER_DONE;
}

