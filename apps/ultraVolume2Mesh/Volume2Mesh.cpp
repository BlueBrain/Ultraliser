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
    Ultraliser::Argument volumePath(
                "--volume-path",
                ARGUMENT_TYPE::STRING,
                "The full path to the volume.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&volumePath);

    Ultraliser::Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the meshes or data will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Ultraliser::Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "Just a prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the mesh file "
                "will be used.");
    args->addArgument(&prefix);

    Ultraliser::Argument isoValue(
                "--iso-value",
                ARGUMENT_TYPE::INTEGER,
                "The iso value where the volume will get segmented, default 127",
                ARGUMENT_PRESENCE::OPTIONAL,
                "127");
    args->addArgument(&isoValue);

    Ultraliser::Argument fullRangeIsoValue(
                "--full-range-iso-value",
                ARGUMENT_TYPE::BOOL,
                "If the voxel contains any value, then use it."
                "If this option is set the --iso-value option is ignored.");
    args->addArgument(&fullRangeIsoValue);

    Ultraliser::Argument writeHistogram(
                "--write-histogram",
                ARGUMENT_TYPE::BOOL,
                "Write the histogram of the volume into a text file.");
    args->addArgument(&writeHistogram);

    Ultraliser::Argument zeroPaddingVoxels(
                "--zero-paddgin-voxels",
                ARGUMENT_TYPE::INTEGER,
                "The number of zero-padding voxels that will be appended to "
                "the volume to avoid any clipping artifacts, default 0",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0");
    args->addArgument(&zeroPaddingVoxels);

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
    options->volumePath = args->getStringValue(&volumePath);
    options->isoValue = args->getIntegrValue(&isoValue);
    options->fullRangeIsoValue = args->getBoolValue(&fullRangeIsoValue);
    options->writeHistogram = args->getBoolValue(&writeHistogram);
    options->zeroPaddingVoxels = args->getIntegrValue(&zeroPaddingVoxels);
    options->outputDirectory = args->getStringValue(&outputDirectory);
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

    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!",
                  options->outputDirectory.c_str());
    }

    if (!(options->exportOBJ || options->exportOFF))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-off]");
    }

    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::Directory::getName(options->volumePath);
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
              "This tool reconstructs a watertight mesh from a given volume."
              "The volume is given in .img/.hdr format."
              "The reconstructed mesh can be optimized to create a mesh with "
              "nicer topology and less tessellation.");

    // Parse the arguments and get the values
    Options* options = parseArguments(&args);

    // Construct a volume from the file
    Ultraliser::Volume* volume;
    Ultraliser::Volume* loadedVolume = new Ultraliser::Volume(
                options->volumePath, Ultraliser::VolumeGrid::TYPE::BYTE);

    std::stringstream prefix;
    if (options->fullRangeIsoValue)
        prefix << options->outputPrefix;
    else
        prefix << options->outputPrefix << "_" << options->isoValue;


    if (options->writeHistogram)
    {
        // Create the histogram
        std::vector<uint64_t> histogram =
                Ultraliser::Volume::createHistogram(loadedVolume);

        // Write the histogram to a file
        const std::string path = prefix.str() + std::string(".histogram");
        Ultraliser::File::writeIntegerDistributionToFile(
                    path, histogram);

    }

    if (options->fullRangeIsoValue)
    {
        // Construct a bit volume with a specific iso value
        volume = Ultraliser::Volume::constructFullRangeVolume(
                    loadedVolume,
                    options->zeroPaddingVoxels);

    }
    else
    {
        // Construct a bit volume with a specific iso value
        volume = Ultraliser::Volume::constructIsoValueVolume(
                    loadedVolume, I2UI8(options->isoValue),
                    options->zeroPaddingVoxels);
    }

    // Free the loaded volume
    loadedVolume->~Volume();



    // Project the mask
    volume->project(prefix.str(),
                    options->projectXY, options->projectZY);

    // Write the volume
    volume->writeVolumes(prefix.str(),
                         options->writeBitVolume,
                         options->writeByteVolume);

    // Write the stacks
    volume->writeStacks(options->outputDirectory, prefix.str(),
                        options->stackXY, options->stackZY);

    // Reconstruct a mesh from the mask with DMC
    Ultraliser::DualMarchingCubes* dmc =
            new Ultraliser::DualMarchingCubes(volume);
    Ultraliser::Mesh* mesh = dmc->generateMesh();

    // Free the voulme
    volume->~Volume();

    // If a scale factor is given, not 1.0, scale the mesh
    if (!(Ultraliser::isEqual(options->xScale, 1.f) &&
          Ultraliser::isEqual(options->xScale, 1.f) &&
          Ultraliser::isEqual(options->xScale, 1.f)))
    {
        // Scale the mesh
        mesh->scale(options->xScale, options->yScale, options->zScale);
    }

    // Export the original mesh
    mesh->exportMesh(prefix.str(),
                     options->exportOBJ,
                     options->exportOFF, false);

    // Optimize the mesh if required
    if (options->optimizeMesh)
    {
        // Optimize
        mesh->optimize(options->smoothingIterations,
                       options->smoothingFactor);

        // Export the original mesh
        prefix << "_optimized";
        mesh->exportMesh(prefix.str(),
                         options->exportOBJ,
                         options->exportOFF,
                         false);
    }

    ULTRALISER_DONE;
}

