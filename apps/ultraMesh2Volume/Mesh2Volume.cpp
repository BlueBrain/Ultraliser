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

namespace Ultraliser
{

Options* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< Args > args = std::make_unique <Args>(argc, argv,
              "This tool reconstructs a volume from an input mesh.");

    Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The full path to the input mesh that will be remeshed using voxelization.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputMesh);

    Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the results will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "A prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the morphology file will be used.");
    args->addArgument(&prefix);

    Argument boundsFile(
                "--bounds-file",
                ARGUMENT_TYPE::STRING,
                "A file that defines the bounding box that will be voxelized and meshed."
                "This option is used to select a region of interest from the space to voxelize.");
    args->addArgument(&boundsFile);

    Argument volumeResolution(
                "--resolution",
                ARGUMENT_TYPE::INTEGER,
                "The base resolution of the volume, default 512."
                "This resolution is set to the larget dimension of the bounding box of the input "
                "dataset, and the resolution of the other dimensions are computed accordingly.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "512");
    args->addArgument(&volumeResolution);

    Argument autoResolution(
                "--auto-resolution",
                ARGUMENT_TYPE::BOOL,
                "Sets the resolution of the volume based on the mesh dimensions.");
    args->addArgument(&autoResolution);

    Argument voxelsPerMicron(
                "--voxels-per-micron",
                ARGUMENT_TYPE::INTEGER,
                "Number of voxels per micron in case --auto-resolution is used, default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
    args->addArgument(&voxelsPerMicron);

    Argument edgeGap(
                "--edge-gap",
                ARGUMENT_TYPE::FLOAT,
                "Some little extra space to avoid edges intersection, default 0.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.0");
    args->addArgument(&edgeGap);

    Argument projectXY(
                "--project-xy",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the Z-axis.");
    args->addArgument(&projectXY);

    Argument projectXZ(
                "--project-xz",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the Y-axis.");
    args->addArgument(&projectXZ);

    Argument projectZY(
                "--project-zy",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the X-axis.");
    args->addArgument(&projectZY);

    Argument projectColorCoded(
                "--project-color-coded",
                ARGUMENT_TYPE::BOOL,
                "Generate color-coded projections of the volume to help debugging it.");
    args->addArgument(&projectColorCoded);

    Argument stackXY(
                "--stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the Z-axis of the volume.");
    args->addArgument(&stackXY);

    Argument stackXZ(
                "--stack-xz",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the Y-axis of the volume.");
    args->addArgument(&stackXZ);

    Argument stackZY(
                "--stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the X-axis of the volume.");
    args->addArgument(&stackZY);

    Argument writeBitVolume(
                "--write-bit-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a bit volume, where each voxel is stored in a single bit.");
    args->addArgument(&writeBitVolume);

    Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in a single byte.");
    args->addArgument(&writeByteVolume);

    Argument writeNRRDVolume(
                "--write-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Create an NRRD volume that is compatible with VTK.");
    args->addArgument(&writeNRRDVolume);

    Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume format to perform the voxelization: [bit, byte, voxel]. "
                "By default, it is a bit volume to reduce the memory foot print.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    args->addArgument(&volumeType);

    Argument useSolidVoxelization(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Use solid voxelization to fill the interior of the surface volume.");
    args->addArgument(&useSolidVoxelization);

    Argument VoxelizationAxis(
                "--vozelization-axis",
                ARGUMENT_TYPE::STRING,
                "The axis where solid voxelization operation will be performed. "
                "Use one of the following options [x, y, z, or xyz]. "
                "If you use x or y or z the voxelization will happen on a single axis, "
                "otherwise, using xyz will perform the solid voxelization along the three main "
                "axes of the volume to avoid filling any loops in the morphology."
                "By default, the Z-axis solid voxelization with xyz is applied if the --solid "
                "flag is set.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "z");
    args->addArgument(&VoxelizationAxis);

    Argument writeStatistics(
                "--stats",
                ARGUMENT_TYPE::BOOL,
                "Write the statistics.");
    args->addArgument(&writeStatistics);

    Argument writeDistributions(
                "--dists",
                ARGUMENT_TYPE::BOOL,
                "Write the distributions.");
    args->addArgument(&writeDistributions);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    /// Get all the options
    // Input / output
    options->inputMesh = args->getStringValue(&inputMesh);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->prefix = args->getStringValue(&prefix);

    // Bounds file for ROI
    options->boundsFile = args->getStringValue(&boundsFile);

    // Volume attributes
    options->autoResolution = args->getBoolValue(&autoResolution);
    options->voxelsPerMicron = args->getUnsignedIntegrValue(&voxelsPerMicron);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->volumeType = args->getStringValue(&volumeType);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->useSolidVoxelization = args->getBoolValue(&useSolidVoxelization);
    options->VoxelizationAxis = Volume::getSolidVoxelizationAxis(args->getStringValue(&VoxelizationAxis));

    // Volume export, file format
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->writeNRRDVolume = args->getBoolValue(&writeNRRDVolume);

    // Statistics and distributions
    options->writeStatistics = args->getBoolValue(&writeStatistics);
    options->writeDistributions = args->getBoolValue(&writeDistributions);


    LOG_TITLE("Creating Context");

    /// Validate the arguments
    if (!Ultraliser::File::exists(options->inputMesh))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMesh.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    /// Validate the arguments
    if (!Ultraliser::File::exists(options->inputMesh))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMesh.c_str());
    }

    if (!(options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "volume: [--write-bit-volume, --write-byte-volume, --write-nrrd-volume]");
    }

    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the volume will be computed on the fly");
        options->boundsFile = EMPTY;
    }
    else
    {
        LOG_WARNING("The bounding box of the volume will be loaded from [ %s ]",
                    options->boundsFile.c_str());
    }

    // If no prefix is given, use the file name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::File::getName(options->inputMesh);
    }

    // Construct the prefixes once and for all
    options->outputPrefix =
            options->outputDirectory + "/" + options->prefix;
    options->meshPrefix =
            options->outputDirectory + "/" + MESHES_DIRECTORY +  "/" + options->prefix;
    options->volumePrefix =
            options->outputDirectory + "/" + VOLUMES_DIRECTORY +  "/" + options->prefix;
    options->projectionPrefix =
            options->outputDirectory + "/" + PROJECTIONS_DIRECTORY +  "/" + options->prefix;
    options->statisticsPrefix =
            options->outputDirectory + "/" + STATISTICS_DIRECTORY +  "/" + options->prefix;
    options->distributionsPrefix =
            options->outputDirectory + "/" + DISTRIBUTIONS_DIRECTORY +  "/" + options->prefix;

    // Create the respective directories
    createRespectiveDirectories(options);

    LOG_TITLE("Ultralizing");
    LOG_SUCCESS("Output Directory [ %s ]", options->outputDirectory.c_str());

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Load the mesh
    Mesh* inputMesh = new Mesh(options->inputMesh);

    // Write the statistics of the original mesh
    if (options->writeStatistics)
        inputMesh->printStats(INPUT_STRING, &options->statisticsPrefix);

    // Write the statistics of the original mesh
    if (options->writeDistributions)
        inputMesh->writeDistributions(INPUT_STRING, &options->statisticsPrefix);

    // Get relaxed bounding box to build the volume
    Ultraliser::Vector3f pMinInput, pMaxInput;
    inputMesh->computeBoundingBox(pMinInput, pMaxInput);

    // Extend the bounding box a little bit to avoid edge issues
    Ultraliser::Vector3f inputBB = pMaxInput - pMinInput;

    // Get the largest dimension
    const float largestDimension = inputBB.getLargestDimension();

    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_SUCCESS("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    Volume *volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Surface voxelization
    volume->surfaceVoxelization(inputMesh, true, true);

    // Free the input mesh
    delete inputMesh;

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->VoxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);
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

