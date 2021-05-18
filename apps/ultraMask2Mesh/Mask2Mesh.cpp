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
    // Arguments
    std::unique_ptr< Args > args = std::make_unique <Args>(argc, argv,
              "This tool reconstructs a watertight mesh from a .tiff mask "
              "extracted from an EM stack. The .tiff mask is given as a series "
              "of .tif images in a single directory. The reconstructed mesh "
              "can be optimized to create a mesh with nicer topology and less "
              "tessellation.");

    Argument maskDirectory(
                "--mask-directory",
                ARGUMENT_TYPE::STRING,
                "The full path to directory that contains the mask.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&maskDirectory);

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

    Argument maskWidth(
                    "--mask-width",
                    ARGUMENT_TYPE::INTEGER,
                    "The width of the mask.",
                    ARGUMENT_PRESENCE::MANDATORY,
                    "0");
        args->addArgument(&maskWidth);

    Argument maskHeight(
                    "--mask-height",
                    ARGUMENT_TYPE::INTEGER,
                    "The height of the mask.",
                    ARGUMENT_PRESENCE::MANDATORY,
                    "0");
        args->addArgument(&maskHeight);

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

    Argument exportVolumeMesh(
                "--export-volume-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export a mesh that represents the volume where each voxel will "
                "be represented by a cube.");
    args->addArgument(&exportVolumeMesh);

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

    Argument optimizeMesh(
                "--optimize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh using the default optimization strategy.");
    args->addArgument(&optimizeMesh);

    Argument adaptiveOptimization(
                "--adaptive-optimization",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh using the adaptive optimization strategy.");
    args->addArgument(&adaptiveOptimization);

    Argument optimizationIterations(
                "--optimization-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to optimize the resulting mesh, default value 1. "
                "If this value is set to 0, the optimization process will be ignored.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1");
    args->addArgument(&optimizationIterations);

    Argument smoothingIterations(
                "--smooth-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh, default 1.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1");
    args->addArgument(&smoothingIterations);

    Argument flatFactor(
                "--flat-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseFlat function, default value is 0.05.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
    args->addArgument(&flatFactor);

    Argument denseFactor(
                "--dense-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseDense function, default value is 4.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5.0");
    args->addArgument(&denseFactor);

    Argument laplacianFilter(
                "--laplacian-filter",
                ARGUMENT_TYPE::BOOL,
                "Use Laplacian filteration to remove the gird artifacts.");
    args->addArgument(&laplacianFilter);

    Argument laplacianIterations(
                "--laplacian-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh with "
                "Laplacian filter, default 3.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "3");
    args->addArgument(&laplacianIterations);

    Argument ignoreSelfIntersections(
                "--ignore-self-intersections",
                ARGUMENT_TYPE::BOOL,
                "Ignore, and take no action if the mesh has self intersecting faces. This process "
                "will speed up the generation of the final mesh, but the output mesh has no "
                "guarntees to be perfectly watertight.");
    args->addArgument(&ignoreSelfIntersections);

    Argument ignoreDMCMesh(
                "--ignore-dmc-mesh",
                ARGUMENT_TYPE::BOOL,
                "Ignore the resulting mesh from the DMC operation.");
    args->addArgument(&ignoreDMCMesh);

    Argument ignoreOptimizedNonWatertightMesh(
                "--ignore-optimized-non-watertight-mesh",
                ARGUMENT_TYPE::BOOL,
                "Ignore the resulting mesh from the optimization process without removing self "
                "intersections.");
    args->addArgument(&ignoreOptimizedNonWatertightMesh);

    Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OBJ file.");
    args->addArgument(&exportOBJ);

    Argument exportPLY(
                "--export-ply",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .PLY file.");
    args->addArgument(&exportPLY);

    Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OFF file.");
    args->addArgument(&exportOFF);

    Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .STL file.");
    args->addArgument(&exportSTL);

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


    Argument preservePartitions(
                "--preserve-partitions",
                ARGUMENT_TYPE::BOOL,
                "Keeps all the partitions of the mesh in the optimized one.");
    args->addArgument(&preservePartitions);

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
    options->inputMaskDirectory = args->getStringValue(&maskDirectory);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->prefix = args->getStringValue(&prefix);

    // Mask dimensions
    options->maskWidth = args->getIntegrValue(&maskWidth);
    options->maskHeight = args->getIntegrValue(&maskHeight);

    // Volume export, file format
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->writeNRRDVolume = args->getBoolValue(&writeNRRDVolume);
    options->exportVolumeMesh = args->getBoolValue(&exportVolumeMesh);

    // Mesh exports, file formats
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportPLY = args->getBoolValue(&exportPLY);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->exportSTL = args->getBoolValue(&exportSTL);

    // Projections
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectXZ = args->getBoolValue(&projectXZ);
    options->projectZY = args->getBoolValue(&projectZY);
    options->projectColorCoded = args->getBoolValue(&projectColorCoded);

    // Stacks
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackXZ = args->getBoolValue(&stackXZ);
    options->stackZY = args->getBoolValue(&stackZY);

    // Mesh scale factors
    options->xScaleFactor = args->getFloatValue(&xScale);
    options->yScaleFactor = args->getFloatValue(&yScale);
    options->zScaleFactor = args->getFloatValue(&zScale);

    // Mesh optimization attributes
    options->optimizeMeshHomogenous = args->getBoolValue(&optimizeMesh);
    options->optimizeMeshAdaptively = args->getBoolValue(&adaptiveOptimization);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->optimizationIterations = args->getUnsignedIntegrValue(&optimizationIterations);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->flatFactor = args->getFloatValue(&flatFactor);
    options->denseFactor = args->getFloatValue(&denseFactor);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->preservePartitions = args->getBoolValue(&preservePartitions);
    options->useLaplacian = args->getBoolValue(&laplacianFilter);
    options->laplacianIterations = args->getIntegrValue(&laplacianIterations);

    // Suppression flags
    options->ignoreDMCMesh = args->getBoolValue(&ignoreDMCMesh);
    options->ignoreSelfIntersections = args->getBoolValue(&ignoreSelfIntersections);
    options->ignoreOptimizedNonWatertightMesh = args->getBoolValue(&ignoreOptimizedNonWatertightMesh);

    // Statistics and distributions
    options->writeStatistics = args->getBoolValue(&writeStatistics);
    options->writeDistributions = args->getBoolValue(&writeDistributions);

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

    if (options->ignoreDMCMesh && options->ignoreSelfIntersections)
    {
        LOG_ERROR("No meshes will be created since you ignored the meshes "
                  "resulting from the DMC stage and also did not use the "
                  "optimization flag to produce an optimized mesh. Enable the "
                  "optimization flag --optimize-mesh to create an optimized "
                  "mesh or remove the --ignore-self-intersections flag.");
    }


    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::Directory::getName(options->inputMaskDirectory);
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

    // DMC mesh output
    if (!options->ignoreDMCMesh)
        generateDMCMeshArtifacts(mesh, options);

    // Laplacian smoorhing
    if (options->useLaplacian)
        applyLaplacianOperator(mesh, options);

    // Optimize the mesh and create a watertight mesh
    if (options->optimizeMeshHomogenous || options->optimizeMeshAdaptively)
        optimizeMesh(mesh, options);

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

