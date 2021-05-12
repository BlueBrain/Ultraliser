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

#include "Args.h"
#include <Ultraliser.h>
#include <system/System.h>

namespace Ultraliser
{

/**
 * @brief parseArguments
 * Parse the arguments of the tool.
 *
 * @param args
 * The command line input arguments given by the user.
 *
 * @return
 * An object of the user defined options.
 */
Options* parseArguments(Args* args)
{
    Argument inputMorphology(
                "--morphology",
                ARGUMENT_TYPE::STRING,
                "The full path to the vascular morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputMorphology);

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
                "This option is used to select a region of interest from the vascular morphology "
                "to voxelize.");
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

    Argument optimizeMesh(
                "--optimize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh.");
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

    Argument preservePartitions(
                "--preserve-partitions",
                ARGUMENT_TYPE::BOOL,
                "Keeps all the partitions of the mesh in the optimized one.");
    args->addArgument(&preservePartitions);

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

    Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OFF file.");
    args->addArgument(&exportOFF);

    Argument exportPLY(
                "--export-ply",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .PLY file.");
    args->addArgument(&exportPLY);

    Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .STL file.");
    args->addArgument(&exportOFF);

    Argument writeStatistics(
                "--statistics",
                ARGUMENT_TYPE::BOOL,
                "Write data and operation statistics.");
    args->addArgument(&writeStatistics);

    Argument writeDistributions(
                "--distributions",
                ARGUMENT_TYPE::BOOL,
                "Write distributions of attributes of the morphology and the resulting meshes.");
    args->addArgument(&writeDistributions);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    /// Get all the options
    // Input / output
    options->inputMorphology = args->getStringValue(&inputMorphology);
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
    options->VoxelizationAxis =
            Volume::getSolidVoxelizationAxis(args->getStringValue(&VoxelizationAxis));

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

    // Mesh optimization attributes
    options->optimizeMesh = args->getBoolValue(&optimizeMesh);
    options->optimizeMeshAdaptively = args->getBoolValue(&adaptiveOptimization);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->optimizationIterations = args->getUnsignedIntegrValue(&optimizationIterations);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->flatFactor = args->getFloatValue(&flatFactor);
    options->denseFactor = args->getFloatValue(&denseFactor);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->preservePartitions = args->getBoolValue(&preservePartitions);

    // Suppression flags
    options->ignoreDMCMesh = args->getBoolValue(&ignoreDMCMesh);
    options->ignoreSelfIntersections = args->getBoolValue(&ignoreSelfIntersections);
    options->ignoreOptimizedNonWatertightMesh = args->getBoolValue(&ignoreOptimizedNonWatertightMesh);


    // Statistics and distributions
    options->writeStatistics = args->getBoolValue(&writeStatistics);
    options->writeDistributions = args->getBoolValue(&writeDistributions);

    /// Validate the arguments
    if (!File::exists(options->inputMorphology))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMorphology.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    if (!(options->exportOBJ || options->exportOFF))
    {
        LOG_ERROR("The user must specify at least one output format of the mesh to export: "
                  "[--export-obj, --export-off, --export-off, --export-stl]");
    }

    // If no prefix is given, use the file name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = File::getName(options->inputMorphology);
    }

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;

    // Construct the prefixes once and for all
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
    LOG_WARNING("Output Directory [ %s ]", options->outputDirectory.c_str());

    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the input will be computed on the fly");
        options->boundsFile = EMPTY;
    }
    else
    {
        LOG_WARNING("The bounding box of the volume will be loaded from [ %s ]",
                    options->boundsFile.c_str());
    }

    // Return the executable options
    return options;
}

/**
 * @brief generateDMCMeshArtifacts
 * Write
 * @param dmcMesh
 * @param options
 */
void generateDMCMeshArtifacts(const Mesh *dmcMesh, const Options* options)
{
    // Write the statistics of the DMC mesh
    if (options->writeStatistics)
        dmcMesh->printStats(DMC_STRING, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        dmcMesh->writeDistributions(DMC_STRING, &options->distributionsPrefix);

    // Export the DMC mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        dmcMesh->exportMesh(options->meshPrefix + DMC_SUFFIX,
                            options->exportOBJ, options->exportPLY,
                            options->exportOFF, options->exportSTL);
}

/**
 * @brief createMeshWithNoSelfIntersections
 * Process the two manifold mesh and create a watertight mesh using the AdvancedMesh processing.
 *
 * @param manifoldMesh
 * An input mesh that is guaranteed to be two-manifold but have some self intersections.
 * @param options
 * User-defined options given to the executable.
 */
void createMeshWithNoSelfIntersections(const Mesh* manifoldMesh, const Options* options)
{
    // Create an advanced mesh to process the manifold mesh and make it watertight if it has any
    // self intersections
    std::unique_ptr< Ultraliser::AdvancedMesh > toBeWatertightMesh = std::make_unique
        < Ultraliser::AdvancedMesh > (manifoldMesh->getVertices(),
                                      manifoldMesh->getNumberVertices(),
                                      manifoldMesh->getTriangles(),
                                      manifoldMesh->getNumberTriangles());

    if (options->preservePartitions)
    {
        // Split the mesh into partitions
        std::vector < Ultraliser::AdvancedMesh* > partitions =
                toBeWatertightMesh->splitPartitions();

        // Ensure watertightness for the rest of the partitions
        for (auto mesh : partitions)
            mesh->ensureWatertightness();

        // Ensures that the mesh is truly two-manifold with no self intersections
        toBeWatertightMesh->ensureWatertightness();

        // Merge back after checking the watertightness
        toBeWatertightMesh->appendMeshes(partitions);

        // Free
        for (auto mesh : partitions)
            mesh->~AdvancedMesh();
    }
    else
    {
        // Ensures that the mesh is truly two-advanced with no self intersections
        toBeWatertightMesh->ensureWatertightness();
    }

    // Print the mesh statistcs
    if (options->writeStatistics)
        toBeWatertightMesh->printStats(WATERTIGHT_STRING, &options->statisticsPrefix);

    // Print the mesh distributions
    if (options->writeDistributions)
        toBeWatertightMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

    // Export the repaired mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        // Prefix
        const std::string prefix = options->meshPrefix + WATERTIGHT_SUFFIX;

        // Export
        toBeWatertightMesh->exportMesh(prefix,
                                       options->exportOBJ, options->exportPLY,
                                       options->exportOFF, options->exportSTL);
    }
}

/**
 * @brief optimizeMesh
 * Optimize the resulting mesh from the DMC algorithm.
 *
 * @param dmcMesh
 * The resulting mesh from the DMC algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void optimizeMesh(Mesh *dmcMesh, const Options* options)
{
    // Further adaptive optimization
    if (options->optimizeMeshAdaptively)
    {
        dmcMesh->optimizeAdaptively(options->optimizationIterations, options->smoothingIterations,
                                    options->flatFactor, options->denseFactor);
    }
    else
    {
        // Default optimization
        if (options->optimizeMesh)
        {
            dmcMesh->optimize(options->optimizationIterations,
                              options->smoothingIterations,
                              options->denseFactor);
        }
    }

    if (!options->ignoreOptimizedNonWatertightMesh)
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            dmcMesh->printStats(OPTIMIZED_STRING, &options->statisticsPrefix);

        // Print the mesh statistcs
        if (options->writeDistributions)
            dmcMesh->writeDistributions(OPTIMIZED_STRING, &options->distributionsPrefix);

        // Export the mesh
        if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
        {
            // Prefix
            const std::string prefix = options->meshPrefix + OPTIMIZED_SUFFIX;

            // Export
            dmcMesh->exportMesh(prefix,
                                options->exportOBJ, options->exportPLY,
                                options->exportOFF, options->exportSTL);
        }
    }

    // Fix self-intersections if any
    if (!options->ignoreSelfIntersections)
    {
        createMeshWithNoSelfIntersections(dmcMesh, options);
    }
}

void generateVolumeArtifacts(const Volume* volume, const Options* options)
{
    // Projecting the volume to validate its content
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        // Project the volume
        volume->project(options->projectionPrefix,
                        options->projectXY, options->projectXZ, options->projectZY,
                        options->projectColorCoded);
    }

    // Write the volume
    if (options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume)
    {
        // Write the volume
        volume->writeVolumes(options->volumePrefix,
                             options->writeBitVolume,
                             options->writeByteVolume,
                             options->writeNRRDVolume);
    }

    // Write the stacks
    if (options->stackXY || options->stackXZ || options->stackZY)
    {
        // Output directory
        std::string outputDirectory = options->outputDirectory + "/" + STACKS_SIRECTORY;

        // Write the stacks
        volume->writeStacks(outputDirectory, options->prefix,
                            options->stackXY, options->stackXZ, options->stackZY);
    }

    // Export volume mesh
    if (options->exportVolumeMesh)
    {
        // Export the mesh
        volume->exportToMesh(options->meshPrefix,
                             options->exportOBJ, options->exportPLY,
                             options->exportOFF, options->exportSTL);
    }

    // Print the volume statistics
    if (options->writeStatistics)
        volume->printStats(options->prefix, &options->statisticsPrefix);
}

/**
 * @brief run
 * Entry function to run the tool.
 *
 * @param argc
 * Argument count
 *
 * @param argv
 * Arguments array
 *
 * @return
 * EXIT_SUCCESS
 */
int run(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool reconstructs a watertight polygonal mesh from an input vasculature "
              "morphology. The generated mesh can be also optimized to reduce the number of "
              "triangles while preserving the volume. "
              "The output mesh is guaranteed in all cases to be two-advanced with no self-intersecting "
              "faces unless the --ignore-self-intersections flag is enabled.");

    // Parse the arguments and get the values
    auto options = parseArguments(&args);

    // Read the file into a morphology structure
    auto vasculatureMorphology = readVascularMorphology(options->inputMorphology);

    if (options->writeStatistics)
        vasculatureMorphology->printStats(options->prefix, &options->statisticsPrefix);

    if (options->writeDistributions)
        vasculatureMorphology->printDistributions(&options->distributionsPrefix);

    // Get relaxed bounding box to build the volume
    Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
    vasculatureMorphology->getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    // Calculate the volume resolution based on the largest dimension in the morphology
    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                VolumeGrid::getType(options->volumeType));

    // Voxelize morphology
    volume->surfaceVoxelizeVasculatureMorphologyParallel(vasculatureMorphology);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->VoxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Generate the DMC mesh, scaled and translated to the original location
    auto dmcMesh = DualMarchingCubes::generateMeshFromVolume(volume);
    dmcMesh->scaleAndTranslate(inputCenter, inputBB);

    // Free the volume, we do not need it anymore
    volume->~Volume();

    // Generate the artifacts of the optimized mesh
    if (!options->ignoreDMCMesh)
        generateDMCMeshArtifacts(dmcMesh, options);

    // Optimize the mesh
    if (options->optimizeMesh || options->optimizeMeshAdaptively)
    {
        optimizeMesh(dmcMesh, options);
    }

    // Free the DMC mesh
    dmcMesh->~Mesh();

    ULTRALISER_DONE;
}
}

int main(int argc , const char** argv)
{
    return Ultraliser::run(argc, argv);
}
