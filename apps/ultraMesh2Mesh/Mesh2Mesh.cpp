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

namespace Ultraliser
{

/**
 * @brief createRespectiveDirectories
 * Create the directory tree where the artifacts will be genrated.
 *
 * @param options
 * Tool options
 */
void createRespectiveDirectories(const Options* options)
{
    // Meshes directory
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << MESHES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Volumes directory
    if (options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << VOLUMES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Projections directory
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << PROJECTIONS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Stacks directory
    if (options->stackXY || options->stackXZ || options->stackZY)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << STACKS_SIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Statistics directory
    if (options->writeStatistics)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << STATISTIC_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }
}

/**
 * @brief parseArguments
 * Parse the arguments of the tool.
 *
 * @param args
 * The command line input arguments given by the user.
 *
 * @return
 * User defined options object.
 */
Options* parseArguments(Args* args)
{
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

    Argument writeNRRDVolume(
                "--write-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Create an NRRD volume that is compatible with VTK.");
    args->addArgument(&writeNRRDVolume);

    Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in a single byte.");
    args->addArgument(&writeByteVolume);

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
                "Optimize the reconstructed mesh.");
    args->addArgument(&optimizeMesh);

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
                "Number of iterations to smooth the reconstructed mesh, default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
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
                "4.0");
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

    Argument writeStatistics(
                "--stats",
                ARGUMENT_TYPE::BOOL,
                "Write the statistics.");
    args->addArgument(&writeStatistics);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    // Get all the options
    options->inputMesh = args->getStringValue(&inputMesh);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->boundsFile = args->getStringValue(&boundsFile);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->autoResolution = args->getBoolValue(&autoResolution);
    options->voxelsPerMicron = args->getUnsignedIntegrValue(&voxelsPerMicron);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->writeNRRDVolume = args->getBoolValue(&writeNRRDVolume);
    options->useSolidVoxelization = args->getBoolValue(&useSolidVoxelization);
    options->VoxelizationAxis =
            Volume::getSolidVoxelizationAxis(args->getStringValue(&VoxelizationAxis));
    options->volumeType = args->getStringValue(&volumeType);
    options->optimizeMesh = args->getBoolValue(&optimizeMesh);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->optimizationIterations = args->getUnsignedIntegrValue(&optimizationIterations);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->flatFactor = args->getFloatValue(&flatFactor);
    options->denseFactor = args->getFloatValue(&denseFactor);
    options->ignoreDMCMesh = args->getBoolValue(&ignoreDMCMesh);
    options->ignoreSelfIntersections = args->getBoolValue(&ignoreSelfIntersections);
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportPLY = args->getBoolValue(&exportPLY);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->exportSTL = args->getBoolValue(&exportSTL);
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectXZ = args->getBoolValue(&projectXZ);
    options->projectZY = args->getBoolValue(&projectZY);
    options->projectColorCoded = args->getBoolValue(&projectColorCoded);
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackXZ = args->getBoolValue(&stackXZ);

    options->stackZY = args->getBoolValue(&stackZY);
    options->prefix = args->getStringValue(&prefix);
    options->writeStatistics = args->getBoolValue(&writeStatistics);
    options->useLaplacian = args->getBoolValue(&laplacianFilter);
    options->laplacianIterations = args->getIntegrValue(&laplacianIterations);

    /// Validate the arguments
    if (!Ultraliser::File::exists(options->inputMesh))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMesh.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!",
                  options->outputDirectory.c_str());
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

    // Create the respective directories
    createRespectiveDirectories(options);

    LOG_TITLE("Ultralizing");
    LOG_SUCCESS("Output Directory [ %s ]", options->outputDirectory.c_str());

    // Return the executable options
    return options;
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
    Ultraliser::AdvancedMesh* toBeWatertightMesh = new
        Ultraliser::AdvancedMesh
            (manifoldMesh->getVertices(), manifoldMesh->getNumberVertices(),
             manifoldMesh->getTriangles(), manifoldMesh->getNumberTriangles());

    std::vector < Ultraliser::AdvancedMesh* > partitions = toBeWatertightMesh->splitPartitions();

    // Ensure watertightness for the rest of the partitions
    for (auto mesh : partitions)
        mesh->ensureWatertightness();

    // Ensures that the mesh is truly two-advanced with no self intersections
    toBeWatertightMesh->ensureWatertightness();

    // Merge back after checking the watertightness
    toBeWatertightMesh->appendMeshes(partitions);

    // Free
    for (auto mesh : partitions)
        mesh->~AdvancedMesh();

    // Print the mesh statistcs
    if (options->writeStatistics)
    {
        // Prefix
        const std::string prefix = options->outputDirectory + "/" + STATISTIC_DIRECTORY  + "/";

        // Statistics
        toBeWatertightMesh->printMeshStats(WATERTIGHT_STRING, &prefix);
    }

    // Export the repaired mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        // Prefix
        const std::string prefix = options->outputDirectory + "/" + MESHES_DIRECTORY + "/" +
                options->prefix + WATERTIGHT_SUFFIX;

        // Export
        toBeWatertightMesh->exportMesh(prefix,
                                       options->exportOBJ, options->exportPLY,
                                       options->exportOFF, options->exportSTL);
    }

    // Free
    toBeWatertightMesh->~AdvancedMesh();
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
    // Optimize the mesh adaptively
    dmcMesh->optimizeAdaptively(options->optimizationIterations, options->smoothingIterations,
                                options->flatFactor, options->denseFactor);

    // Print the mesh statistcs
    if (options->writeStatistics)
    {
        // Prefix
        const std::string prefix = options->outputDirectory + "/" + STATISTIC_DIRECTORY +  "/";

        // Statistics
        dmcMesh->printMeshStats(OPTIMIZED_STRING, &prefix);
    }

    // Export the mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        // Prefix
        const std::string prefix = options->outputDirectory + "/" + MESHES_DIRECTORY + "/" +
                options->prefix + OPTIMIZED_SUFFIX;

        // Export
        dmcMesh->exportMesh(prefix,
                            options->exportOBJ, options->exportPLY,
                            options->exportOFF, options->exportSTL);
    }

    // Fix self-intersections if any
    if (!options->ignoreSelfIntersections)
    {
        createMeshWithNoSelfIntersections(dmcMesh, options);
    }
}

/**
 * @brief scaleAndTranslateGeneratedMesh
 * Scale the translate the generated mesh to fit the dimensions of the input mesh.
 *
 * @param generatedMesh
 * The generated mesh from the DMC algorithm.
 * @param inputCenter
 * The center of the input mesh.
 * @param inputBB
 * The bounding box of the input mesh.
 */
void scaleAndTranslateGeneratedMesh(Mesh *generatedMesh,
                                    const Vector3f &inputCenter,
                                    const Vector3f &inputBB)
{
    // Center the reconstructed mesh at the origin
    generatedMesh->centerAtOrigin();

    // Compute the bounding box of the created mesh
    Ultraliser::Vector3f pMaxGenerated, pMinGenerated;
    generatedMesh->computeBoundingBox(pMinGenerated, pMaxGenerated);

    // Compute the scale needed
    const Ultraliser::Vector3f generatedBB = pMaxGenerated - pMinGenerated;
    const Ultraliser::Vector3f scale = inputBB / generatedBB;

    // Scale the mesh
    generatedMesh->scale(scale.x(), scale.y(), scale.z());

    // Translate it back to the original center of the input mesh
    generatedMesh->translate(inputCenter);
}

/**
 * @brief writeDMCMesh
 * Writes the resulting mesh from the DMC algorithm.
 *
 * @param dmcMesh
 * THe generetd mesh from the DMC algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void writeDMCMesh(const Mesh *dmcMesh, const Options* options)
{
    // Write the statistics of the DMC mesh
    if (options->writeStatistics)
    {
        // Prefix
        const std::string prefix =
                options->outputDirectory + "/" + STATISTIC_DIRECTORY +  "/" + options->prefix;

        // Print statistics
        dmcMesh->printMeshStats(DMC_STRING, &prefix);
    }

    // Export the DMC mesh
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        // Prefix
        const std::string prefix = options->outputDirectory + "/" + MESHES_DIRECTORY + "/" +
                options->prefix + DMC_SUFFIX;

        // Export the mesh
        dmcMesh->exportMesh(prefix,
                            options->exportOBJ, options->exportPLY,
                            options->exportOFF, options->exportSTL);
    }
}

/**
 * @brief runMesh2Mesh
 * @param argc
 * @param argv
 */
void runMesh2Mesh(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool reconstructs a watertight polygonal mesh from an input "
              "non-watertight mesh. The generated mesh can be also optimized to "
              "reduce the number of triangles while preserving the volume. "
              "The output mesh is guaranteed in all cases to be two-advanced "
              "with no self-intersecting faces unless the "
              "--ignore-self-intersections flag is enabled.");

    // Parse the arguments and get the tool options
    auto options = parseArguments(&args);

    // Load the mesh and construct the mesh object
    Mesh* inputMesh = new Mesh(options->inputMesh);

    // Write the statistics of the original mesh
    if (options->writeStatistics)
    {
        // Prefix
        const std::string prefix =
                options->outputDirectory + "/" + STATISTIC_DIRECTORY + "/" + options->prefix;

        // Print statistics
        inputMesh->printMeshStats(INPUT_STRING, &prefix);
    }

    // Get relaxed bounding box to build the volume
    Ultraliser::Vector3f pMinInput, pMaxInput;
    inputMesh->computeBoundingBox(pMinInput, pMaxInput);

    // Extend the bounding box a little bit to avoid edge issues
    Ultraliser::Vector3f inputBB = pMaxInput - pMinInput;
    Ultraliser::Vector3f inputCenter = pMinInput + 0.5 * (pMaxInput - pMinInput);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    uint64_t resolution;
    if (options->autoResolution)
    {
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    }
    else
    {
        resolution = options->volumeResolution;
    }
    LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    // Construct the volume
    Volume *volume = new Volume(pMinInput, pMaxInput, resolution, options->edgeGap,
                                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Surface voxelization
    volume->surfaceVoxelization(inputMesh, true, true);

    // Free the input mesh
    inputMesh->~Mesh();

    // Enable solid voxelization
    if (options->useSolidVoxelization)
    {
        volume->solidVoxelization(options->VoxelizationAxis);
    }

    // Projecting the volume to validate its content
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        // Prefix
        const std::string prefix =
                options->outputDirectory + "/" + PROJECTIONS_DIRECTORY +  "/" + options->prefix;

        // Project the volume
        volume->project(prefix,
                        options->projectXY, options->projectXZ, options->projectZY,
                        options->projectColorCoded);
    }

    // Write the volume
    if (options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume)
    {
        // Prefix
        const std::string prefix =
                options->outputDirectory + "/" + VOLUMES_DIRECTORY +  "/" + options->prefix;

        // Write the volume
        volume->writeVolumes(prefix,
                             options->writeBitVolume,
                             options->writeByteVolume,
                             options->writeNRRDVolume);
    }

    // Write the statistics of the reconstructed volume
    if (options->writeStatistics)
    {
        // Prefix
        const std::string prefix =
                options->outputDirectory + "/" + STATISTIC_DIRECTORY +  "/" + options->prefix;

        // Print the volume statistics
        volume->printVolumeStats(VOLUME_STRING, &prefix);
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

    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr<Ultraliser::DualMarchingCubes> dmc =
            std::make_unique<Ultraliser::DualMarchingCubes>(volume);

    // Generate the mesh using the DMC algorithm
    Mesh* generatedMesh = dmc->generateMesh();

    // Free the volume
    volume->~Volume();

    // DMC mesh output
    if (!options->ignoreDMCMesh)
    {
        writeDMCMesh(generatedMesh, options);
    }

    // Laplacian smoorhing
    if (options->useLaplacian)
    {
        // Apply the Laplacian filter
        generatedMesh->applyLaplacianSmooth(options->laplacianIterations,
                                            0.2, 0.1);
        // Prefix
        const std::string prefix  = options->outputDirectory + "/" + MESHES_DIRECTORY + "/" +
                options->prefix + LAPLACIAN_SUFFIX;

        // Export the mesh
        generatedMesh->exportMesh(prefix,
                                  options->exportOBJ,
                                  options->exportPLY,
                                  options->exportOFF,
                                  options->exportSTL);

        // Print the mesh statistcs
        if (options->writeStatistics)
            generatedMesh->printMeshStats(LAPLACIAN_STRING, &options->outputPrefix);
    }

    // Scane and translate the generated mesh to fit the origin mesh
    scaleAndTranslateGeneratedMesh(generatedMesh, inputCenter, inputBB);

    // Optimize the mesh
    if (options->optimizeMesh)
    {
        optimizeMesh(generatedMesh, options);
    }
}

}

int main(int argc , const char** argv)
{
    TIMER_SET;

    Ultraliser::runMesh2Mesh(argc, argv);

    LOG_STATUS_IMPORTANT("Ultralization Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    ULTRALISER_DONE;
}
