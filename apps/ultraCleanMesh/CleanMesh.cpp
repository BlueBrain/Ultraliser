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
              "This tool takes an input mesh that is relatively valid, i.e "
              "with no crappy geometries, but might have holes, non advanced "
              "edges and vertices or have self-intersections. The resulting "
              "mesh will be watertight, and could be optimized too."
              "This tool is better than using the voxelization-based remeshing "
              "in ultraMesh2Mesh in terms of time and space complexity.");

    Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The full path to the mesh.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputMesh);

    Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the volume will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "Just a prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the mesh file "
                "will be used.");
    args->addArgument(&prefix);

    Argument boundsFile(
                "--bounds-file",
                ARGUMENT_TYPE::STRING,
                "A file that defines the bounding box that will be voxelized.");
    args->addArgument(&boundsFile);

    Argument voxelizeMesh(
                "--voxelize-input-mesh",
                ARGUMENT_TYPE::BOOL,
                "Voxelize the input mesh.");
    args->addArgument(&voxelizeMesh);

    Argument volumeResolution(
                "--resolution",
                ARGUMENT_TYPE::INTEGER,
                "The basic resolution of the volume, default 512.",
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
                "Number of voxels per micron in case auto resolution is used, "
                "default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "3");
    args->addArgument(&voxelsPerMicron);

    Argument edgeGap(
                "--edge-gap",
                ARGUMENT_TYPE::FLOAT,
                "Some little extra space to avoid edges intersection, "
                "default 0.0 for automatic selection.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.0");
    args->addArgument(&edgeGap);

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
                "Create a bit volume, where each voxel is stored in "
                "a single bit.");
    args->addArgument(&writeBitVolume);

    Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in "
                "a single byte.");
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
                "Specify a volume format to perform the voxelization: "
                "[bit, byte, voxel]. By default, it is a bit volume.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    args->addArgument(&volumeType);

    Argument solid(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Create a solid volume where the interior is filled.");
    args->addArgument(&solid);

    Argument preservePartitions(
                "--preserve-partitions",
                ARGUMENT_TYPE::BOOL,
                "Keeps all the partitions of the mesh in the optimized one.");
    args->addArgument(&preservePartitions);

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
    // Parse the arguments and get the values
    auto options = parseArguments(argc, argv);

    // Load the mesh
    std::unique_ptr<Ultraliser::AdvancedMesh> inputMesh =
            std::make_unique<Ultraliser::AdvancedMesh>(options->inputMesh);

    // Write the statistics of the original mesh
    if (options->writeStatistics)
        inputMesh->printStats(INPUT_STRING, &options->statisticsPrefix);

    // Distributions
    if (options->writeDistributions)
        inputMesh->writeDistributions(INPUT_STRING, &options->distributionsPrefix);

    // Clean mesh and ensure watertightness before optimization to avoid failure
    inputMesh->ensureWatertightness();

    // Voxelize the mesh
    if (options->voxelizeMesh)
    {
        // Get relaxed bounding box to build the volume
        Ultraliser::AdvancedPoint pMinInput, pMaxInput;
        inputMesh->getBoundingBox(pMinInput, pMaxInput);

        // Extend the bounding box a little bit to avoid edge issues
        Ultraliser::AdvancedPoint inputBB = pMaxInput - pMinInput;

        // Get the largest dimension
        float largestDimension = inputBB.x;
        if (inputBB.y > largestDimension)
            largestDimension = inputBB.y;
        if (inputBB.z > largestDimension)
            largestDimension = inputBB.z;

        uint64_t resolution;
        if (options->autoResolution)
            resolution = uint64_t(options->voxelsPerMicron * largestDimension);
        else
            resolution = options->volumeResolution;
        LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

        Ultraliser::Volume* volume = new Ultraliser::Volume(
                    Ultraliser::Vector3f(pMinInput.x, pMinInput.y, pMinInput.z),
                    Ultraliser::Vector3f(pMaxInput.x, pMaxInput.y, pMaxInput.z),
                    resolution, options->edgeGap,
                    Ultraliser::VolumeGrid::getType(options->volumeType));

        // Surface voxelization
        volume->surfaceVoxelization(inputMesh.get());

        // Enable solid voxelization
        if (options->useSolidVoxelization)
            volume->solidVoxelization();

        // Generate the volume artifacts based on the given options
        generateVolumeArtifacts(volume, options);

        // Destructor
        delete volume;
    }

    Ultraliser::Vertex* vertexArray;
    Ultraliser::Triangle* triangleArray;
    uint64_t numberVertices;
    uint64_t numberTriangles;
    inputMesh->getVerticesAndTrianglesArray(vertexArray, triangleArray,
                                            numberVertices, numberTriangles);

    if (options->optimizationIterations > 0)
    {
        // Construct the optimization mesh
        Ultraliser::Mesh* optimizationMesh = new Ultraliser::Mesh(numberVertices, numberTriangles);
        for (uint64_t i = 0; i < numberVertices; ++i)
            optimizationMesh->_vertices[i] = vertexArray[i];
        for (uint64_t i = 0; i < numberTriangles; ++i)
            optimizationMesh->_triangles[i] = triangleArray[i];

        optimizationMesh->optimizeAdaptively(options->optimizationIterations,
                                             options->smoothingIterations,
                                             options->flatFactor,
                                             options->denseFactor);

        // Fix the mesh if it has any self intersections
        std::unique_ptr<Ultraliser::AdvancedMesh> watertightMesh =
                std::make_unique<Ultraliser::AdvancedMesh>
                (optimizationMesh->getVertices(),
                 optimizationMesh->getNumberVertices(),
                 optimizationMesh->getTriangles(),
                 optimizationMesh->getNumberTriangles());

        // Free the input mesh
        delete optimizationMesh;

        // Ensures that the mesh is truly two-advanced with no self intersections
        watertightMesh->ensureWatertightness();

        // Print the mesh statistcs
        if (options->writeStatistics)
            watertightMesh->printStats(WATERTIGHT_STRING, &options->statisticsPrefix);

        // Distributions
        if (options->writeDistributions)
            watertightMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->outputPrefix + MANIFOLD_SUFFIX;
        watertightMesh->exportMesh(filePrefix,
                                 options->exportOBJ,
                                 options->exportPLY,
                                 options->exportOFF,
                                 options->exportSTL);
    }
    else
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            inputMesh->printStats(WATERTIGHT_STRING, &options->statisticsPrefix);

        // Distributions
        if (options->writeDistributions)
            inputMesh->writeDistributions(WATERTIGHT_STRING, &options->distributionsPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->outputPrefix + MANIFOLD_SUFFIX;
        inputMesh->exportMesh(filePrefix,
                                 options->exportOBJ,
                                 options->exportPLY,
                                 options->exportOFF,
                                 options->exportSTL);
    }
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
