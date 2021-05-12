/*******************************************************************************
 * Copyright (c) 2016 - 2020
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

Options* parseArguments(Args* args)
{
    Ultraliser::Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The full path to the mesh.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputMesh);

    Ultraliser::Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the volume will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Ultraliser::Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "Just a prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the mesh file "
                "will be used.");
    args->addArgument(&prefix);

    Ultraliser::Argument boundsFile(
                "--bounds-file",
                ARGUMENT_TYPE::STRING,
                "A file that defines the bounding box that will be voxelized.");
    args->addArgument(&boundsFile);

    Ultraliser::Argument voxelizeMesh(
                "--voxelize-input-mesh",
                ARGUMENT_TYPE::BOOL,
                "Voxelize the input mesh.");
    args->addArgument(&voxelizeMesh);

    Ultraliser::Argument volumeResolution(
                "--resolution",
                ARGUMENT_TYPE::INTEGER,
                "The basic resolution of the volume, default 512.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "512");
    args->addArgument(&volumeResolution);

    Ultraliser::Argument autoResolution(
                "--auto-resolution",
                ARGUMENT_TYPE::BOOL,
                "Sets the resolution of the volume based on the mesh dimensions.");
    args->addArgument(&autoResolution);

    Ultraliser::Argument voxelsPerMicron(
                "--voxels-per-micron",
                ARGUMENT_TYPE::INTEGER,
                "Number of voxels per micron in case auto resolution is used, "
                "default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "3");
    args->addArgument(&voxelsPerMicron);

    Ultraliser::Argument edgeGap(
                "--edge-gap",
                ARGUMENT_TYPE::FLOAT,
                "Some little extra space to avoid edges intersection, "
                "default 0.0 for automatic selection.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.0");
    args->addArgument(&edgeGap);

    Ultraliser::Argument projectXY(
                "--project-xy",
                ARGUMENT_TYPE::BOOL,
                "Project an XY projection of the volume.");
    args->addArgument(&projectXY);

    Ultraliser::Argument projectZY(
                "--project-zy",
                ARGUMENT_TYPE::BOOL,
                "Project an ZY projection of the volume.");
    args->addArgument(&projectZY);

    Ultraliser::Argument stackXY(
                "--stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the XY direction.");
    args->addArgument(&stackXY);

    Ultraliser::Argument stackZY(
                "--stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the ZY direction.");
    args->addArgument(&stackZY);

    Ultraliser::Argument writeBitVolume(
                "--write-bit-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a bit volume, where each voxel is stored in "
                "a single bit.");
    args->addArgument(&writeBitVolume);

    Ultraliser::Argument writeNRRDVolume(
                "--write-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Create an NRRD volume that is compatible with VTK.");
    args->addArgument(&writeNRRDVolume);

    Ultraliser::Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in "
                "a single byte.");
    args->addArgument(&writeByteVolume);

    Ultraliser::Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume format to perform the voxelization: "
                "[bit, byte, voxel]. By default, it is a bit volume.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    args->addArgument(&volumeType);

    Ultraliser::Argument solid(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Create a solid volume where the interior is filled.");
    args->addArgument(&solid);

    Ultraliser::Argument optimizationIterations(
                "--optimization-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to optimize the given mesh, "
                "default value 1. If this value is set to 0, the optimization "
                "process will be ignored.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1");
    args->addArgument(&optimizationIterations);

    Ultraliser::Argument smoothingIterations(
                "--smooth-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh, "
                "default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
    args->addArgument(&smoothingIterations);

    Ultraliser::Argument flatFactor(
                "--flat-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseFlat function. "
                "Default value is 0.05.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
    args->addArgument(&flatFactor);

    Ultraliser::Argument denseFactor(
                "--dense-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseDense function. "
                "Default value is 4.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "4.0");
    args->addArgument(&denseFactor);

    Ultraliser::Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OBJ file.");
    args->addArgument(&exportOBJ);

    Ultraliser::Argument exportPLY(
                "--export-ply",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .PLY file.");
    args->addArgument(&exportPLY);

    Ultraliser::Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OFF file.");
    args->addArgument(&exportOFF);

    Ultraliser::Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .STL file.");
    args->addArgument(&exportSTL);

    Ultraliser::Argument writeStatistics(
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
    options->voxelizeMesh = args->getBoolValue(&voxelizeMesh);
    options->solid = args->getBoolValue(&solid);
    options->volumeType = args->getStringValue(&volumeType);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->autoResolution = args->getBoolValue(&autoResolution);
    options->voxelsPerMicron = args->getUnsignedIntegrValue(&voxelsPerMicron);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectZY = args->getBoolValue(&projectZY);
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackZY = args->getBoolValue(&stackZY);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->writeNRRDVolume = args->getBoolValue(&writeNRRDVolume);
    options->optimizationIterations = args->getUnsignedIntegrValue(&optimizationIterations);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->flatFactor = args->getFloatValue(&flatFactor);
    options->denseFactor = args->getFloatValue(&denseFactor);
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportPLY = args->getBoolValue(&exportPLY);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->exportSTL = args->getBoolValue(&exportSTL);
    options->prefix = args->getStringValue(&prefix);
    options->writeStatistics = args->getBoolValue(&writeStatistics);


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

    if (!(options->exportOBJ || options->exportOFF || options->exportSTL))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-off, --export-stl]");
    }

    LOG_TITLE("Ultralizing");
    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_SUCCESS("The bounding box of the volume will be computed "
                    "on the fly!");
        options->boundsFile = EMPTY;
    }
    else
    {
        LOG_SUCCESS("The bounding box of the volume will be loaded from [ %s ]",
                    options->boundsFile.c_str());
    }

    // If no prefix is given, use the file name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Ultraliser::File::getName(options->inputMesh);
    }

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;
    LOG_SUCCESS("Output Directory [ %s ]", options->outputDirectory.c_str());

    // Return the executable options
    return options;
}

int main(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool takes an input mesh that is relatively valid, i.e "
              "with no crappy geometries, but might have holes, non advanced "
              "edges and vertices or have self-intersections. The resulting "
              "mesh will be watertight, and could be optimized too."
              "This tool is better than using the voxelization-based remeshing "
              "in ultraMesh2Mesh in terms of time and space complexity.");

    // Parse the arguments and get the values
    Options* options = parseArguments(&args);

    // Load the mesh
    std::unique_ptr<Ultraliser::AdvancedMesh> inputMesh =
            std::make_unique<Ultraliser::AdvancedMesh>(options->inputMesh);

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
        LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution,
                    largestDimension);

        Ultraliser::Volume* volume = new Ultraliser::Volume(
                    Ultraliser::Vector3f(pMinInput.x, pMinInput.y, pMinInput.z),
                    Ultraliser::Vector3f(pMaxInput.x, pMaxInput.y, pMaxInput.z),
                    resolution, options->edgeGap,
                    Ultraliser::VolumeGrid::getType(options->volumeType));

        // Surface voxelization
        volume->surfaceVoxelization(inputMesh.get());

        // Enable solid voxelization
        if (options->solid)
            volume->solidVoxelization();

        // Projecting the volume to validate its content
        volume->project(options->outputPrefix,
                        options->projectXY, options->projectZY);

        // Write the volume
        volume->writeVolumes(options->outputPrefix,
                             options->writeBitVolume,
                             options->writeByteVolume,
                             options->writeNRRDVolume);

        // Print the volume statistics
        if (options->writeStatistics)
            volume->printStats("Volume", &options->outputPrefix);

        // Write the stacks
        volume->writeStacks(options->outputDirectory, options->prefix,
                            options->stackXY, options->stackZY);

        // Destructor
        volume->~Volume();
    }

    // Write the statistics of the original mesh
    if (options->writeStatistics)
        inputMesh->printStats("input", &options->outputPrefix);

    Ultraliser::Vertex* vertexArray;
    Ultraliser::Triangle* triangleArray;
    uint64_t numberVertices;
    uint64_t numberTriangles;
    inputMesh->getVerticesAndTrianglesArray(vertexArray, triangleArray,
                                            numberVertices, numberTriangles);

    if (options->optimizationIterations > 0)
    {
        // Construct the optimization mesh
        Ultraliser::Mesh* optimizationMesh =
                new Ultraliser::Mesh(numberVertices, numberTriangles);
        for (uint64_t i = 0; i < numberVertices; ++i)
            optimizationMesh->_vertices[i] = vertexArray[i];
        for (uint64_t i = 0; i < numberTriangles; ++i)
            optimizationMesh->_triangles[i] = triangleArray[i];

        optimizationMesh->optimizeAdaptively(options->optimizationIterations,
                                             options->smoothingIterations,
                                             options->flatFactor,
                                             options->denseFactor);

        // Fix the mesh if it has any self intersections
        std::unique_ptr<Ultraliser::AdvancedMesh> advancedMesh =
                std::make_unique<Ultraliser::AdvancedMesh>
                (optimizationMesh->getVertices(),
                 optimizationMesh->getNumberVertices(),
                 optimizationMesh->getTriangles(),
                 optimizationMesh->getNumberTriangles());

        // Free the input mesh
        optimizationMesh->~Mesh();

        // Ensures that the mesh is truly two-advanced with no
        // self intersections
        advancedMesh->ensureWatertightness();

        // Print the mesh statistcs
        if (options->writeStatistics)
            advancedMesh->printStats("watertight",
                                         &options->outputPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->outputPrefix + MANIFOLD_SUFFIX;
        advancedMesh->exportMesh(filePrefix,
                                 options->exportOBJ,
                                 options->exportPLY,
                                 options->exportOFF,
                                 options->exportSTL);
    }
    else
    {
        // Print the mesh statistcs
        if (options->writeStatistics)
            inputMesh->printStats("watertight", &options->outputPrefix);

        // Export the repaired mesh
        std::string filePrefix = options->outputPrefix + MANIFOLD_SUFFIX;
        inputMesh->exportMesh(filePrefix,
                                 options->exportOBJ,
                                 options->exportPLY,
                                 options->exportOFF,
                                 options->exportSTL);
    }

    ULTRALISER_DONE;
}
