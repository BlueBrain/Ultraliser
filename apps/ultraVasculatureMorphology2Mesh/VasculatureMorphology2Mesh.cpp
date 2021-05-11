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
#include "Args.h"

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

    Argument projectColorCodedProjections(
                "--color-coded-projection",
                ARGUMENT_TYPE::BOOL,
                "Generate color-coded projections of the volume to help debugging it.");
    args->addArgument(&projectColorCodedProjections);

    Argument stackXY(
                "--stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the Z-axis direction.");
    args->addArgument(&stackXY);

    Argument stackXZ(
                "--stack-xz",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the Y-axis direction.");
    args->addArgument(&stackXZ);

    Argument stackZY(
                "--stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Create an image stack along the X-axis direction.");
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

    Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume format to perform the voxelization: "
                "[bit, byte, voxel]. By default, it is a bit volume.",
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

    Argument smoothingIterations(
                "--smooth-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh, default 10.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "10");
    args->addArgument(&smoothingIterations);

    Argument smoothingFactor(
                "--smooth-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor used to remove unnecessary geometry from the "
                "reconstructed mesh, by default 10.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "10");
    args->addArgument(&smoothingFactor);

    Argument ignoreSelfIntersections(
                "--ignore-self-intersections",
                ARGUMENT_TYPE::BOOL,
                "Ignore, and take no action if the mesh has self intersecting faces. "
                "This process will speed up the generation of the final mesh, but the output mesh "
                "has no guarntees to be perfectly watertight.");
    args->addArgument(&ignoreSelfIntersections);

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

    Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .STL file.");
    args->addArgument(&exportOFF);

    Argument writeStatistics(
                "--stats",
                ARGUMENT_TYPE::BOOL,
                "Write data and operation statistics.");
    args->addArgument(&writeStatistics);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    // Get all the options
    options->inputMorphology = args->getStringValue(&inputMorphology);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->boundsFile = args->getStringValue(&boundsFile);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->useSolidVoxelization = args->getBoolValue(&useSolidVoxelization);
    options->VoxelizationAxis =
            Volume::getSolidVoxelizationAxis(args->getStringValue(&VoxelizationAxis));
    options->volumeType = args->getStringValue(&volumeType);
    options->optimizeMesh = args->getBoolValue(&optimizeMesh);
    options->smoothingFactor = args->getFloatValue(&smoothingFactor);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->ignoreSelfIntersections = args->getBoolValue(&ignoreSelfIntersections);
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->exportSTL = args->getBoolValue(&exportOFF);
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectXZ = args->getBoolValue(&projectXZ);
    options->projectZY = args->getBoolValue(&projectZY);
    options->projectColorCodedProjections = args->getBoolValue(&projectColorCodedProjections);
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackXZ = args->getBoolValue(&stackXZ);
    options->stackZY = args->getBoolValue(&stackZY);
    options->prefix = args->getStringValue(&prefix);
    options->writeStatistics = args->getBoolValue(&writeStatistics);

    /// Validate the arguments
    if (!File::exists(options->inputMorphology))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMorphology.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!",
                  options->outputDirectory.c_str());
    }

    if (!(options->exportOBJ || options->exportOFF))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-off]");
    }

    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the volume will be computed "
                    "on the fly");
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
        options->prefix = File::getName(options->inputMorphology);
    }

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;

    LOG_TITLE("Ultralizing");
    LOG_STATUS("Output Directory [ %s ]", options->outputDirectory.c_str());

    // Return the executable options
    return options;
}

VasculatureMorphology* readVascularMorphology(std::string morphologyFilePath)
{

    if (String::subStringFound(morphologyFilePath, std::string(".h5")) ||
            String::subStringFound(morphologyFilePath, std::string(".H5")))
    {
        // Read the file
        auto reader = std::make_unique<VasculatureH5Reader>(morphologyFilePath);

        // Get a pointer to the morphology to start using it
        return reader->getMorphology();
    }
    else if (String::subStringFound(morphologyFilePath, std::string(".vmv")) ||
             String::subStringFound(morphologyFilePath, std::string(".VMV")))
    {
        // Read the file
        auto reader = std::make_unique<VasculatureVMVReader>(morphologyFilePath);

        // Get a pointer to the morphology to start using it
        return reader->getMorphology();
    }
    else
        LOG_ERROR("Unrecognized file format");
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
              "The output mesh is guaranteed in all cases to be two-advanced "
              "with no self-intersecting faces unless the "
              "--ignore-self-intersections flag is enabled.");

    // Parse the arguments and get the values
    auto options = parseArguments(&args);

    // Read the file
    auto vasculatureMorphology = readVascularMorphology(options->inputMorphology);

    // Get relaxed bounding box to build the volume
    Vector3f pMinInput, pMaxInput, boundsInput, centerInput;
    vasculatureMorphology->getBoundingBox(pMinInput, pMaxInput,
                                          boundsInput, centerInput);

    if (options->writeStatistics)
        vasculatureMorphology->printMorphologyStats("Morphology", &options->outputPrefix);

    // Extend the bounding box a little bit to avoid edge issues
    Vector3f inputBB = pMaxInput - pMinInput;
    Vector3f inputCenter = pMinInput + 0.5 * ( pMaxInput - pMinInput);

    uint64_t resolution = uint64_t(inputBB.x());
    if (inputBB.y() > resolution)
        resolution = uint64_t(inputBB.y());
    if (inputBB.z() > resolution)
        resolution = uint64_t(inputBB.z());

    std::unique_ptr< Volume > volume =
            std::make_unique<Volume>(pMinInput, pMaxInput,
                                     options->volumeResolution, options->edgeGap,
                                     VolumeGrid::getType(options->volumeType));

    // Voxelize morphology
    volume->surfaceVoxelizeVasculatureMorphology(vasculatureMorphology, false, false);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->VoxelizationAxis);

    // Projecting the volume to validate its content
    volume->project(options->outputPrefix,
                    options->projectXY, options->projectXZ, options->projectZY,
                    options->projectColorCodedProjections);

    // Write the volume
    volume->writeVolumes(options->outputPrefix,
                         options->writeBitVolume,
                         options->writeByteVolume);

    // Print the volume statistics
    if (options->writeStatistics)
        volume->printVolumeStats("Volume", &options->outputPrefix);

    // Write the stacks
    volume->writeStacks(options->outputDirectory, options->prefix,
                        options->stackXY, options->stackXZ, options->stackZY);

    // Reconstruct a watertight mesh from the volume with DMC
    auto dmc = new DualMarchingCubes(volume.get());

    // Generate the mesh
    auto generatedMesh = dmc->generateMesh();

    // Center the reconstructed mesh at the origin
    generatedMesh->centerAtOrigin();

    // Compute the bounding box of the created mesh
    Vector3f pMaxGenerated, pMinGenerated;
    generatedMesh->computeBoundingBox(pMinGenerated, pMaxGenerated);

    const auto generatedBB = pMaxGenerated - pMinGenerated;
    const auto scale = inputBB / generatedBB;

    // Scale the mesh
    generatedMesh->scale(scale.x(), scale.y(), scale.z());

    // Translate it back to the original center of the input mesh
    generatedMesh->translate(inputCenter);

    // Print the mesh statistcs
    if (options->writeStatistics)
        generatedMesh->printMeshStats("dmc", &options->outputPrefix);

    // Export it as a .OBJ mesh, and optionally an .OFF mesh
    std::string prefix = options->outputPrefix + "-dmc";
    generatedMesh->exportMesh(prefix, options->exportOBJ, options->exportOFF);

    // Optimize the mesh
    if (options->optimizeMesh)
    {
        // Set the prefix
        prefix = options->outputPrefix + OPTIMIZED_SUFFIX;

        //        // Optimize the mesh using the default parameters
        //        generatedMesh->optimize(options->smoothingIterations,
        //                                options->smoothingFactor);

        //        // Print the mesh statistcs
        //        if (options->writeStatistics)
        //            generatedMesh->printMeshStats("optimized", &options->outputPrefix);

        //        // Export the mesh
        //        generatedMesh->exportMesh(prefix, options->exportOBJ, options->exportOFF);


        // Set the prefix
        prefix = options->outputPrefix + OPTIMIZED_SUFFIX;

        // Optimization, if required
        for (int i = 0; i < 5; ++i)
        {
            // Coarse flat
            generatedMesh->coarseFlat(0.05);

            // Smooth
            // generatedMesh->smooth(15, 150, options->smoothingIterations);

            // Coarse dense
            generatedMesh->coarseDense(4.0);

            // Smooth
            // generatedMesh->smooth(15, 150, options->smoothingIterations);

            // Smooth normals
            generatedMesh->smoothNormals();
        }

        // Smooth
        generatedMesh->smooth(15, 150, options->smoothingIterations);

        // Smooth normals
        generatedMesh->smoothNormals();

        // Print the mesh statistcs
        if (options->writeStatistics)
            generatedMesh->printMeshStats("optimized", &options->outputPrefix);

        // Export the mesh
        generatedMesh->exportMesh(prefix,
                                  options->exportOBJ,
                                  options->exportOFF);

        // Fix self-intersections if any
        if (!options->ignoreSelfIntersections)
        {
            // Fix the mesh if it has any self intersections
            std::unique_ptr<AdvancedMesh> advancedMesh =
                    std::make_unique<AdvancedMesh>
                    (generatedMesh->getVertices(),
                     generatedMesh->getNumberVertices(),
                     generatedMesh->getTriangles(),
                     generatedMesh->getNumberTriangles());

            // Ensures that the mesh is truly two-advanced with no
            // self intersections
            advancedMesh->ensureWatertightness();

            // Print the mesh statistcs
            if (options->writeStatistics)
                advancedMesh->printMeshStats("watertight",
                                             &options->outputPrefix);

            // Export the repaired mesh
            std::string fileName = prefix + MANIFOLD_SUFFIX + OBJ_EXTENSION;
            advancedMesh->exportOBJ(fileName.c_str());
        }
    }

    // Free the DMC mesh
    generatedMesh->~Mesh();

    ULTRALISER_DONE;
}
}

int main(int argc , const char** argv)
{
    return Ultraliser::run(argc, argv);
}
