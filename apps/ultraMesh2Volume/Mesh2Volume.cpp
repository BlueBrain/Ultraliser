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
                "Sets the resolution of the volume based on the mesh "
                "dimensions.");
    args->addArgument(&autoResolution);

    Ultraliser::Argument voxelsPerMicron(
                "--voxels-per-micron",
                ARGUMENT_TYPE::INTEGER,
                "Number of voxels per micron in case --auto-resolution is used, "
                "default 3.",
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

    Ultraliser::Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Create a byte volume, where each voxel is stored in "
                "a single byte.");
    args->addArgument(&writeByteVolume);

    Ultraliser::Argument writeNRRDVolume(
                "--write-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Create an NRRD volume that is compatible with VTK.");
    args->addArgument(&writeNRRDVolume);

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

    Ultraliser::Argument reconstructMesh(
                "--reconstruct-mesh",
                ARGUMENT_TYPE::BOOL,
                "Reconstruct an output watertight mesh from the resulting "
                "volume and export it to .OBJ file.");
    args->addArgument(&reconstructMesh);

    Ultraliser::Argument optimizeMesh(
                "--optimize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh.");
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

    Ultraliser::Argument ignoreSelfIntersections(
                "--ignore-self-intersections",
                ARGUMENT_TYPE::BOOL,
                "Ignore, and take no action if the mesh has self intersecting "
                "faces. This process will speed up the generation of the final "
                "mesh, but the output mesh has no guarntees to be perfectly "
                "watertight.");
    args->addArgument(&ignoreSelfIntersections);

    Ultraliser::Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OBJ file.");
    args->addArgument(&exportOBJ);

    Ultraliser::Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OFF file.");
    args->addArgument(&exportOFF);

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
    options->autoResolution = args->getBoolValue(&autoResolution);
    options->voxelsPerMicron = args->getUnsignedIntegrValue(&voxelsPerMicron);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->writeNRRDVolume = args->getBoolValue(&writeNRRDVolume);
    options->solid = args->getBoolValue(&solid);
    options->volumeType = args->getStringValue(&volumeType);
    options->reconstructMesh = args->getBoolValue(&reconstructMesh);
    options->optimizeMesh = args->getBoolValue(&optimizeMesh);
    options->smoothingFactor = args->getFloatValue(&smoothingFactor);
    options->smoothingIterations = args->getUnsignedIntegrValue(&smoothingIterations);
    options->ignoreSelfIntersections = args->getBoolValue(&ignoreSelfIntersections);
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->projectXY = args->getBoolValue(&projectXY);
    options->projectZY = args->getBoolValue(&projectZY);
    options->stackXY = args->getBoolValue(&stackXY);
    options->stackZY = args->getBoolValue(&stackZY);
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

    if (!(options->writeBitVolume || options->writeByteVolume))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "volume: [--write-bit-volume, --write-byte-volume]");
    }

    if (options->optimizeMesh)
    {
        if (!options->reconstructMesh)
        {
            LOG_ERROR("The flag --reconstruct-mesh must be set to create a mesh "
                      "that we can optimize.");
        }
    }

    if (options->reconstructMesh)
    {
        if (!(options->exportOBJ || options->exportOFF))
        {
            LOG_ERROR("The user must specify at least one output format of the "
                      "mesh to export: [--export-obj, --export-off]");
        }
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
        options->prefix = Ultraliser::File::getName(options->inputMesh);
    }

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;

    LOG_TITLE("Ultralizing");
    LOG_STATUS("Output Directory [ %s ]", options->outputDirectory.c_str());

    // Return the executable options
    return options;
}

int main(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool reconstructs a volume from an input mesh.");

    // Parse the arguments and get the values
    Options* options = parseArguments(&args);

    // Load the mesh
    Ultraliser::Mesh* inputMesh = new Ultraliser::Mesh(options->inputMesh);

    // Write the statistics of the original mesh
    if (options->writeStatistics)
        inputMesh->printMeshStats("input", &options->outputPrefix);

    // Get relaxed bounding box to build the volume
    Ultraliser::Vector3f pMinInput, pMaxInput;
    inputMesh->computeBoundingBox(pMinInput, pMaxInput);


    // Extend the bounding box a little bit to avoid edge issues
    Ultraliser::Vector3f inputBB = pMaxInput - pMinInput;
    Ultraliser::Vector3f inputCenter = pMinInput + 0.5 * (pMaxInput - pMinInput);

    // Get the largest dimension
    float largestDimension = inputBB.x();
    if (inputBB.y() > largestDimension)
        largestDimension = inputBB.y();
    if (inputBB.z() > largestDimension)
        largestDimension = inputBB.z();

    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_WARNING("Volume resolution [%d], Largest dimension [%f]", resolution,
                largestDimension);

    Ultraliser::Volume* volume = new Ultraliser::Volume(
                pMinInput, pMaxInput, resolution,
                options->edgeGap,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Surface voxelization
    volume->surfaceVoxelization(inputMesh, true, false);

    // Free the mesh
    inputMesh->~Mesh();

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
        volume->printVolumeStats("Volume", &options->outputPrefix);

    // Write the stacks
    volume->writeStacks(options->outputDirectory, options->prefix,
                        options->stackXY, options->stackZY);

    if (options->reconstructMesh)
    {
        // Reconstruct a watertight mesh from the volume with DMC
        Ultraliser::DualMarchingCubes* dmc =
                new Ultraliser::DualMarchingCubes(volume);

        // Generate the mesh
        Ultraliser::Mesh* generatedMesh = dmc->generateMesh();

        // Free the volume
        volume->~Volume();

        // Center the reconstructed mesh at the origin
        generatedMesh->centerAtOrigin();

        // Compute the bounding box of the created mesh
        Ultraliser::Vector3f pMaxGenerated, pMinGenerated;
        generatedMesh->computeBoundingBox(pMinGenerated, pMaxGenerated);

        const Ultraliser::Vector3f generatedBB = pMaxGenerated - pMinGenerated;
        const Ultraliser::Vector3f scale = inputBB / generatedBB;

        // Scale the mesh
        generatedMesh->scale( scale.x(), scale.y(), scale.z());

        // Translate it back to the original center of the input mesh
        generatedMesh->translate(inputCenter);

        std::string prefix = options->outputPrefix + "-reconstructed";
        if (options->writeStatistics)
            generatedMesh->printMeshStats("DMC", &prefix);

        // Export it as a .OBJ mesh, and optionally an .OFF mesh
        generatedMesh->exportMesh(prefix, options->exportOBJ, options->exportOFF,
                                  false);

        // Optimize the mesh
        if (options->optimizeMesh)
        {
            // Set the prefix
            prefix = options->outputPrefix + OPTIMIZED_SUFFIX;

            // Optimize the mesh using the default parameters
            generatedMesh->optimize(options->smoothingIterations,
                              options->smoothingFactor);

            if (options->writeStatistics)
            {
                // Print the mesh statistcs
                generatedMesh->printMeshStats("Optimized", &prefix);
            }

            // Export it as a .OBJ mesh
            generatedMesh->exportMesh(prefix,
                                options->exportOBJ,
                                options->exportOFF, false);

            // Fix self0intersections if any
            if (!options->ignoreSelfIntersections)
            {
                // Fix the mesh if it has any self intersections
                std::unique_ptr<Ultraliser::AdvancedMesh> advancedMesh =
                        std::make_unique<Ultraliser::AdvancedMesh>
                        (generatedMesh->getVertices(),
                         generatedMesh->getNumberVertices(),
                         generatedMesh->getTriangles(),
                         generatedMesh->getNumberTriangles());

                // Ensures that the mesh is truly two-advanced with no
                // self intersections
                advancedMesh->ensureWatertightness();

                if (options->writeStatistics)
                {
                    // Print the mesh statistcs
                    advancedMesh->printMeshStats("Manifold", &prefix);
                }

                // Export the repaired mesh
                std::string fileName = prefix + MANIFOLD_SUFFIX + OBJ_EXTENSION;
                advancedMesh->exportOBJ(fileName.c_str());
            }
        }

        // Free the DMC mesh
        generatedMesh->~Mesh();
    }

    ULTRALISER_DONE;
}
