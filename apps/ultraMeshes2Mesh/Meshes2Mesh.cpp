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
    Ultraliser::Argument inputDirectory(
                "--input-directory",
                ARGUMENT_TYPE::STRING,
                "Input directory where the input meshes are located.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputDirectory);

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
                "Sets the resolution of the volume based on the mesh dimensions.");
    args->addArgument(&autoResolution);

    Ultraliser::Argument edgeGap(
                "--edge-gap",
                ARGUMENT_TYPE::FLOAT,
                "Some little extra space to avoid edges intersection, "
                "default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
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
    options->inputDirectory = args->getStringValue(&inputDirectory);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->boundsFile = args->getStringValue(&boundsFile);
    options->volumeResolution = args->getUnsignedIntegrValue(&volumeResolution);
    options->autoResolution = args->getBoolValue(&autoResolution);
    options->edgeGap = args->getFloatValue(&edgeGap);
    options->writeBitVolume = args->getBoolValue(&writeBitVolume);
    options->writeByteVolume = args->getBoolValue(&writeByteVolume);
    options->solid = args->getBoolValue(&solid);
    options->volumeType = args->getStringValue(&volumeType);
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
    if (!Ultraliser::Directory::exists(options->inputDirectory))
        LOG_ERROR("The directory [ %s ] does NOT exist! ",
                  options->inputDirectory.c_str());

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
        LOG_ERROR("The directory [ %s ] does NOT exist!",
                  options->outputDirectory.c_str());

    if (!(options->exportOBJ || options->exportOFF))
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-off]");

    if (options->boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the volume will be computed "
                    "on the fly");
        options->boundsFile = EMPTY;
    }
    else
        LOG_WARNING("The bounding box of the volume will be loaded from [ %s ]",
                    options->boundsFile.c_str());

    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
        options->prefix =
                Ultraliser::Directory::getName(options->inputDirectory);

    // Construct the output prefix
    options->outputPrefix = options->outputDirectory + "/" + options->prefix;

    LOG_TITLE("Ultralizing");
    LOG_STATUS_IMPORTANT("Output Directory [ %s ]\n",
                         options->outputDirectory.c_str());

    // Return the executable options
    return options;
}

void computeBoundingBox(const Options* options,
                        std::vector< std::string > meshFiles,
                        Ultraliser::Vector3f& pMax, Ultraliser::Vector3f& pMin )
{
    if (options->boundsFile == EMPTY )
    {
        LOG_STATUS("Computing Bounding Box" );
        TIMER_SET;

        // Vectors containing all the pMin and pMax of all the objects
        std::vector< Ultraliser::Vector3f > pMinVector, pMaxVector;

        // Resize to perform in parallel
        pMinVector.resize(meshFiles.size());
        pMaxVector.resize(meshFiles.size());

        LOOP_STARTS("Loading Meshes");
        size_t loadedMeshCount = 0;
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp parallel for
#endif
        for( size_t iMesh = 0; iMesh < meshFiles.size(); iMesh++ )
        {
            // Create and load the mesh from the file
            std::string meshName= meshFiles[iMesh];
            std::string meshFile = meshName;
            if (options->inputDirectory != EMPTY )
                meshFile = options->inputDirectory + "/" + meshName;

            if (Ultraliser::File::exists(meshFile))
            {

#ifdef ULTRALISER_USE_OPENMP
                #pragma omp atomic
                ++loadedMeshCount;
                if (omp_get_thread_num() == 0)
                    LOOP_PROGRESS(loadedMeshCount,  meshFiles.size());
#else
                ++loadedMeshCount;
                LOOP_PROGRESS(loadedMeshCount,  meshFiles.size());
#endif

                // Load the mesh
                Ultraliser::Mesh* mesh = new Ultraliser::Mesh (meshFile);

                // Compute its bounding box
                Ultraliser::Vector3f pMinMesh, pMaxMesh;
                mesh->computeBoundingBox(pMinMesh, pMaxMesh);
                pMinVector[iMesh] = pMinMesh;
                pMaxVector[iMesh] = pMaxMesh;
                mesh->~Mesh();
            }
            else
            {
                LOG_WARNING("Ignoring Mesh: [ %s ]", meshFile.c_str());
            }
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);

        if (loadedMeshCount == 0 )
            LOG_ERROR("No Loaded Meshes");
        else
            LOG_DETAIL("Loaded Meshes: [%zu/%zu]",
                                 loadedMeshCount, meshFiles.size());


        // Compute the bounding box of the group
        pMax.x() = std::numeric_limits< float >::min();
        pMax.y() = std::numeric_limits< float >::min();
        pMax.z() = std::numeric_limits< float >::min();

        pMin.x() = std::numeric_limits< float >::max();
        pMin.z() = std::numeric_limits< float >::max();
        pMin.y() = std::numeric_limits< float >::max();

        LOOP_STARTS("Computing Bounding Box");
        TIMER_RESET;
        for(uint64_t iMesh = 0; iMesh < meshFiles.size(); ++iMesh)
        {
            LOOP_PROGRESS(iMesh,  meshFiles.size());

            Ultraliser::Vector3f pMinObject = pMinVector[ iMesh ];
            Ultraliser::Vector3f pMaxObject = pMaxVector[ iMesh ];

            if (pMinObject.x() < pMin.x()) pMin.x() = pMinObject.x();
            if (pMinObject.y() < pMin.y()) pMin.y() = pMinObject.y();
            if (pMinObject.z() < pMin.z()) pMin.z() = pMinObject.z();

            if (pMaxObject.x() > pMax.x()) pMax.x() = pMaxObject.x();
            if (pMaxObject.y() > pMax.y()) pMax.y() = pMaxObject.y();
            if (pMaxObject.z() > pMax.z()) pMax.z() = pMaxObject.z();
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOG_STATUS_IMPORTANT("Loading Bounding Box from [ %s ]",
                             options->boundsFile.c_str());

        // Verify the bounding box file
        if (Ultraliser::File::exists(options->boundsFile))
            Ultraliser::File::parseBoundsFile(options->boundsFile, pMin, pMax);
        else
            LOG_ERROR("No Bounding Box File is Provided !");
    }
}

int main(int argc , const char** argv)
{
    // Arguments
    Args args(argc, argv,
              "This tool reconstructs a watertight polygonal mesh from an set of "
              "input non-watertight meshes. The generated mesh can be also "
              "optimized to reduce the number of triangles while preserving "
              "the volume. "
              "The output mesh is guaranteed in all cases to be two-advanced "
              "with no self-intersecting faces unless the "
              "--ignore-self-intersections flag is enabled."
              "The input meshes could only be of .obj or .ply extensions, "
              "otherwise, will be ignored during the loading."
              "If the extent bounding box is not given via --bounds-file, the "
              "aggregate bounding box of all the loaded meshes will be computed "
              "on the fly and used AS IS." );

    // Parse the arguments and get the values
    Options* options = parseArguments(&args);

    // A list of all the mesh files in the directory
    std::vector< std::string > meshFiles;

    // Query the input directory to see how many input meshes are there
    Ultraliser::Directory::locateMeshFiles(options->inputDirectory,
                                           meshFiles);

    // If the size of the list is zero, exit
    if(meshFiles.size() == 0)
        LOG_ERROR("No meshes were found in the input directory [ %s ]",
                  options->inputDirectory.c_str());

    // Compute the bounding box
    Ultraliser::Vector3f pMaxInput, pMinInput;
    computeBoundingBox(options, meshFiles, pMaxInput, pMinInput);

    // Keep thr original values of the bounding box
    Ultraliser::Vector3f inputBB = pMaxInput - pMinInput;
    Ultraliser::Vector3f inputCenter = pMinInput + 0.5 * ( pMaxInput - pMinInput);

    Ultraliser::Volume* volume = new Ultraliser::Volume(
                pMinInput, pMaxInput, options->volumeResolution, 10,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Surface voxelization
    volume->surfaceVoxelization(options->inputDirectory, meshFiles);

    // Enable solid voxelization
    if (options->solid)
        volume->solidVoxelization();

    // Projecting the volume to validate its content
    volume->project(options->outputPrefix,
                    options->projectXY, options->projectZY);

    // Write the volume
    volume->writeVolumes(options->outputPrefix,
                         options->writeBitVolume,
                         options->writeByteVolume);

    // Print the volume statistics
    if (options->writeStatistics)
        volume->printVolumeStats("Volume", &options->outputPrefix);

    // Write the stacks
    volume->writeStacks(options->outputDirectory, options->prefix,
                        options->stackXY, options->stackZY);

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
    {
        // Print the mesh statistcs
        generatedMesh->printMeshStats("DMC", &prefix);
    }

    // Export it as a .OBJ mesh, and optionally an .OFF mesh
    generatedMesh->exportMesh(prefix,
                              options->exportOBJ,
                              options->exportOFF,
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

            // Remove the reconstructed mesh
            generatedMesh->~Mesh();

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
            advancedMesh->exportOBJ(fileName.c_str()) ;
        }
    }

    // Free the DMC mesh
    generatedMesh->~Mesh();

    ULTRALISER_DONE;
}
