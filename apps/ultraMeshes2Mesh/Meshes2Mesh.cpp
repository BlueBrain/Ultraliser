/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#include <Ultraliser.h>
#include <AppCommon.h>
#include <AppArguments.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tool generates an output mesh from a group of different input meshes."
              "The tool can naively append a group of meshes together into a single mesh object "
              "with multuple partitions if you enable the flag [--use-simple-mesh-join], or use "
              "volume reconstruction to reconstruct a mesh via volume reconstruction. The latter "
              "approach creates a high quality mesh that is guaranteed to be watertight even if all "
              "the input meshes are NOT watertights." );

    args->addInputMeshesDirectoryArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addMeshJoiningArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshesDirectoryArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyBoudsFileArgument();
    options->verifyMeshesPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void reconstructMeshWithJointOperation(const std::vector< std::string >& meshFiles,
                                       const AppOptions* options)
{
    TIMER_SET;

    // A new joint mesh that will be used to append all the meshes into a single mesh.
    Mesh* jointMesh = new Mesh();

    LOOP_STARTS("Merging Meshes")
    for (size_t i = 0; i < meshFiles.size(); ++i)
    {
        // Load the mesh
        std::string meshPath = options->inputMeshesDirectory + "/" + meshFiles[i];
        auto mesh = new Mesh(meshPath, false);

        // Joint
        jointMesh->append(mesh);

        // Release
        delete mesh;

        LOOP_PROGRESS(i, meshFiles.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_STATUS_IMPORTANT("Joining Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    // Export the mesh
    jointMesh->exportMesh(options->meshPrefix + JOINT_SUFFIX,
                          options->exportOBJ, options->exportPLY,
                          options->exportOFF, options->exportSTL);

    // Free the mesh
    delete jointMesh;
}

void reconstructMeshWithVolumeReconstruction(const std::vector< std::string >& meshFiles,
                                             const AppOptions* options)
{
    // Compute the bounding box
    Vector3f pMaxInput, pMinInput;
    computeBoundingBoxForMeshes(options->boundsFile, options->inputMeshesDirectory,
                                meshFiles, pMaxInput, pMinInput,
                                options->xScaleFactor, options->yScaleFactor, options->zScaleFactor);

    // Keep thr original values of the bounding box
    Vector3f inputBB = pMaxInput - pMinInput;

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    size_t resolution;
    if (options->scaledResolution)
        resolution = static_cast< size_t >(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_SUCCESS("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    auto volume = new Volume(
                pMinInput, pMaxInput, resolution, options->edgeGap,
                VolumeGrid::getType(options->volumeType));

    // Surface voxelization for all the mesh files
    volume->surfaceVoxelization(options->inputMeshesDirectory, meshFiles);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Extract the mesh from the volume again
    auto mesh = reconstructMeshFromVolume(volume, options);

    // Free the volume, it is not needed any further
    delete volume;

    // Generate the mesh artifacts
    generateReconstructedMeshArtifacts(mesh, options);

    // Free
    delete mesh;
    delete options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc, argv);

    // A list of all the mesh files in the directory
    std::vector< std::string > meshFiles;

    // Query the input directory to see how many input meshes are there
    Directory::locateMeshFiles(options->inputMeshesDirectory, meshFiles);

    // If the size of the list is zero, exit
    if(meshFiles.size() == 0)
    {
        LOG_ERROR("No meshes were found in the input directory [ %s ]",
                  options->inputMeshesDirectory.c_str());
    }

    if (options->simpleMeshJoin)
        reconstructMeshWithJointOperation(meshFiles, options);
    else
        reconstructMeshWithVolumeReconstruction(meshFiles, options);
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
