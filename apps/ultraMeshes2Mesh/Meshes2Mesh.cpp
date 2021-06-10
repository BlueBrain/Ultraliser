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
#include <AppArguments.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
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

    args->addInputMeshesDirectoryArguments();
    args->addOutputArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshesDirectoryArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyBoudsFileArgument();
    options->verifyMeshesPrefixArgument();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
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
    auto options = parseArguments(argc, argv);

    // A list of all the mesh files in the directory
    std::vector< std::string > meshFiles;

    // Query the input directory to see how many input meshes are there
    Ultraliser::Directory::locateMeshFiles(options->inputMeshesDirectory, meshFiles);

    // If the size of the list is zero, exit
    if(meshFiles.size() == 0)
        LOG_ERROR("No meshes were found in the input directory [ %s ]",
                  options->inputMeshesDirectory.c_str());

    // Compute the bounding box
    Vector3f pMaxInput, pMinInput;
    computeBoundingBoxForMeshes(options->boundsFile, options->inputMeshesDirectory,
                                meshFiles, pMaxInput, pMinInput);

    // Keep thr original values of the bounding box
    Vector3f inputBB = pMaxInput - pMinInput;
    Vector3f inputCenter = pMinInput + 0.5 * ( pMaxInput - pMinInput);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();

    uint64_t resolution;
    if (options->autoResolution)
        resolution = uint64_t(options->voxelsPerMicron * largestDimension);
    else
        resolution = options->volumeResolution;
    LOG_SUCCESS("Volume resolution [%d], Largest dimension [%f]", resolution, largestDimension);

    auto volume = new Volume(
                pMinInput, pMaxInput, resolution, options->edgeGap,
                Ultraliser::VolumeGrid::getType(options->volumeType));

    // Surface voxelization for all the mesh files
    volume->surfaceVoxelization(options->inputMeshesDirectory, meshFiles);

    // Enable solid voxelization
    if (options->useSolidVoxelization)
        volume->solidVoxelization(options->voxelizationAxis);

    // Generate the volume artifacts based on the given options
    generateVolumeArtifacts(volume, options);

    // Generate the reconstructed mesh from the marching cubes algorithm
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume);

    // Free the volume, it is not needed any further
    delete volume;

    // Generate the mesh artifacts
    generateMarchingCubesMeshArtifacts(mesh, options);

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
