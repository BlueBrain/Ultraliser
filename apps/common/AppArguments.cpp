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
#include "AppArguments.h"

namespace Ultraliser
{

AppArguments::AppArguments(const int& argc , const char** argv, const std::string& help)
{
    _args = new Args(argc, argv, help);

    // Please fill the _options in the corresponding directories
    _options = new Options();
}

void AppArguments::addInputMeshArguments()
{
    Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the input mesh.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputMesh);
    _options->inputMesh = _args->getStringValue(&inputMesh);
}

void AppArguments::addInputMeshesDirectoryArguments()
{
    Argument inputMeshesDirectory(
                "--input-directory",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the directory where the input meshes are located.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputMeshesDirectory);
    _options->inputMeshesDirectory = _args->getStringValue(&inputMeshesDirectory);
}

void AppArguments::addInputMorphologyArguments()
{
    Argument inputMorphology(
                "--morphology",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the input morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputMorphology);
    _options->inputMorphology = _args->getStringValue(&inputMorphology);
}

void AppArguments::addInputMaskDirectoryArguments()
{
    Argument maskDirectory(
                "--mask-directory",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the directory where the mask stack is located.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&maskDirectory);
    _options->inputMaskDirectory = _args->getStringValue(&maskDirectory);
}

void AppArguments::addMaskArguments()
{
    Argument maskWidth(
                "--mask-width",
                ARGUMENT_TYPE::INTEGER,
                "The width of the mask.",
                ARGUMENT_PRESENCE::MANDATORY,
                "0");
    _args->addArgument(&maskWidth);
    _options->maskWidth = _args->getIntegrValue(&maskWidth);

    Argument maskHeight(
                "--mask-height",
                ARGUMENT_TYPE::INTEGER,
                "The height of the mask.",
                ARGUMENT_PRESENCE::MANDATORY,
                "0");
    _args->addArgument(&maskHeight);
    _options->maskHeight = _args->getIntegrValue(&maskWidth);
}

void AppArguments::addInputVolumeArguments()
{
    Argument inputVolume(
                "--volume",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the volume.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputVolume);
    _options->inputVolume = _args->getStringValue(&inputVolume);
}

void AppArguments::addInputVolumeParametersArguments()
{
    Argument isoValue(
                "--iso-value",
                ARGUMENT_TYPE::INTEGER,
                "The iso value where the volume will get segmented. Default 127.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "127");
    _args->addArgument(&isoValue);
    _options->isoValue = _args->getIntegrValue(&isoValue);

    Argument fullRangeIsoValue(
                "--full-range-iso-value",
                ARGUMENT_TYPE::BOOL,
                "If the voxel contains any value other than zero, then use it."
                "If this option is set the --iso-value option is ignored.");
    _args->addArgument(&fullRangeIsoValue);
    _options->fullRangeIsoValue = _args->getIntegrValue(&fullRangeIsoValue);

    Argument writeHistogram(
                "--write-histogram",
                ARGUMENT_TYPE::BOOL,
                "Write the histogram of the volume into a text file.");
    _args->addArgument(&writeHistogram);
    _options->writeHistogram = _args->getBoolValue(&isoValue);

    Argument zeroPaddingVoxels(
                "--zero-paddgin-voxels",
                ARGUMENT_TYPE::INTEGER,
                "The number of zero-padding voxels that will be appended to "
                "the volume to avoid any clipping artifacts, default 0",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0");
    _args->addArgument(&zeroPaddingVoxels);
    _options->zeroPaddingVoxels = _args->getIntegrValue(&zeroPaddingVoxels);
}

void AppArguments::addOutputArguments()
{
    Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the directory where the results (artifacts) will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&outputDirectory);
    _options->outputDirectory = _args->getStringValue(&outputDirectory);

    Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "A prefix that will be used to label the output files. "
                "If this is not given, the name of the input file/directory will be considered.");
    _args->addArgument(&prefix);
    _options->prefix = _args->getStringValue(&prefix);
}

void AppArguments::addVoxelizationArguments()
{
    Argument boundsFile(
                "--bounds-file",
                ARGUMENT_TYPE::STRING,
                "A file that defines the bounding box or ROI that will be voxelized and meshed."
                "This option is used to select a specifc region of interest from the space to voxelize.");
    _args->addArgument(&boundsFile);
    _options->boundsFile = _args->getStringValue(&boundsFile);

    Argument edgeGap(
                "--edge-gap",
                ARGUMENT_TYPE::FLOAT,
                "Some little extra space to avoid edges intersection. Default 0.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.0");
    _args->addArgument(&edgeGap);
    _options->edgeGap = _args->getFloatValue(&edgeGap);

    Argument volumeResolution(
                "--resolution",
                ARGUMENT_TYPE::INTEGER,
                "The base resolution of the volume. Default 512."
                "This resolution is set to the larget dimension of the bounding box of the input "
                "dataset, and the resolution of the other dimensions are computed accordingly.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "512");
    _args->addArgument(&volumeResolution);
    _options->volumeResolution = _args->getUnsignedIntegrValue(&volumeResolution);

    Argument autoResolution(
                "--auto-resolution",
                ARGUMENT_TYPE::BOOL,
                "Sets the resolution of the volume based on the mesh dimensions.");
    _args->addArgument(&autoResolution);
    _options->autoResolution = _args->getBoolValue(&autoResolution);

    Argument voxelsPerMicron(
                "--voxels-per-micron",
                ARGUMENT_TYPE::INTEGER,
                "Number of voxels per micron in case --auto-resolution is used. Default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
    _args->addArgument(&voxelsPerMicron);
    _options->voxelsPerMicron = _args->getUnsignedIntegrValue(&voxelsPerMicron);

    Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume format to perform the voxelization: [bit, byte, voxel]. "
                "By default, it is a bit volume to reduce the memory foot print.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    _args->addArgument(&volumeType);
    _options->volumeType = _args->getStringValue(&volumeType);

    Argument useSolidVoxelization(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Use solid voxelization to fill the interior of the volume shell.");
    _args->addArgument(&useSolidVoxelization);
    _options->useSolidVoxelization = _args->getBoolValue(&useSolidVoxelization);

    Argument voxelizationAxis(
                "--vozelization-axis",
                ARGUMENT_TYPE::STRING,
                "The axis where the solid voxelization operation will be performed. "
                "Use one of the following options [x, y, z, or xyz]. "
                "If you use x or y or z the voxelization will happen along a single axis, "
                "otherwise, using xyz will perform the solid voxelization along the three main "
                "axes of the volume to avoid filling any loops in the morphology."
                "By default, the Z-axis solid voxelization with xyz is applied if the --solid "
                "flag is set.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "z");
    _args->addArgument(&voxelizationAxis);
    _options->voxelizationAxis = Volume::getSolidvoxelizationAxis(
                _args->getStringValue(&voxelizationAxis));
}

void AppArguments::addVolumeProjectionArguments()
{
    Argument projectXY(
                "--project-xy",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the Z-axis into a gray-scale image.");
    _args->addArgument(&projectXY);
    _options->projectXY = _args->getBoolValue(&projectXY);

    Argument projectXZ(
                "--project-xz",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the Y-axis into a gray-scale image.");
    _args->addArgument(&projectXZ);
    _options->projectXZ = _args->getBoolValue(&projectXZ);

    Argument projectZY(
                "--project-zy",
                ARGUMENT_TYPE::BOOL,
                "Project the volume along the X-axis into a gray-scale image.");
    _args->addArgument(&projectZY);
    _options->projectZY = _args->getBoolValue(&projectZY);

    Argument projectColorCoded(
                "--project-color-coded",
                ARGUMENT_TYPE::BOOL,
                "Generate color-coded projections of the volume to help debugging it.");
    _args->addArgument(&projectColorCoded);
    _options->projectColorCoded = _args->getBoolValue(&projectColorCoded);
}

void AppArguments::addVolumeExportArguments()
{
    Argument writeBitVolume(
                "--write-bit-volume",
                ARGUMENT_TYPE::BOOL,
                "Generate a bit volume, where each voxel is stored in a single bit.");
    _args->addArgument(&writeBitVolume);
    _options->writeBitVolume = _args->getBoolValue(&writeBitVolume);

    Argument writeByteVolume(
                "--write-byte-volume",
                ARGUMENT_TYPE::BOOL,
                "Generate a byte volume, where each voxel is stored in a single byte.");
    _args->addArgument(&writeByteVolume);
    _options->writeByteVolume = _args->getBoolValue(&writeByteVolume);

    Argument writeNRRDVolume(
                "--write-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Generate an NRRD volume that is compatible with VTK.");
    _args->addArgument(&writeNRRDVolume);
    _options->writeNRRDVolume = _args->getBoolValue(&writeNRRDVolume);

    Argument exportVolumeMesh(
                "--export-volume-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export a mesh that represents the volume where each voxel will be a cube.");
    _args->addArgument(&exportVolumeMesh);
    _options->exportVolumeMesh = _args->getBoolValue(&exportVolumeMesh);

    Argument exportVolumeBoundingBoxMesh(
                "--export-volume-bounding-box-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export a mesh that represents the bounding box of the volume."
                "This mesh is primarily used for debugging purposes.");
    _args->addArgument(&exportVolumeBoundingBoxMesh);
    _options->exportVolumeBoundingBoxMesh = _args->getBoolValue(&exportVolumeBoundingBoxMesh);

    Argument exportVolumeGridMesh(
                "--export-volume-grid-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export a mesh that represents the volumetric grid used to voxelize the mesh."
                "This mesh is primarily used for debugging purposes.");
    _args->addArgument(&exportVolumeGridMesh);
    _options->exportVolumeGridMesh = _args->getBoolValue(&exportVolumeGridMesh);
}

void AppArguments::addStacksArguments()
{
    Argument stackXY(
                "--stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the Z-axis of the volume.");
    _args->addArgument(&stackXY);
    _options->stackXY = _args->getBoolValue(&stackXY);


    Argument stackXZ(
                "--stack-xz",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the Y-axis of the volume.");
    _args->addArgument(&stackXZ);
    _options->stackXZ = _args->getBoolValue(&stackXZ);

    Argument stackZY(
                "--stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the X-axis of the volume.");
    _args->addArgument(&stackZY);
    _options->stackZY = _args->getBoolValue(&stackZY);
}

void AppArguments::addMeshVoxelizationArgument()
{
    Argument voxelizeMesh(
                "--voxelize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Voxelize the given mesh and produce its artifacts.");
    _args->addArgument(&voxelizeMesh);
    _options->voxelizeMesh = _args->getBoolValue(&voxelizeMesh);
}

void AppArguments::addVolumeArguments()
{
    addVoxelizationArguments();
    addVolumeProjectionArguments();
    addStacksArguments();
    addVolumeExportArguments();
}

void AppArguments::addMeshOptimizationArguments()
{
    Argument preservePartitions(
                "--preserve-partitions",
                ARGUMENT_TYPE::BOOL,
                "Keeps all the partitions of the mesh if the input mesh contains more than one.");
    _args->addArgument(&preservePartitions);
    _options->preservePartitions = _args->getBoolValue(&preservePartitions);

    Argument optimizeMesh(
                "--optimize-mesh",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh using the default optimization strategy.");
    _args->addArgument(&optimizeMesh);
    _options->optimizeMeshHomogenous = _args->getBoolValue(&optimizeMesh);

    Argument adaptiveOptimization(
                "--adaptive-optimization",
                ARGUMENT_TYPE::BOOL,
                "Optimize the reconstructed mesh using the adaptive optimization strategy.");
    _args->addArgument(&adaptiveOptimization);
    _options->optimizeMeshAdaptively = _args->getBoolValue(&adaptiveOptimization);

    Argument optimizationIterations(
                "--optimization-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to optimize the resulting mesh. Default value 1. "
                "If this value is set to 0, the optimization process will be ignored.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
    _args->addArgument(&optimizationIterations);
    _options->optimizationIterations = _args->getUnsignedIntegrValue(&optimizationIterations);

    Argument smoothingIterations(
                "--smooth-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh, Default 5.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5");
    _args->addArgument(&smoothingIterations);
    _options->smoothingIterations = _args->getUnsignedIntegrValue(&smoothingIterations);

    Argument flatFactor(
                "--flat-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseFlat function. Default 0.05.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
    _args->addArgument(&flatFactor);
    _options->flatFactor = _args->getFloatValue(&flatFactor);

    Argument denseFactor(
                "--dense-factor",
                ARGUMENT_TYPE::FLOAT,
                "A factor that is used for the coarseDense function. Default 5.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "5.0");
    _args->addArgument(&denseFactor);
    _options->denseFactor = _args->getFloatValue(&denseFactor);
}

void AppArguments::addLaplacianOperatorArguments()
{
    Argument ignoreLaplacianSmoothing(
                "--ignore-laplacian-filter",
                ARGUMENT_TYPE::BOOL,
                "Ignore Laplacian smoothing.");
    _args->addArgument(&ignoreLaplacianSmoothing);
    _options->ignoreLaplacianSmoothing = _args->getBoolValue(&ignoreLaplacianSmoothing);

    Argument laplacianIterations(
                "--laplacian-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh with Laplacian filter. "
                "Default 3.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "3");
    _args->addArgument(&laplacianIterations);
    _options->laplacianIterations = _args->getIntegrValue(&laplacianIterations);
}

void AppArguments::addMeshExportArguments()
{
    Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the reconstructed mesh to an .OBJ file.");
    _args->addArgument(&exportOBJ);
    _options->exportOBJ = _args->getBoolValue(&exportOBJ);

    Argument exportPLY(
                "--export-ply",
                ARGUMENT_TYPE::BOOL,
                "Export the reconstructed mesh to an .PLY file.");
    _args->addArgument(&exportPLY);
    _options->exportPLY = _args->getBoolValue(&exportPLY);


    Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the reconstructed mesh to an .OFF file.");
    _args->addArgument(&exportOFF);
    _options->exportOFF = _args->getBoolValue(&exportOFF);

    Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the reconstructed mesh to an .STL file.");
    _args->addArgument(&exportSTL);
    _options->exportSTL = _args->getBoolValue(&exportSTL);
}

void AppArguments::addMeshScaleArguments()
{
    Argument xScaleFactor(
                "--x-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the X-axis, Default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    _args->addArgument(&xScaleFactor);
    _options->xScaleFactor = _args->getIntegrValue(&xScaleFactor);

    Argument yScaleFactor(
                "--y-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Y-axis. Default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    _args->addArgument(&yScaleFactor);
    _options->yScaleFactor = _args->getIntegrValue(&yScaleFactor);

    Argument zScaleFactor(
                "--z-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Z-axis. Default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    _args->addArgument(&zScaleFactor);
    _options->zScaleFactor = _args->getIntegrValue(&zScaleFactor);
}

void AppArguments::addMeshArguments()
{
    addMeshOptimizationArguments();
    addLaplacianOperatorArguments();
    addMeshScaleArguments();
    addMeshExportArguments();
}

void AppArguments::addDataArguments()
{
    Argument writeStatistics(
                "--stats",
                ARGUMENT_TYPE::BOOL,
                "Write the statistics of the resulting meshes/volumes/morphologies.");
    _args->addArgument(&writeStatistics);
    _options->writeStatistics = _args->getBoolValue(&writeStatistics);

    Argument writeDistributions(
                "--dists",
                ARGUMENT_TYPE::BOOL,
                "Write the distributions of the resulting meshes/volumes/morphologies.");
    _args->addArgument(&writeDistributions);
    _options->writeDistributions = _args->getBoolValue(&writeDistributions);
}

void AppArguments::addSuppressionArguments()
{
    Argument ignoreSelfIntersections(
                "--ignore-self-intersections",
                ARGUMENT_TYPE::BOOL,
                "Ignore, and take no action if the mesh has self intersecting faces. This process "
                "will speed up the generation of the final mesh, but the output mesh has no "
                "guarntees to be perfectly watertight.");
    _args->addArgument(&ignoreSelfIntersections);
    _options->ignoreSelfIntersections = _args->getBoolValue(&ignoreSelfIntersections);

    Argument writeMarchingCubeMesh(
                "--write-marching-cubes-mesh",
                ARGUMENT_TYPE::BOOL,
                "Write the resulting mesh (and its artifacts) from the marching cubes operation to "
                "debug it.");
    _args->addArgument(&writeMarchingCubeMesh);
    _options->writeMarchingCubeMesh = _args->getBoolValue(&writeMarchingCubeMesh);

    Argument ignoreOptimizedNonWatertightMesh(
                "--ignore-optimized-non-watertight-mesh",
                ARGUMENT_TYPE::BOOL,
                "Ignore the resulting mesh from the optimization process without removing self "
                "intersections.");
    _args->addArgument(&ignoreOptimizedNonWatertightMesh);
    _options->ignoreOptimizedNonWatertightMesh =
            _args->getBoolValue(&ignoreOptimizedNonWatertightMesh);
}

Options* AppArguments::getOptions()
{
    // Parse the arguments
    _args->parse();

    return _options;
}

}

