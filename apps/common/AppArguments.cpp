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
#include "AppArguments.h"

namespace Ultraliser
{

AppArguments::AppArguments(const int& argc , const char** argv, const std::string& help)
{
    _args = new Args(argc, argv, help);

    // Please fill the _options in the corresponding directories
    _options = new AppOptions();
}

void AppArguments::addInputMeshArguments()
{
    Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the input mesh. "
                "Supported mesh types: .OBJ, .PLY, .STL, .OFF, .H5"
                "The .H5 meshes are specified by the MICrONS Explorer.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputMesh);
    _options->inputMeshPath = _args->getStringValue(&inputMesh);
}

void AppArguments::addTargetMeshArguments()
{
    Argument targetMesh(
                "--target-mesh",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the target mesh, with which the input mesh will be mapped to",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&targetMesh);
    _options->targetMeshPath = _args->getStringValue(&targetMesh);
}

void AppArguments::addInputMeshesDirectoryArguments()
{
    Argument inputMeshesDirectory(
                "--input-directory",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the directory where a group of input meshes are located to "
                "be loaded in Ultraliser to create an alternative output watertight mesh.",
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
    _options->inputMorphologyPath = _args->getStringValue(&inputMorphology);
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
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&maskWidth);
    _options->maskWidth = _args->getIntegrValue(&maskWidth);

    Argument maskHeight(
                "--mask-height",
                ARGUMENT_TYPE::INTEGER,
                "The height of the mask.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&maskHeight);
    _options->maskHeight = _args->getIntegrValue(&maskHeight);
}

void AppArguments::addMorphologyExtractionArguments()
{
    Argument bboxWidth(
                "--bbox-width",
                ARGUMENT_TYPE::FLOAT,
                "The width of bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxWidth);
    _options->bboxWidth = _args->getFloatValue(&bboxWidth);

    Argument bboxHeight(
                "--bbox-height",
                ARGUMENT_TYPE::FLOAT,
                "The height of bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxHeight);
    _options->bboxHeight = _args->getFloatValue(&bboxHeight);

    Argument bboxDepth(
                "--bbox-depth",
                ARGUMENT_TYPE::FLOAT,
                "The depth of bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxDepth);
    _options->bboxDepth = _args->getFloatValue(&bboxDepth);

    Argument bboxCenterX(
                "--bbox-center-x",
                ARGUMENT_TYPE::FLOAT,
                "The X-center of the bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxCenterX);
    _options->bboxCenterX = _args->getFloatValue(&bboxCenterX);

    Argument bboxCenterY(
                "--bbox-center-y",
                ARGUMENT_TYPE::FLOAT,
                "The Y-center of the bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxCenterY);
    _options->bboxCenterY = _args->getFloatValue(&bboxCenterY);

    Argument bboxCenterZ(
                "--bbox-center-z",
                ARGUMENT_TYPE::FLOAT,
                "The Y-center of the bounding box of the region to be extracted from the morphology.",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&bboxCenterZ);
    _options->bboxCenterZ = _args->getFloatValue(&bboxCenterZ);
}

void AppArguments::addInputVolumeArguments()
{
    Argument inputVolume(
                "--volume",
                ARGUMENT_TYPE::STRING,
                "The absolute path to the volume (prefix).",
                ARGUMENT_PRESENCE::MANDATORY);
    _args->addArgument(&inputVolume);
    _options->inputVolumePath = _args->getStringValue(&inputVolume);
}

void AppArguments::addInputVolumeParametersArguments()
{
    std::string info = "The --iso-option can be one of the following options: ";
    info += "[(value), values, min, max, range, nonzero] ";
    info += "* If --iso-option=value is selected, the value is defined by the option "
            "--iso-value=[VALUE] ";
    info += "* If --iso-option=values is selected, the values are loaded from a given ASCII file "
            "as a space-separated list of values. The ASCII file is given using the option "
            "--iso-values-file=[ISO_VALUES_FILE]. Note that the values don't have to be a "
            "continuous range.";
    info += "* If --iso-option=min is selected, all the range above or equal to a given value will "
            "be used to segment the volume. The minimum value is defined by the option "
            "--min-value=[VALUE].";
    info += "* If --iso-option=max is selected, all the range below or equal to a given value will "
            "be used to segment the volume. The maximum value is defined by the option "
            "--max-value=[VALUE].";
    info += "* If the --iso-option=range is selected, then all the values between two values will "
            "be used to segment the volume. The minimum and maximum values will be defined by the "
            "options --min-value=[VALUE] and --max-value=[VALUE].";
    info += "* If the --iso-option=nonzero is selected, then all the non-zero voxels will be used "
            "to calculate the iso-volume.";

    Argument isoOption(
                "--iso-option",
                ARGUMENT_TYPE::STRING,
                info,
                ARGUMENT_PRESENCE::OPTIONAL,
                "value");
    _args->addArgument(&isoOption);
    _options->isoOption = _args->getStringValue(&isoOption);

    Argument isoValue(
                "--iso-value",
                ARGUMENT_TYPE::INTEGER,
                "The iso-value where the volume will get segmented. Default 128.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "128");
    _args->addArgument(&isoValue);
    _options->isoValue = _args->getIntegrValue(&isoValue);

    Argument minIsoValue(
                "--min-value",
                ARGUMENT_TYPE::INTEGER,
                "The minimum value used to segment the volume. Default 1.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1");
    _args->addArgument(&minIsoValue);
    _options->minIsoValue = _args->getIntegrValue(&minIsoValue);

    Argument maxIsoValue(
                "--max-value",
                ARGUMENT_TYPE::INTEGER,
                "The maximum value used to segment the volume. Default 255.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "255");
    _args->addArgument(&maxIsoValue);
    _options->maxIsoValue = _args->getIntegrValue(&maxIsoValue);

    Argument isovaluesFile(
                "--isovalues-file",
                ARGUMENT_TYPE::STRING,
                "A file containing a list of values to extract the volume.",
                ARGUMENT_PRESENCE::OPTIONAL);
    _args->addArgument(&isovaluesFile);
    _options->isovaluesFile = _args->getStringValue(&isovaluesFile);

    Argument writeHistogram(
                "--write-histogram",
                ARGUMENT_TYPE::BOOL,
                "Write the histogram of the volume into a text file.");
    _args->addArgument(&writeHistogram);
    _options->writeHistogram = _args->getBoolValue(&writeHistogram);
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
                "If this is not given, the name of the input file/directory will be considered.",
                ARGUMENT_PRESENCE::OPTIONAL);
    _args->addArgument(&prefix);
    _options->prefix = _args->getStringValue(&prefix);
}

void AppArguments::addAutoParametersArguments()
{
    Argument autoParameters(
                "--auto-parameters",
                ARGUMENT_TYPE::BOOL,
                "Use the parameters that are computed automatically based on the statistical "
                "analysis of the input data.");
    _args->addArgument(&autoParameters);
    _options->autoParameters = _args->getBoolValue(&autoParameters);
}

void AppArguments::addSolidVoxelizationArguments()
{
    Argument useSolidVoxelization(
                "--solid",
                ARGUMENT_TYPE::BOOL,
                "Use solid voxelization to fill the interior of the volume shell.");
    _args->addArgument(&useSolidVoxelization);
    _options->useSolidVoxelization = _args->getBoolValue(&useSolidVoxelization);

    Argument voxelizationAxis(
                "--voxelization-axis",
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
                "Some little extra space to avoid edges intersection. Default 0.05.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
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

    Argument scaledResolution(
                "--scaled-resolution",
                ARGUMENT_TYPE::BOOL,
                "Sets the resolution of the volume based on the mesh dimensions.");
    _args->addArgument(&scaledResolution);
    _options->scaledResolution = _args->getBoolValue(&scaledResolution);

    Argument voxelsPerMicron(
                "--voxels-per-micron",
                ARGUMENT_TYPE::FLOAT,
                "Number of voxels per micron. If this flag is set, then resolution will be ignored.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0");
    _args->addArgument(&voxelsPerMicron);
    _options->voxelsPerMicron = _args->getFloatValue(&voxelsPerMicron);

    Argument volumeType(
                "--volume-type",
                ARGUMENT_TYPE::STRING,
                "Specify a volume format to perform the voxelization: [bit, byte, voxel]. "
                "By default, it is a bit volume to reduce the memory foot print.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "bit");
    _args->addArgument(&volumeType);
    _options->volumeType = _args->getStringValue(&volumeType);

    addSolidVoxelizationArguments();
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
    Argument exportBitVolume(
                "--export-bit-volume",
                ARGUMENT_TYPE::BOOL,
                "Export an Ultraliser-specific bit volume, where each voxel is stored in 1 bit. "
                "The header and data are stored in a single file with the extention .vol.");
    _args->addArgument(&exportBitVolume);
    _options->exportBitVolume = _args->getBoolValue(&exportBitVolume);

    Argument exportUnsignedVolume(
                "--export-unsigned-volume",
                ARGUMENT_TYPE::BOOL,
                "Export an Ultraliser-specific unsigned volume, where each voxel is stored either "
                "in 1, 2, 3 or 4 bytes depending on the type of the volume. "
                "The data and header are stored in a single file with the extention .vol.");
    _args->addArgument(&exportUnsignedVolume);
    _options->exportUnsignedVolume = _args->getBoolValue(&exportUnsignedVolume);

    Argument exportFloatVolume(
                "--export-float-volume",
                ARGUMENT_TYPE::BOOL,
                "Export an Ultraliser-specific float volume, where each voxel is stored in float. "
                "The data and header are stored in a single file with the extention .vol.");
    _args->addArgument(&exportFloatVolume);
    _options->exportFloatVolume = _args->getBoolValue(&exportFloatVolume);

    Argument exportRawVolume(
                "--export-raw-volume",
                ARGUMENT_TYPE::BOOL,
                "Export a raw volume, where each voxel is stored in 1, 2, 3 or 4 bytes depending "
                "on the volume type."
                "The resulting files are: .img file (contains data) and .hdr file (meta-data)");
    _args->addArgument(&exportRawVolume);
    _options->exportRawVolume = _args->getBoolValue(&exportRawVolume);

    Argument exportNRRDVolume(
                "--export-nrrd-volume",
                ARGUMENT_TYPE::BOOL,
                "Export a .nrrd volume that is compatible with VTK and can be loaded with "
                "Paraview for visualization purposes."
                "The resulting file contains the header and the data.");
    _args->addArgument(&exportNRRDVolume);
    _options->exportNRRDVolume = _args->getBoolValue(&exportNRRDVolume);

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
    Argument exportStackXY(
                "--export-stack-xy",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the Z-axis of the volume.");
    _args->addArgument(&exportStackXY);
    _options->exportStackXY = _args->getBoolValue(&exportStackXY);


    Argument exportStackXZ(
                "--export-stack-xz",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the Y-axis of the volume.");
    _args->addArgument(&exportStackXZ);
    _options->exportStackXZ = _args->getBoolValue(&exportStackXZ);

    Argument exportStackZY(
                "--export-stack-zy",
                ARGUMENT_TYPE::BOOL,
                "Generate an image stack along the X-axis of the volume.");
    _args->addArgument(&exportStackZY);
    _options->exportStackZY = _args->getBoolValue(&exportStackZY);
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

void AppArguments::addMeshJoiningArguments()
{
    Argument simpleMeshJoin(
                "--use-simple-mesh-join",
                ARGUMENT_TYPE::BOOL,
                "If this flag is set, the resulting mesh will be based on a mesh joint operation "
                "that simply merges all the given meshes into a single mesh object with multiple "
                "partition. "
                "Note that is this flag is set, all the mesh reconstruction parameters are ignored.");
    _args->addArgument(&simpleMeshJoin);
    _options->simpleMeshJoin = _args->getBoolValue(&simpleMeshJoin);
}

void AppArguments::addMeshExtractionArguments()
{
    Argument isosurfaceTechnique(
                "--isosurface-technique",
                ARGUMENT_TYPE::STRING,
                "Specify a technique to extract the isosurface from the volume: [mc, dmc]. "
                "By default, it is dmc (Dual Marching Cubes)",
                ARGUMENT_PRESENCE::OPTIONAL,
                "dmc");
    _args->addArgument(&isosurfaceTechnique);
    _options->isosurfaceTechnique = _args->getStringValue(&isosurfaceTechnique);
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

    Argument minDihedralAngle(
                "--min-dihedral-angle",
                ARGUMENT_TYPE::FLOAT,
                "The required minimum dihedral angle. Default 0.1",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.1");
    _args->addArgument(&minDihedralAngle);
    _options->minDihedralAngle = _args->getFloatValue(&minDihedralAngle);
}

void AppArguments::addLaplacianOperatorArguments()
{
    Argument laplacianIterations(
                "--laplacian-iterations",
                ARGUMENT_TYPE::INTEGER,
                "Number of iterations to smooth the reconstructed mesh with Laplacian filter. "
                "Default 10.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "10");
    _args->addArgument(&laplacianIterations);
    _options->laplacianIterations = _args->getIntegrValue(&laplacianIterations);
}

void AppArguments::addMorphologyExportArguments()
{
    Argument exportMorphology(
                "--export-morphology",
                ARGUMENT_TYPE::BOOL,
                "Export the resampled morphology skeleton used to create the mesh in (.swc) format.");
    _args->addArgument(&exportMorphology);
    _options->exportMorphology = _args->getBoolValue(&exportMorphology);
}

void AppArguments::addMeshExportArguments()
{
    Argument exportOBJ(
                "--export-obj-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export the resulting mesh(es) to Wavefront format (.obj).");
    _args->addArgument(&exportOBJ);
    _options->exportOBJ = _args->getBoolValue(&exportOBJ);

    Argument exportPLY(
                "--export-ply-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export the resulting mesh(es) to the Stanford triangle format (.ply).");
    _args->addArgument(&exportPLY);
    _options->exportPLY = _args->getBoolValue(&exportPLY);

    Argument exportOFF(
                "--export-off-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export the resulting mesh(es) to the object file format (.off).");
    _args->addArgument(&exportOFF);
    _options->exportOFF = _args->getBoolValue(&exportOFF);

    Argument exportSTL(
                "--export-stl-mesh",
                ARGUMENT_TYPE::BOOL,
                "Export the resulting mesh(es) to the  stereolithography CAD format (.stl).");
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
    _options->xScaleFactor = _args->getFloatValue(&xScaleFactor);

    Argument yScaleFactor(
                "--y-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Y-axis. Default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    _args->addArgument(&yScaleFactor);
    _options->yScaleFactor = _args->getFloatValue(&yScaleFactor);

    Argument zScaleFactor(
                "--z-scale",
                ARGUMENT_TYPE::FLOAT,
                "Scaling factor for the mesh along the Z-axis. Default 1.0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "1.0");
    _args->addArgument(&zScaleFactor);
    _options->zScaleFactor = _args->getFloatValue(&zScaleFactor);
}

void AppArguments::addMeshArguments()
{
    addMeshExtractionArguments();
    addMeshOptimizationArguments();
    addLaplacianOperatorArguments();
    addMeshScaleArguments();
    addMeshExportArguments();
}

void AppArguments::addProcessingArguments()
{
    Argument threads(
                "--threads",
                ARGUMENT_TYPE::INTEGER,
                "Number of cores used to process the parallel chunks in the code. "
                "If this value is set to 0, all the cores available in the system will be used."
                "The default valus is 0.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0");
    _args->addArgument(&threads);
    _options->threads = _args->getUnsignedIntegrValue(&threads);
}

void AppArguments::addAstrocyteSpecificArguments()
{
    Argument exportAstrocyteAtOrigin(
                "--export-at-origin",
                ARGUMENT_TYPE::BOOL,
                "Export the astrocyte mesh at the origin.");
    _args->addArgument(&exportAstrocyteAtOrigin);
    _options->exportAstrocyteAtOrigin = _args->getBoolValue(&exportAstrocyteAtOrigin);
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

void AppArguments::addPackingAlgorithmArguments()
{
    Argument packingAlgorithm(
                "--packing-algorithm",
                ARGUMENT_TYPE::STRING,
                "The packing algorithm used to create the proxy mesh. "
                "Options: [polylines, polylines-with-spheres, sdf]. Default [polylines].",
                ARGUMENT_PRESENCE::OPTIONAL,
                POLYLINE_PACKING);
    _args->addArgument(&packingAlgorithm);
    _options->packingAlgorithm = _args->getStringValue(&packingAlgorithm);
}

void AppArguments::addSuppressionArguments()
{
    Argument ignoreMarchingCubesMesh(
                "--ignore-marching-cubes-mesh",
                ARGUMENT_TYPE::BOOL,
                "If this flag is set, the mesh reconstructed with the marching cubes algorithm "
                "will not be written to disk.");
    _args->addArgument(&ignoreMarchingCubesMesh);
    _options->ignoreMarchingCubesMesh = _args->getBoolValue(&ignoreMarchingCubesMesh);

    Argument ignoreLaplacianMesh(
                "--ignore-laplacian-mesh",
                ARGUMENT_TYPE::BOOL,
                "If this flag is set, the mesh resulting from the application of the Laplacian "
                "operator will be ignored and will not be written to disk.");
    _args->addArgument(&ignoreLaplacianMesh);
    _options->ignoreLaplacianMesh = _args->getBoolValue(&ignoreLaplacianMesh);

    Argument ignoreOptimizedMesh(
                "--ignore-optimized-mesh",
                ARGUMENT_TYPE::BOOL,
                "If this flag is set, the optimized mesh will not be written to disk.");
    _args->addArgument(&ignoreOptimizedMesh);
    _options->ignoreOptimizedMesh = _args->getBoolValue(&ignoreOptimizedMesh);

    Argument ignoreWatertightMesh(
                "--ignore-optimized-mesh",
                ARGUMENT_TYPE::BOOL,
                "If this flag is set, the watertight mesh will not be written to disk.");
    _args->addArgument(&ignoreWatertightMesh);
    _options->ignoreWatertightMesh = _args->getBoolValue(&ignoreWatertightMesh);
}

void AppArguments::addNeuronMorphologyBranchOrderArguments()
{
    Argument axonBranchOrder(
                "--axon-branch-order",
                ARGUMENT_TYPE::INTEGER,
                "The maximum branch order applied in the neuron morphology axon mesh reconstruction.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "-1");
    _args->addArgument(&axonBranchOrder);
    int branchOrder = _args->getIntegrValue(&axonBranchOrder);
    if (branchOrder < 0)
        branchOrder = INT_MAX;
    _options->axonBranchOrder = branchOrder;

    Argument basalBranchOrder(
                "--basal-branch-order",
                ARGUMENT_TYPE::INTEGER,
                "The maximum branch order applied in the neuron morphology basal dendrites mesh reconstruction.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "-1");
    _args->addArgument(&basalBranchOrder);
    branchOrder = _args->getIntegrValue(&basalBranchOrder);
    if (branchOrder < 0)
        branchOrder = INT_MAX;
    _options->basalBranchOrder = branchOrder;
    
    Argument apicalBranchOrder(
                "--apical-branch-order",
                ARGUMENT_TYPE::INTEGER,
                "The maximum branch order applied in the neuron morphology apical dendrites mesh reconstruction.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "-1");
    _args->addArgument(&apicalBranchOrder);
    branchOrder = _args->getIntegrValue(&apicalBranchOrder);
    if (branchOrder < 0)
        branchOrder = INT_MAX;
    _options->apicalBranchOrder = branchOrder;
}

void AppArguments::addMorphologyAdjustmentParameters()
{
    Argument minSampleRadius(
                "--min-sample-radius",
                ARGUMENT_TYPE::FLOAT,
                "A user-defined value used to set the radius of any sample whose actual value is "
                "below this given threshold. "
                "The default value is 0.05 micron.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
    _args->addArgument(&minSampleRadius);
    _options->minSampleRadius = _args->getFloatValue(&minSampleRadius);

    Argument minROIRadius(
                "--roi-radius",
                ARGUMENT_TYPE::FLOAT,
                "The minimum radius of the regions of interest ROIs to select the  "
                "The default value is 0.05 micron.",
                ARGUMENT_PRESENCE::OPTIONAL,
                "0.05");
    _args->addArgument(&minSampleRadius);
    _options->minSampleRadius = _args->getFloatValue(&minSampleRadius);


    Argument autoParameters(
                "--auto-parameters",
                ARGUMENT_TYPE::BOOL,
                "When the auto-parameters flag is set, the morphology per-processing, voxelization, "
                "meshing and optimization parmeters will be automatically set to produce the most "
                "optimized mesh.");
    _args->addArgument(&autoParameters);
    _options->autoParameters = _args->getBoolValue(&autoParameters);
}

AppOptions* AppArguments::getOptions()
{
    // Parse the arguments
    _args->parse();

    return _options;
}

}

