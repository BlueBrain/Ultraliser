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

#include "AppOptions.h"
#include <common/Common.h>

namespace Ultraliser
{

void AppOptions::verifyProcessingArguments()
{
#ifdef ULTRALISER_USE_OPENMP
    if (threads == 0)
    {
         omp_set_num_threads(omp_get_max_threads());
         LOG_SUCCESS("Running the application with [%d] thread(s)!", omp_get_max_threads());
    }
    else
    {
        if (threads  > omp_get_max_threads())
        {
            omp_set_num_threads(omp_get_max_threads());
            LOG_SUCCESS("Running the application with [%d] thread(s)!", omp_get_max_threads());
        }
        else
        {
            omp_set_num_threads(threads);
            LOG_SUCCESS("Running the application with [%d] thread(s)!", threads);
        }
    }
#else
    LOG_WARNING("OpenMP NOT Found! Using 1 core to run the application!");
#endif
}

void AppOptions::verifyInputMeshArgument()
{
    if (!File::exists(inputMeshPath))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", inputMeshPath.c_str());
    }

    std::string path = inputMeshPath;
    String::toLower(path);
    if (!(String::subStringFound(path, ".obj") || String::subStringFound(path, ".ply") ||
          String::subStringFound(path, ".stl") || String::subStringFound(path, ".off") ||
          String::subStringFound(path, ".h5")))
    {
        LOG_ERROR("Unsupported mesh type [%s]!", inputMeshPath.c_str());
    }
}

void AppOptions::verifyInputMeshesDirectoryArgument()
{
    if (!Directory::exists(inputMeshesDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist! ", inputMeshesDirectory.c_str());
    }
}

void AppOptions::verifyInputMorphologyArgument()
{
    if (!File::exists(inputMorphologyPath))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", inputMorphologyPath.c_str());
    }
}

void AppOptions::verifyInputMaskDirectoryArgument()
{
    if (!Directory::exists(inputMaskDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist! ", inputMaskDirectory.c_str());
    }
}

void AppOptions::verifyOutputDirectoryArgument()
{
    // Try to make the output directory
    mkdir(outputDirectory.c_str(), 0777);
    if (!Directory::exists(outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", outputDirectory.c_str());
    }
}

void AppOptions::verifyBoudsFileArgument()
{
    if (boundsFile == NO_DEFAULT_VALUE)
    {
        LOG_WARNING("The bounding box of the volume will be computed on the fly");
        boundsFile = EMPTY;
    }
    else
    {
        LOG_WARNING("The bounding box of the volume will be loaded from [ %s ]", boundsFile.c_str());
    }
}

void AppOptions::verifyMeshPrefixArgument()
{
    // If no prefix is given, use the file name
    if (prefix == NO_DEFAULT_VALUE)
    {
        prefix = File::getName(inputMeshPath);
    }
}

void AppOptions::verifyMeshesPrefixArgument()
{
    // If no prefix is given, use the file name
    if (prefix == NO_DEFAULT_VALUE)
    {
        prefix = Directory::getName(inputMeshesDirectory);
    }
}

void AppOptions::verifyMorphologyPrefixArgument()
{
    // If no prefix is given, use the file name
    if (prefix == NO_DEFAULT_VALUE)
    {
        prefix = File::getName(inputMorphologyPath);
    }
}

void AppOptions::verifyMaskPrefixArgument()
{
    // If no prefix is given, use the file name
    if (prefix == NO_DEFAULT_VALUE)
    {
        prefix = Directory::getName(inputMaskDirectory);
    }
}

void AppOptions::verifyVolumePrefixArgument()
{
    // If no prefix is given, use the file name
    if (prefix == NO_DEFAULT_VALUE)
    {
        prefix = Ultraliser::File::getName(inputVolumePath);
    }
}

void AppOptions::verifyMaskDimensionsArguments()
{
    if (maskWidth == 0 || maskHeight == 0)
    {
        LOG_ERROR("Mask dimensions cannot be zero: [%d x %d]", maskWidth, maskHeight);
    }
}

void AppOptions::verifyMorphologyExtractionArguments()
{
    if (bboxWidth == 0 || bboxHeight == 0 || bboxDepth == 0)
    {
        LOG_ERROR("Bounding box dimensions cannot be zero: [%f x %f x %f]",
                  bboxWidth, bboxHeight, bboxDepth);
    }
}

void AppOptions::verifyPackingAlgorithmArgument()
{
    if (!(packingAlgorithm == POLYLINE_PACKING ||
          packingAlgorithm == POLYLINE_SPHERE_PACKING ||
          packingAlgorithm == SDF_PACKING))
    {
        LOG_ERROR("The given packing algorithm [%s] is not correct! "
                  "Please use one of the following packing algorithms: "
                  "[%s, %s or %s]",
                  packingAlgorithm.c_str(),
                  POLYLINE_PACKING.c_str(), POLYLINE_SPHERE_PACKING.c_str(), SDF_PACKING.c_str());
    }
}

void AppOptions::verifyMeshExportLogic()
{

}

void AppOptions::verifyIsoSurfaceExtractionArgument()
{
    if (!(isosurfaceTechnique == DMC_STRING || isosurfaceTechnique == MC_STRING))
    {
        LOG_ERROR("You must specify a correct isosurface extraction technique: [dmc or mc]");
    }
}

void AppOptions::verifyIsoOptionArgument()
{
    if (isoOption == "isovalue")
    {
        LOG_SUCCESS("Isovalue [%zu] will be used to segment the volume", isoValue);
    }
    else if (isoOption == ISOVALUES_STRING)
    {
        if (isovaluesFile == NO_DEFAULT_VALUE)
        {
            LOG_ERROR("The isovalues-file is NOT given. Please specify a text file containing all "
                      "the iso values required to segment the volume.");
        }

        // Check if the file exists
        if (!File::exists(isovaluesFile))
        {
            LOG_ERROR("The file [ %s ] does NOT exist!", isovaluesFile.c_str());
        }
    }
    else if (isoOption == MIN_ISOVALUE_STRING)
    {
        LOG_SUCCESS("Isovalues [%zu-inf] will be used to segment the volume", minIsoValue);
    }
    else if (isoOption == MAX_ISOVALUE_STRING)
    {
        LOG_SUCCESS("Isovalues [0-%zu] will be used to segment the volume", maxIsoValue);
    }
    else if (isoOption == ISOVALUE_RANGE_STRING)
    {
        if (minIsoValue > maxIsoValue)
        {
            LOG_ERROR("The minimum isovalue CANNOT be greater than the maximum isovalue");
        }
        else
        {
            LOG_SUCCESS("Isovalues [%zu-%zu] range will be used to segment the volume",
                        minIsoValue, maxIsoValue);
        }
    }
    else if (isoOption == NON_ZERO_STRING)
    {
        LOG_WARNING("All the NON-zero voxels in the volume will be used to segment the volume");
    }
    else
    {
        LOG_ERROR("The iso-option [%s] CANNOT be recognized, see the options from the help");
    }
}

void AppOptions::verifyMeshExportArguments()
{
    // Exporting formats, at least one of them must be there
    if (!(exportOBJ || exportPLY || exportOFF || exportSTL))
    {
        LOG_ERROR("You must specify at least one valid mesh format to export:"
                  "[--export-obj-mesh, --export-ply-mesh, --export-off-mesh, --export-stl-mesh]");
    }

    // Then ensure that there is at least one mesh will be exported.
    verifyMeshExportLogic();
}

void AppOptions::verifyVascularMorphologyExportArguments()
{
    // Exporting formats, at least one of them must be there
    if (exportVMV)
    {
        LOG_ERROR("You must specify at least one valid morphology format to export:"
                  "[--export-vmv-morphology]");
    }
}

void AppOptions::verifyNeuronalMorphologyExportArguments()
{
    // Exporting formats, at least one of them must be there
    if (!exportBranches)
    {
        LOG_ERROR("You must specify at least one valid morphology format to export:"
                  "[--export-morphology-branches]");
    }
}

void AppOptions::verifyVolumeExportArguments()
{
    if (!(exportBitVolume || exportRawVolume || exportNRRDVolume))
    {
        LOG_ERROR("You must specify at least one output format of the volume: "
                  "[--export-bit-volume, --export-byte-volume, --export-nrrd-volume]");
    }
}

void AppOptions::createRespectiveDirectories()
{
    // Meshes directory
    if (exportOBJ || exportPLY || exportOFF || exportSTL)
    {
        std::stringstream path;
        path << outputDirectory << "/" << MESHES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Volumes directory
    if (exportBitVolume || exportRawVolume || exportNRRDVolume)
    {
        std::stringstream path;
        path << outputDirectory << "/" << VOLUMES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Projections directory
    if (projectXY || projectXZ || projectZY)
    {
        std::stringstream path;
        path << outputDirectory << "/" << PROJECTIONS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Stacks directory
    if (exportStackXY || exportStackXZ || exportStackZY)
    {
        std::stringstream path;
        path << outputDirectory << "/" << STACKS_SIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Morphology directory
    if (exportVMV || exportBranches)
    {
        std::stringstream path;
        path << outputDirectory << "/" << MORPHOLOGIES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Statistics directory
    if (writeStatistics)
    {
        std::stringstream path;
        path << outputDirectory << "/" << STATISTICS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Distributions directory
    if (writeDistributions)
    {
        std::stringstream path;
        path << outputDirectory << "/" << DISTRIBUTIONS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }
}

void AppOptions::initializeContext()
{
    // Construct the prefixes once and for all
    outputPrefix =
            outputDirectory + "/" + prefix;
    morphologyPrefix =
            outputDirectory + "/" + MORPHOLOGIES_DIRECTORY +  "/" + prefix;
    meshPrefix =
            outputDirectory + "/" + MESHES_DIRECTORY +  "/" + prefix;
    volumePrefix =
            outputDirectory + "/" + VOLUMES_DIRECTORY +  "/" + prefix;
    projectionPrefix =
            outputDirectory + "/" + PROJECTIONS_DIRECTORY +  "/" + prefix;
    statisticsPrefix =
            outputDirectory + "/" + STATISTICS_DIRECTORY +  "/" + prefix;
    distributionsPrefix =
            outputDirectory + "/" + DISTRIBUTIONS_DIRECTORY +  "/" + prefix;

    // Create the respective directories
    createRespectiveDirectories();

    LOG_TITLE("Ultralizing");
    LOG_SUCCESS("Output Directory [ %s ]", outputDirectory.c_str());
}

}
