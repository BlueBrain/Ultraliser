/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

namespace Ultraliser
{

void AppOptions::verifyInputMeshArgument()
{
    if (!File::exists(inputMeshPath))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", inputMeshPath.c_str());
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

void AppOptions::verifyMaskDimensionsArguments()
{
    if (maskWidth == 0 || maskHeight == 0)
    {
        LOG_ERROR("Mask dimensions cannot be zero: [%d x %d]", maskWidth, maskHeight);
    }
}

void AppOptions::verifyMeshExportLogic()
{

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

void AppOptions::verifyVolumeExportArguments()
{
    if (!(exportBitVolume || exportByteVolume || exportNRRDVolume))
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
    if (exportBitVolume || exportByteVolume || exportNRRDVolume)
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
