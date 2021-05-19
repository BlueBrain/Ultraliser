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

#include "RunTime.h"
#include "Options.hh"
#include <common/Common.h>

namespace Ultraliser
{

void createRespectiveDirectories(const Options* options)
{
    // Meshes directory
    if (options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << MESHES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Volumes directory
    if (options->writeBitVolume || options->writeByteVolume || options->writeNRRDVolume)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << VOLUMES_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Projections directory
    if (options->projectXY || options->projectXZ || options->projectZY)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << PROJECTIONS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Stacks directory
    if (options->stackXY || options->stackXZ || options->stackZY)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << STACKS_SIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Statistics directory
    if (options->writeStatistics)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << STATISTICS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }

    // Distributions directory
    if (options->writeDistributions)
    {
        std::stringstream path;
        path << options->outputDirectory << "/" << DISTRIBUTIONS_DIRECTORY;
        mkdir(path.str().c_str(), 0777);
    }
}

void initializeContext(Options* options)
{
    // Construct the prefixes once and for all
    options->outputPrefix =
            options->outputDirectory + "/" + options->prefix;
    options->meshPrefix =
            options->outputDirectory + "/" + MESHES_DIRECTORY +  "/" + options->prefix;
    options->volumePrefix =
            options->outputDirectory + "/" + VOLUMES_DIRECTORY +  "/" + options->prefix;
    options->projectionPrefix =
            options->outputDirectory + "/" + PROJECTIONS_DIRECTORY +  "/" + options->prefix;
    options->statisticsPrefix =
            options->outputDirectory + "/" + STATISTICS_DIRECTORY +  "/" + options->prefix;
    options->distributionsPrefix =
            options->outputDirectory + "/" + DISTRIBUTIONS_DIRECTORY +  "/" + options->prefix;

    // Create the respective directories
    createRespectiveDirectories(options);

    LOG_TITLE("Ultralizing");
    LOG_SUCCESS("Output Directory [ %s ]", options->outputDirectory.c_str());
}

}
