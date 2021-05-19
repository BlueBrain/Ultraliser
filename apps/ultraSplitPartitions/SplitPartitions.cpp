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

Options* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This application splits an input mesh with multiple partitions into multiple "
              "meshes with single partitions.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addMeshExportArguments();
    args->addDataArguments();

    // Get all the options
    Options* options = args->getOptions();

    LOG_TITLE("Creating Context");

    /// Validate the arguments
    if (!File::exists(options->inputMesh))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMesh.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Directory::exists(options->outputDirectory))
    {
        LOG_ERROR("The directory [ %s ] does NOT exist!", options->outputDirectory.c_str());
    }

    // Exporting formats, at least one of them must be there
    if (!(options->exportOBJ || options->exportPLY || options->exportOFF || options->exportSTL))
    {
        LOG_ERROR("The user must specify at least one output format of the "
                  "mesh to export: [--export-obj, --export-ply, --export-off, --export-stl]");
    }

    // If no prefix is given, use the directory name
    if (options->prefix == NO_DEFAULT_VALUE)
    {
        options->prefix = Directory::getName(options->inputMaskDirectory);
    }

    // Initialize context
    initializeContext(options);

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc , argv);

    // Load the mesh
    std::unique_ptr<Ultraliser::AdvancedMesh> inputMesh =
            std::make_unique<Ultraliser::AdvancedMesh>(options->inputMesh);

    // Split the partitions
    std::vector < Ultraliser::AdvancedMesh* > partitions = inputMesh->splitPartitions();

    // Export the partitions
    for (uint64_t i = 0; i < partitions.size(); ++i)
    {
        // Export the repaired mesh
        std::stringstream filePrefix;
        filePrefix << options->meshPrefix << i <<  MANIFOLD_SUFFIX;
        partitions[i]->exportMesh(filePrefix.str(),
                                  options->exportOBJ, options->exportPLY,
                                  options->exportOFF, options->exportSTL);
    }
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
