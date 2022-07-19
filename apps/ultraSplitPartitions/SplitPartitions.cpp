/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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
              "This application splits an input mesh with multiple partitions into multiple "
              "meshes with single partitions.");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addMeshExportArguments();
    args->addDataArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshExportArguments();
    options->verifyMeshPrefixArgument();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc , argv);

    // Load the mesh
    std::unique_ptr<Ultraliser::AdvancedMesh> inputMesh =
            std::make_unique<Ultraliser::AdvancedMesh>(options->inputMeshPath);

    // Split the partitions
    std::vector < Ultraliser::AdvancedMesh* > partitions = inputMesh->splitPartitions();

    // Export the partitions
    for (size_t i = 0; i < partitions.size(); ++i)
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
