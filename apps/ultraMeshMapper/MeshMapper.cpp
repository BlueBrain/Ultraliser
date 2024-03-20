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
    std::unique_ptr< AppArguments > args = std::make_unique < AppArguments >(
        argc, argv, COPYRIGHT
        "This application maps a source mesh to a destination mesh");

    std::cout << "22\n";
    args->addMeshMappingArguments();
    args->addOutputArguments();
    args->addMeshExportArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyMeshMappingArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyMeshPrefixArgument();
    options->verifyMeshExportArguments();

    // Initialize context
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc, argv);

    std::cout << "11\n";

    // Load the source and destination meshes
    auto sourceMesh = new Mesh(options->sourceMesh);
    auto destinationMesh = new Mesh(options->destinationMesh);

    std::vector< Vector3f > destinationMeshCloud;
    auto destinationMeshVertices = destinationMesh->_vertices;
    auto destinationMeshNumberVertices = destinationMesh->getNumberVertices();
    for (size_t i = 0; i < destinationMeshNumberVertices; ++i)
    {
        auto point = destinationMeshVertices[i];
        destinationMeshCloud.push_back(Vector3f(point.x(), point.y(), point.z()));
    }

    // Map the source mesh to the destination mesh
    sourceMesh->kdTreeMapping(destinationMeshCloud);

    // Export the source mesh
    sourceMesh->exportMesh(options->meshPrefix + MAPPED_SUFFIX,
                           options->exportOBJ,
                           options->exportPLY,
                           options->exportOFF,
                           options->exportSTL);
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
