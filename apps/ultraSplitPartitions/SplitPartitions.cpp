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

namespace Ultraliser
{

Options* parseArguments(const int& argc , const char** argv)
{
    // Arguments
    std::unique_ptr< Args > args = std::make_unique <Args>(argc, argv,
              "This application splits an input mesh with multiple partitions into multiple "
              "meshes with single partitions.");

    Ultraliser::Argument inputMesh(
                "--mesh",
                ARGUMENT_TYPE::STRING,
                "The full path to the mesh.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&inputMesh);

    Ultraliser::Argument outputDirectory(
                "--output-directory",
                ARGUMENT_TYPE::STRING,
                "Output directory where the volume will be generated.",
                ARGUMENT_PRESENCE::MANDATORY);
    args->addArgument(&outputDirectory);

    Ultraliser::Argument prefix(
                "--prefix",
                ARGUMENT_TYPE::STRING,
                "Just a prefix that will be used to label the output files. "
                "If this is not given by the user, the name of the mesh file "
                "will be used.");
    args->addArgument(&prefix);

    Ultraliser::Argument exportOBJ(
                "--export-obj",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OBJ file.");
    args->addArgument(&exportOBJ);

    Ultraliser::Argument exportPLY(
                "--export-ply",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .PLY file.");
    args->addArgument(&exportPLY);

    Ultraliser::Argument exportOFF(
                "--export-off",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .OFF file.");
    args->addArgument(&exportOFF);

    Ultraliser::Argument exportSTL(
                "--export-stl",
                ARGUMENT_TYPE::BOOL,
                "Export the output mesh to an .STL file.");
    args->addArgument(&exportSTL);

    Argument writeStatistics(
                "--stats",
                ARGUMENT_TYPE::BOOL,
                "Write the statistics.");
    args->addArgument(&writeStatistics);

    Argument writeDistributions(
                "--dists",
                ARGUMENT_TYPE::BOOL,
                "Write the distributions.");
    args->addArgument(&writeDistributions);

    // Parse the command line options
    args->parse();

    // Construct the options
    Options* options = new Options();

    // Get all the options
    options->inputMesh = args->getStringValue(&inputMesh);
    options->outputDirectory = args->getStringValue(&outputDirectory);
    options->prefix = args->getStringValue(&prefix);

    // Mesh exports, file formats
    options->exportOBJ = args->getBoolValue(&exportOBJ);
    options->exportPLY = args->getBoolValue(&exportPLY);
    options->exportOFF = args->getBoolValue(&exportOFF);
    options->exportSTL = args->getBoolValue(&exportSTL);

    // Statistics and distributions
    options->writeStatistics = args->getBoolValue(&writeStatistics);
    options->writeDistributions = args->getBoolValue(&writeDistributions);

    /// Validate the arguments
    if (!Ultraliser::File::exists(options->inputMesh))
    {
        LOG_ERROR("The file [ %s ] does NOT exist! ", options->inputMesh.c_str());
    }

    // Try to make the output directory
    mkdir(options->outputDirectory.c_str(), 0777);
    if (!Ultraliser::Directory::exists(options->outputDirectory))
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
        options->prefix = Ultraliser::Directory::getName(options->inputMaskDirectory);
    }

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
