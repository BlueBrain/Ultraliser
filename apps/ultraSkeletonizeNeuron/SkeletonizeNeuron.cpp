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
#include <data/meshes/simple/IcoSphere.h>
#include <algorithms/skeletonization/NeuronSkeletonizer.h>
#include <geometry/Sphere.h>

namespace Ultraliser
{

AppOptions* parseArguments(const int& argc , const char** argv)
{
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tool reconstructs a watertight polygonal mesh from an input "
              "non-watertight mesh. The generated mesh can be also optimized to "
              "reduce the number of triangles while preserving the volume. "
              "The output mesh is guaranteed in all cases to be two-manifold");

    args->addInputMeshArguments();
    args->addOutputArguments();
    args->addNeuronalMorphologyExportArguments();
    args->addVolumeArguments();
    args->addMeshArguments();
    args->addSuppressionArguments();
    args->addDataArguments();
    args->addProcessingArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMeshArgument();
    options->verifyOutputDirectoryArgument();
    options->verifyBoudsFileArgument();
    options->verifyMeshPrefixArgument();
    options->verifyIsoSurfaceExtractionArgument();
    options->verifyNeuronalMorphologyExportArguments();
    options->verifyProcessingArguments();

    // Initialize context, once everything is in place and all the options are verified
    options->initializeContext();

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the tool options
    auto options = parseArguments(argc, argv);

    // Load the input mesh
    auto inputMesh = loadInputMesh(options);
    // inputMesh->smoothSurface(10);

    auto prefix = options->projectionPrefix;

    // Create the volume from the mesh
    auto solidVolume = createVolumeGrid(inputMesh, options);

    // Surface voxelization
    solidVolume->surfaceVoxelization(inputMesh, false, false, 1.0);
    solidVolume->solidVoxelization(options->voxelizationAxis);
    solidVolume->surfaceVoxelization(inputMesh, false, false, 0.5);

    auto solidRanges = solidVolume->getSolidRanges();

//    for (size_t i = 0; i < solidRanges.size(); ++i)
//    {
//        for (size_t j = 0; j < solidRanges[i].size(); ++j)
//        {
//            if (solidRanges[i][j].size() > 0)
//            std::cout << solidRanges[i][j].size() << ": ";

//            for (size_t k = 0; k < solidRanges[i][j].size(); ++k)
//            {
//                std::cout << "[" << solidRanges[i][j][k]->lower
//                          << ", " << solidRanges[i][j][k]->upper
//                          << "], ";
//            }
//            if (solidRanges[i][j].size() > 0)
//                std::cout << "\n";
//        }
//    }

    size_t xx1 = 0, xx2 = 0;

    std::stringstream stream;
    for (size_t i =0; i < solidVolume->getWidth(); ++i)
    {
        for (size_t j = 0; j < solidVolume->getHeight(); ++j)
        {

            for (size_t k =0; k < solidVolume->getDepth(); ++k)
            {
                if (solidVolume->isFilled(i, j, k))
                {
                    xx1++;
                    stream << i << ", " << j << ", " << k << "\n";
                }
            }
        }
    }

    std::fstream fileStream;
    std::string filePath = "/data/microns-skeletonization-meshes/skeletonization-output-spines/morphologies/normal.txt";
    fileStream.open(filePath, std::ios::out);
    fileStream << stream.str();
    fileStream.close();

    std::stringstream stream2;
    for (size_t i = 0; i < solidRanges.size(); ++i)
    {
        for (size_t j = 0; j < solidRanges[i].size(); ++j)
        {
            for (size_t k = 0; k < solidRanges[i][j].size(); ++k)
            {
                // stream2 << i << ", " << j << "\n";
                for (size_t w = solidRanges[i][j][k]->lower; w <= solidRanges[i][j][k]->upper; w++)
                {
                    xx2++;
                   stream2 << i << ", " << j << ", " << w << "\n";
                }
            }
        }
    }

    std::cout << xx1 << ", " << xx2 << "\n";
    filePath = "/data/microns-skeletonization-meshes/skeletonization-output-spines/morphologies/solid.txt";
    fileStream.open(filePath, std::ios::out);
    fileStream << stream2.str();
    fileStream.close();
    exit(0);


    // Create a skeletonization object
    NeuronSkeletonizer* skeletonizer = new NeuronSkeletonizer(solidVolume, inputMesh);

    skeletonizer->skeletonizeVolume();

    // Project the volume
    solidVolume->project(options->morphologyPrefix, true, true, true);


    // Build the graph
    skeletonizer->constructGraph();

    skeletonizer->exportIndividualBranches(options->morphologyPrefix + "-1");


    skeletonizer->segmentComponents();

    skeletonizer->exportIndividualBranches(options->morphologyPrefix + "-2");


    skeletonizer->exportSWCFile(options->morphologyPrefix);

    skeletonizer->exportSomaMesh(options->meshPrefix,
                                 options->exportOBJ, options->exportPLY,
                                 options->exportOFF, options->exportSTL);
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
