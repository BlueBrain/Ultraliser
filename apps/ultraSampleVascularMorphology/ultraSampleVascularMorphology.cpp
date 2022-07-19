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
    std::unique_ptr< AppArguments > args = std::make_unique <AppArguments>(argc, argv,
              "This tools verifies the connecitivty of the VMV datasets and creates the "
              "connectivity if it is missing.");

    args->addInputMorphologyArguments();
    args->addOutputArguments();
    args->addMorphologyExtractionArguments();

    // Get all the options
    AppOptions* options = args->getOptions();

    LOG_TITLE("Creating Context");

    // Verify the arguments after parsing them and extracting the application options.
    options->verifyInputMorphologyArgument();
    options->verifyMorphologyExtractionArguments();
    options->verifyOutputDirectoryArgument();
    options->verifyMorphologyPrefixArgument();

    // Initialize context
    options->initializeContext();

    // Create the morphologies directory directory
    std::stringstream path;
    path << options->outputDirectory << "/" << MORPHOLOGIES_DIRECTORY;
    mkdir(path.str().c_str(), 0777);

    // Return the executable options
    return options;
}

void run(int argc , const char** argv)
{
    // Parse the arguments and get the values
    auto options = parseArguments(argc, argv);

    // Read the file into a morphology structure
    auto fullMorphology = readVascularMorphology(options->inputMorphologyPath);

    // Compute the morphology bounding box
    Vector3f pMin, pMax, bounds, center;
    fullMorphology->getBoundingBox(pMin, pMax, bounds, center);

    // Extract a subset from the big morphology
    auto extractedMorphology = fullMorphology->extractRegion(center,
                                                             options->bboxWidth,
                                                             options->bboxHeight,
                                                             options->bboxDepth);

    // Export the vascular morphology to VMV file
    extractedMorphology->exportVascularMorphologyVMV(options->morphologyPrefix);
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
