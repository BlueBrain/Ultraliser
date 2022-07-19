/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include <inttypes.h>
#include <stdint.h>
#include <sys/stat.h>
#include <microscopes/MicroscopyConfiguration.h>
#include <microscopes/Microscopes.h>
#include <utilities/File.h>
#include <common/Common.h>
#include <utilities/Timer.h>

namespace Ultraliser
{
namespace PBRT
{

void createPBRTMicroscopyConfiguration(const MicroscopyConfiguration *config,
                                       Neurons neurons,
                                       Vector3f pMax,
                                       Vector3f pMin,
                                       Vector3f pMaxCenter,
                                       Vector3f pMinCenter,
                                       size_t volumeWidth,
                                       size_t volumeHeight,
                                       size_t volumeDepth,
                                       MICROSCOPE microscopeType)
{
    // PBRT parameters
    PBRTParameters paramsXY, paramsZY;
    paramsXY.pMin = pMinCenter;
    paramsZY.pMin = pMinCenter;
    paramsXY.pMax = pMaxCenter;
    paramsZY.pMax = pMaxCenter;
    paramsXY.resolution = config->imageResolution;
    paramsZY.resolution = config->imageResolution;
    paramsXY.lensRadius = config->lensRadius;
    paramsZY.lensRadius = config->lensRadius;

    // Compute the step
    Vector3f bounds = pMaxCenter - pMinCenter;
    paramsXY.step = config->stepScale * bounds.z() / (1.f * volumeDepth);
    paramsZY.step = config->stepScale * bounds.x() / (1.f * volumeWidth);

    // Volume prefix
    paramsXY.volumePrefix = config->prefix;
    paramsZY.volumePrefix = config->prefix;
    if (config->outputDirectory != EMPTY)
    {
        paramsXY.volumePrefix = config->outputDirectory + "/" + config->prefix;
        paramsZY.volumePrefix = config->outputDirectory + "/" + config->prefix;
    }

    std::string pbrtDirectory;
    if (microscopeType == MICROSCOPE::BFM)
        pbrtDirectory = BRIGHTFIELD_DIRECTORY;
    else if (microscopeType == MICROSCOPE::FLUORESCENCE)
        pbrtDirectory = FLUORESCENCE_DIRECTORY;
    else if (microscopeType == MICROSCOPE::LSFM)
        pbrtDirectory = LSFM_DIRECTORY;
    else
        LOG_ERROR("Wrong microscopy configuration");

    // Create the pbrt directory
    if (config->outputDirectory != EMPTY)
        pbrtDirectory = config->outputDirectory + "/" + pbrtDirectory;
    mkdir(pbrtDirectory.c_str(), S_IRWXU);

    // Create configuration file per neuron in the targets.
    if (config->focusAllObjects)
    {
        // Create the BFM fixed focus directory
        std::string directoryFixedFocus = pbrtDirectory + "/" + FOCUS_ALL_OBJECTS;
        mkdir(directoryFixedFocus.c_str(), S_IRWXU);

        // Create the projection directories
        std::string directoryXY = directoryFixedFocus + "/" + XY_DIRECTORY;
        std::string directoryZY = directoryFixedFocus + "/" + ZY_DIRECTORY;
        mkdir(directoryXY.c_str(), S_IRWXU);
        mkdir(directoryZY.c_str(), S_IRWXU);

        // Write the configuration per neuron
        TIMER_SET;
        LOOP_STARTS("All Objects Configurations");
        for (size_t iNeuron = 0; iNeuron < neurons.size(); iNeuron++)
        {
            LOOP_PROGRESS(iNeuron, neurons.size());

            // Create and load the mesh from the file
            std::stringstream file;
            file << neurons[iNeuron].gid << PBRT_EXTENSION;

            // Compute the focal distances, and update the configuration
            if (config->globalCoordinates)
            {
                paramsXY.focalDistance =
                        pMax.z() - neurons[iNeuron].somaPosition.z();
                paramsZY.focalDistance =
                        pMax.x() - neurons[iNeuron].somaPosition.x();
            }
            else
            {
                // The cameras are always located at pMax.z() and pMax.x(),
                // and when there is not transformation, then the neurons (or
                // their somata) are located at the origin (before translating)
                // their bounding box. Therefore, the focal distance should be
                // basically computed as the distance between pMax and the origin.
                paramsXY.focalDistance = pMax.z();
                paramsZY.focalDistance = pMax.x();
            }

            // Update the film prefix
            paramsXY.filmPrefix = directoryXY + "/" + neurons[iNeuron].gid;
            paramsZY.filmPrefix = directoryZY + "/" + neurons[iNeuron].gid;

            // Create the configuration files
            std::string configXY, configZY;
            if (microscopeType == MICROSCOPE::BFM)
            {
                configXY = PBRT::createConfigBrightFieldMicroscopeXY
                        (config->pbrtBFMConfig, paramsXY);
                configZY = PBRT::createConfigBrightFieldMicroscopeZY
                        (config->pbrtBFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::FLUORESCENCE)
            {
                configXY = PBRT::createConfigFluorescenceMicroscopeXY
                        (config->pbrtFMConfig, paramsXY);
                configZY = PBRT::createConfigFluorescenceMicroscopeZY
                        (config->pbrtFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::LSFM)
            {
                configXY = PBRT::createConfigLSFMXY
                        (config->pbrtLSFMConfig, paramsXY);
                configZY = PBRT::createConfigLSFMZY
                        (config->pbrtLSFMConfig, paramsZY);
            }
            else
                LOG_ERROR("Wrong microscopy configuration");

            // Writing the configuration files
            std::string pbrtConfigFile;

            // XY file
            pbrtConfigFile = directoryXY + "/" + file.str();
            std::fstream streamConfigXY;
            streamConfigXY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigXY << configXY;
            streamConfigXY.close();

            // ZY file
            pbrtConfigFile = directoryZY + "/" + file.str();
            std::fstream streamConfigZY;
            streamConfigZY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigZY << configZY;
            streamConfigZY.close();
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Create configuration file for each slice (step) in the volume
    if (config->focusAllSlices)
    {
        // Create the BFM varying focus directory
        std::string directoryVaryingFocus = pbrtDirectory + "/" + FOCUS_ALL_SLICES;
        mkdir(directoryVaryingFocus.c_str(), S_IRWXU);

        // Create the projection directories
        std::string directoryXY = directoryVaryingFocus + "/" + XY_DIRECTORY;
        std::string directoryZY = directoryVaryingFocus + "/" + ZY_DIRECTORY;
        mkdir(directoryXY.c_str(), S_IRWXU);
        mkdir(directoryZY.c_str(), S_IRWXU);

        // Compute the stepping, and then create the files
        float deltaZ = pMaxCenter.z() - pMinCenter.z();
        float stepZ = deltaZ / volumeDepth;

        // Scane through the tissue
        TIMER_SET;
        LOOP_STARTS("All Slices Configurations");
        for (size_t i = 0; i < volumeDepth; i++)
        {
            LOOP_PROGRESS(i, volumeDepth);

            // Compute the focal distances and update the parameters
            float focalDistanceXY = ((i + 1) * stepZ) + 0.5f;
            paramsXY.focalDistance = focalDistanceXY;

            // Film prefix
            std::stringstream prefix;
            prefix << std::setw(5) << std::setfill('0') << i;
            paramsXY.filmPrefix = directoryXY + "/" + prefix.str();

            // Create the configuration files
            std::string configXY;
            if (microscopeType == MICROSCOPE::BFM)
            {
                configXY = PBRT::createConfigBrightFieldMicroscopeXY
                        (config->pbrtBFMConfig, paramsXY);
            }
            else if (microscopeType == MICROSCOPE::FLUORESCENCE)
            {
                configXY = PBRT::createConfigFluorescenceMicroscopeXY
                        (config->pbrtFMConfig, paramsXY);
            }
            else if (microscopeType == MICROSCOPE::LSFM)
            {
                configXY = PBRT::createConfigLSFMXY
                        (config->pbrtLSFMConfig, paramsXY);
            }
            else
                LOG_ERROR("Wrong microscopy configuration");

            // Writing the configuration files
            std::string pbrtConfigFile;
            pbrtConfigFile = directoryXY + "/" + prefix.str() + PBRT_EXTENSION;
            std::fstream streamConfigXY;
            streamConfigXY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigXY << configXY;
            streamConfigXY.close();
        }

        float deltaX = pMaxCenter.x() - pMinCenter.x();
        float stepX = deltaX / volumeWidth;

        // Scane through the tissue
        for (size_t i = 0; i < volumeWidth; i++)
        {
            // Compute the focal distances and update the parameters
            float focalDistanceZY = ((i + 1) * stepX) + 0.5f;
            paramsZY.focalDistance = focalDistanceZY;

            // Film prefix
            std::stringstream prefix;
            prefix << std::setw(5) << std::setfill('0') << i;
            paramsZY.filmPrefix = directoryZY + "/" + prefix.str();

            // Create the configuration files
            std::string configZY;
            if (microscopeType == MICROSCOPE::BFM)
            {
                configZY = PBRT::createConfigBrightFieldMicroscopeZY
                        (config->pbrtBFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::FLUORESCENCE)
            {
                configZY = PBRT::createConfigFluorescenceMicroscopeZY
                        (config->pbrtFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::LSFM)
            {
                configZY = PBRT::createConfigLSFMZY
                        (config->pbrtLSFMConfig, paramsZY);
            }
            else
                LOG_ERROR("Wrong microscopy configuration");

            // Writing the configuration files
            std::string pbrtConfigFile;
            pbrtConfigFile = directoryZY + "/" + prefix.str() + PBRT_EXTENSION;
            std::fstream streamConfigZY;
            streamConfigZY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigZY << configZY;
            streamConfigZY.close();
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Create configuration file for each step in the volume
    if (config->numberFocalPlanes > 0)
    {
        // Create the BFM varying focus directory
        std::stringstream directoryVaryingFocus;
        directoryVaryingFocus << pbrtDirectory << "/"
                              << config->numberFocalPlanes
                              << "-slices";
        mkdir(directoryVaryingFocus.str().c_str(), S_IRWXU);

        // Create the projection directories
        std::string directoryXY = directoryVaryingFocus.str() + "/" + XY_DIRECTORY;
        std::string directoryZY = directoryVaryingFocus.str() + "/" + ZY_DIRECTORY;
        mkdir(directoryXY.c_str(), S_IRWXU);
        mkdir(directoryZY.c_str(), S_IRWXU);

        // Compute the stepping, and then create the files
        float deltaZ = pMaxCenter.z() - pMinCenter.z();
        float stepZ = deltaZ / config->numberFocalPlanes;
        float deltaX = pMaxCenter.x() - pMinCenter.x();
        float stepX = deltaX / config->numberFocalPlanes;

        // Scane through the tissue
        TIMER_SET;
        LOOP_STARTS("N Slices Configurations");
        for (size_t i = 0; i < config->numberFocalPlanes; i++)
        {
            LOOP_PROGRESS(i, config->numberFocalPlanes);

            // Compute the focal distances and update the parameters
            float focalDistanceXY = (i + 1) * stepZ;
            float focalDistanceZY = (i + 1) * stepX;
            paramsXY.focalDistance = focalDistanceXY;
            paramsZY.focalDistance = focalDistanceZY;

            // Film prefix
            std::stringstream prefix;
            prefix << std::setw(5) << std::setfill('0') << i;
            paramsXY.filmPrefix = directoryXY + "/" + prefix.str();
            paramsZY.filmPrefix = directoryZY + "/" + prefix.str();

            // Create the configuration files
            std::string configXY, configZY;
            if (microscopeType == MICROSCOPE::BFM)
            {
                configXY = PBRT::createConfigBrightFieldMicroscopeXY
                        (config->pbrtBFMConfig, paramsXY);
                configZY = PBRT::createConfigBrightFieldMicroscopeZY
                        (config->pbrtBFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::FLUORESCENCE)
            {
                configXY = PBRT::createConfigFluorescenceMicroscopeXY
                        (config->pbrtFMConfig, paramsXY);
                configZY = PBRT::createConfigFluorescenceMicroscopeZY
                        (config->pbrtFMConfig, paramsZY);
            }
            else if (microscopeType == MICROSCOPE::LSFM)
            {
                configXY = PBRT::createConfigLSFMXY
                        (config->pbrtLSFMConfig, paramsXY);
                configZY = PBRT::createConfigLSFMZY
                        (config->pbrtLSFMConfig, paramsZY);
            }
            else
                LOG_ERROR("Wrong microscopy configuration");

            // Writing the configuration files
            std::string pbrtConfigFile;
            pbrtConfigFile = directoryXY + "/" + prefix.str() + PBRT_EXTENSION;
            std::fstream streamConfigXY;
            streamConfigXY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigXY << configXY;
            streamConfigXY.close();

            pbrtConfigFile = directoryZY + "/" + prefix.str() + PBRT_EXTENSION;
            std::fstream streamConfigZY;
            streamConfigZY.open(pbrtConfigFile.c_str(), std::ios::out);
            streamConfigZY << configZY;
            streamConfigZY.close();
        }

        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void createPBRTMicroscopyConfigurations(const std::string &pbrtBFMConfig,
                                        const std::string &pbrtFMConfig,
                                        const std::string &pbrtLSFMConfig,
                                        const MicroscopyConfiguration& config,
                                        const Neurons& neurons,
                                        size_t volumeWidth,
                                        size_t volumeHeight,
                                        size_t volumeDepth,
                                        const Vector3f& pMaxInput,
                                        const Vector3f& pMinInput,
                                        const Vector3f& pMaxCenter,
                                        const Vector3f& pMinCenter)
{
    LOG_TITLE("Microscopy Configuration");

    // Create BFM configuration files
    if (pbrtBFMConfig != EMPTY)
    {
        LOG_STATUS("Brightfield Microscopy");
        createPBRTMicroscopyConfiguration(&config, neurons,
                                          pMaxInput, pMinInput,
                                          pMaxCenter, pMinCenter,
                                          volumeWidth, volumeHeight, volumeDepth,
                                          MICROSCOPE::BFM);
    }

    // Create fluorescence microscopy configuration files
    if(pbrtFMConfig != EMPTY)
    {
        LOG_STATUS("Epi-widefield Fluorescenec Microscopy");
        createPBRTMicroscopyConfiguration(&config, neurons,
                                          pMaxInput, pMinInput,
                                          pMaxCenter, pMinCenter,
                                          volumeWidth, volumeHeight, volumeDepth,
                                          MICROSCOPE::FLUORESCENCE);
    }

    // Create LSFM configuration files
    if(pbrtLSFMConfig != EMPTY)
    {
        LOG_STATUS("Light Sheet Fluoresecnce Microscopy");
        createPBRTMicroscopyConfiguration(&config, neurons,
                                          pMaxInput, pMinInput,
                                          pMaxCenter, pMinCenter,
                                          volumeWidth, volumeHeight, volumeDepth,
                                          MICROSCOPE::LSFM);
    }
}

}
}

