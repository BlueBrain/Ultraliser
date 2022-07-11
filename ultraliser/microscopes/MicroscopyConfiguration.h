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

#pragma once

#include <data/NeuronData.h>
#include <microscopes/PBRTParameters.h>
#include <microscopes/PBRT.h>

namespace Ultraliser
{
namespace PBRT
{

/**
 * @brief The MicroscopyConfiguration struct
 * Paramters that should be used to create the PBRT microscopy configurations.
 */
struct MicroscopyConfiguration
{
    /**
     * @brief prefix
     * Output files prefix.
     */
    std::string prefix;

    /**
     * @brief outputDirectory
     * The output diretory where the files will be written.
     */
    std::string outputDirectory;

    /**
     * @brief pbrtBFMConfig
     * Input BFM microscopy configuration template.
     */
    std::string pbrtBFMConfig;

    /**
     * @brief pbrtFMConfig
     * Input FM microscopy configuration template.
     */
    std::string pbrtFMConfig;

    /**
     * @brief pbrtLSFMConfig
     * Input LSFM microscopy configuration template.
     */
    std::string pbrtLSFMConfig;

    /**
     * @brief imageResolution
     * Image resolution of the microscope.
     */
    uint64_t imageResolution;

    /**
     * @brief numberFocalPlanes
     * Number of focal planes used to slice the volume.
     */
    uint64_t numberFocalPlanes;

    /**
     * @brief focusAllObjects
     * Create configurations that focus on all the objects in the volume.
     */
    bool focusAllObjects;

    /**
     * @brief focusAllSlices
     * Create configurations that focus on all the slices in the volume for a
     * one-to-one mapping with the volume.
     */
    bool focusAllSlices;

    /**
     * @brief lensRadius
     * The lens radius of the microscope.
     */
    float lensRadius;

    /**
     * @brief stepScale
     * A scale factor that is used to scale the stepping through the volume.
     */
    float stepScale;

    /**
     * @brief globalCoordinates
     * Set the neurons in the global coordinates
     */
    bool globalCoordinates;
};

/**
 * @brief createPBRTMicroscopyConfiguration
 * Creates the PBRT configuration files of the simulated microscope.
 *
 * @param config
 * Generic microscope configuration.
 * @param neurons
 * A list of neurons.
 * @param pMax
 * pMax of the original bounding box, before centering around the origin.
 * @param pMin
 * pMin of the original bounding box, before centering around the origin.
 * @param pMaxCenter
 * pMax of the centered bounding box.
 * @param pMinCenter
 * pMin of the centered bounding box.
 * @param volumeWidth
 * Volume width.
 * @param volumeHeight
 * Volume height.
 * @param volumeDepth
 * Volume depth.
 * @param microscopeType
 * The type of the microscope.
 */
void createPBRTMicroscopyConfiguration(const MicroscopyConfiguration* config,
                                       Neurons neurons,
                                       Vector3f pMax,
                                       Vector3f pMin,
                                       Vector3f pMaxCenter,
                                       Vector3f pMinCenter,
                                       uint64_t volumeWidth,
                                       uint64_t volumeHeight,
                                       uint64_t volumeDepth,
                                       MICROSCOPE microscopeType);

/**
 * @brief createPBRTMicroscopyConfigurations
 * Creates the PBRT configuration files of the simulated microscope.
 *
 * @param pbrtBFMConfig
 * Template configuration for the BFM.
 * @param pbrtFMConfig
 * Template configuration for the FM.
 * @param pbrtLSFMConfig
 * Template configuration for the LSFM.
 * @param config
 * Generic microscope configuration.
 * @param neurons
 * A list of neurons.
 * @param volumeWidth
 * Volume width.
 * @param volumeHeight
 * Volume height.
 * @param volumeDepth
 * Volume depth.
 * @param pMaxInput
 * pMax of the original bounding box, before centering around the origin.
 * @param pMinInput
 * pMin of the original bounding box, before centering around the origin.
 * @param pMaxCenter
 * pMax of the centered bounding box.
 * @param pMinCenter
 * pMin of the centered bounding box.
 */
void createPBRTMicroscopyConfigurations(const std::string &pbrtBFMConfig,
                                        const std::string &pbrtFMConfig,
                                        const std::string &pbrtLSFMConfig,
                                        const MicroscopyConfiguration& config,
                                        const Neurons& neurons,
                                        uint64_t volumeWidth,
                                        uint64_t volumeHeight,
                                        uint64_t volumeDepth,
                                        const Vector3f& pMaxInput,
                                        const Vector3f& pMinInput,
                                        const Vector3f& pMaxCenter,
                                        const Vector3f& pMinCenter);

}
}
