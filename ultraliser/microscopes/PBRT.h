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

#include <string>
#include <sstream>
#include <math/Math.h>
#include <microscopes/PBRTParameters.h>
#include <data/NeuronData.h>

namespace Ultraliser
{
namespace PBRT
{

/**
 * @brief The PROJECTION enum
 */
enum PROJECTION
{
    XY,
    ZY
};

/**
 * @brief The MICROSCOPE enum
 */
enum MICROSCOPE
{
    BFM,
    FLUORESCENCE,
    LSFM
};

/**
 * @brief The Resolution struct
 */
struct Resolution
{
    float x;
    float y;
};

/**
 * @brief computeFilmResolution
 * Computes the PBRT film resolution given the scene data.
 *
 * @param projection
 * Microscope projection.
 * @param baseResolution
 * The base resolution of the film.
 * @param pMin
 * Bounding box pMin.
 * @param pMax
 * Bounding box pMax.
 * @return
 * Film resolution.
 */
Resolution computeFilmResolution(PROJECTION projection,
                                 size_t baseResolution,
                                 Vector3f pMin,
                                 Vector3f pMax);

/**
 * @brief createDepthOfFieldFile
 * Creates the PBRT configuration file at a specific depth of field.
 *
 * @param dofDirectory
 * The output directory where the files will be created.
 * @param neuron
 * A list of neurons in the scene.
 * @param pMin
 * Bounding box pMin.
 * @param pMax
 * Bounding box pMax.
 */
void createDepthOfFieldFile(const std::string dofDirectory,
                            const Neuron& neuron,
                            Vector3f pMax, Vector3f pMin);

}
}
