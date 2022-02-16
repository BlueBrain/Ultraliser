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

#ifndef ULTRALISER_PBRT_PBRT_PARAMETERS_H
#define ULTRALISER_PBRT_PBRT_PARAMETERS_H

#include <string>
#include <math/Math.h>

namespace Ultraliser
{
namespace PBRT
{

/**
 * @brief The PBRTParameters struct
 */
struct PBRTParameters
{
    /**
     * @brief pMin
     */
    Vector3f pMin;

    /**
     * @brief pMax
     */
    Vector3f pMax;

    /**
     * @brief focalDistance
     */
    float focalDistance;

    /**
     * @brief resolution
     */
    uint64_t resolution;

    /**
     * @brief filmPrefix
     */
    std::string filmPrefix;

    /**
     * @brief volumePrefix
     */
    std::string volumePrefix;

    /**
     * @brief step
     */
    float step;

    /**
     * @brief lensRadius
     */
    float lensRadius;
};

}
}

#endif // ULTRALISER_PBRT_PBRT_PARAMETERS_H

