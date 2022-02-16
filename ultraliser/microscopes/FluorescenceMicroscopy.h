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

#ifndef PBRT_FLUORESCENCE_MICROSCOPY_H
#define PBRT_FLUORESCENCE_MICROSCOPY_H

#include <math/Math.h>
#include <microscopes/PBRTParameters.h>
#include <microscopes/PBRT.h>

namespace Ultraliser
{
namespace PBRT
{

/**
 * @brief createConfigFluorescenceMicroscopeXY
 * Create PBRT configuration files for the fluoresecnce microscope along the XY axis.
 *
 * @param inputConfig
 * Input template configuration.
 * @param params
 * System parameters.
 * @return
 * An std::string with the microscopy configuration to be written to a file directly.
 */
std::string createConfigFluorescenceMicroscopeXY(const std::string inputConfig,
                                                 const PBRTParameters params);

/**
 * @brief createConfigFluorescenceMicroscopeZY
 * Create PBRT configuration files for the fluoresecnce microscope along the XY axis.
 *
 * @param inputConfig
 * Input template configuration.
 * @param params
 * System parameters.
 * @return
 * An std::string with the microscopy configuration to be written to a file directly.
 */
std::string createConfigFluorescenceMicroscopeZY(const std::string inputConfig,
                                                 const PBRTParameters params);

}
}

#endif // PBRT_FLUORESCENCE_MICROSCOPY_H

