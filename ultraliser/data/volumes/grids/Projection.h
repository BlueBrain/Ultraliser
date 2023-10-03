/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

#include <data/common/ColorMap.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/**
 * @brief saveColorMappedProjection
 * @param prefix
 * @param projection
 * @param projectionWidth
 * @param projectionHeight
 * @param minValue
 * @param maxValue
 * @param colorMap
 */
void saveColorMappedProjection(const std::string &prefix,
                               const double* projection,
                               const size_t &projectionWidth,
                               const size_t &projectionHeight,
                               const double& minValue,
                               const double& maxValue,
                               const std::string &colorMap);

/**
 * @brief saveColorMappedProjectionWithAllColorMaps
 * @param prefix
 * @param projection
 * @param projectionWidth
 * @param projectionHeight
 * @param minValue
 * @param maxValue
 */
void saveColorMappedProjectionWithAllColorMaps(const std::string &prefix,
                                               const double *projection,
                                               const size_t &projectionWidth,
                                               const size_t &projectionHeight,
                                               const double &minValue,
                                               const double &maxValue);

/**
 * @brief writeProjections
 * @param prefix
 * @param projectionImage
 * @param projectionWidth
 * @param projectionHeight
 * @param colorCodedProjections
 * @param verbose
 */
void writeProjections(const std::string &prefix,
                      const double* projectionImage,
                      const size_t &projectionWidth,
                      const size_t &projectionHeight,
                      const bool& colorCodedProjections,
                      const bool &verbose = true);

}

