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

#include "Projection.h"

namespace Ultraliser
{

void saveColorMappedProjection(const std::string &prefix,
                               const float* projection,
                               const uint64_t &projectionWidth,
                               const uint64_t &projectionHeight,
                               const float& minValue,
                               const float& maxValue,
                               const std::string &colorMap)
{

    // Construct the colored projection array
    auto coloredProjection = std::make_unique< Vector3f[] > (projectionWidth * projectionHeight);

    // Fill the color-coded project
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t index = 0; index < projectionWidth * projectionHeight; ++index)
    {
        coloredProjection[index] =
                ColorMap::getRGBColorF(projection[index], minValue, maxValue, colorMap);
    }

    // Write the colormap to a PPM image
    std::stringstream stream;

    stream << prefix << "-" << colorMap;
    Utilities::savePPMColoredImage(stream.str(), coloredProjection.get(),
                                   projectionWidth, projectionHeight);
}


void saveColorMappedProjectionWithAllColorMaps(const std::string &prefix,
                                               const float* projection,
                                               const uint64_t &projectionWidth,
                                               const uint64_t &projectionHeight,
                                               const float& minValue,
                                               const float& maxValue)
{
    // Starts the timer
    TIMER_SET;

    LOOP_STARTS("Color Coded Projections");
    for (uint64_t i = 0; i < ColorMap::MAPS.size(); ++i)
    {
        LOOP_PROGRESS(i, ColorMap::MAPS.size());
        saveColorMappedProjection(prefix, projection, projectionWidth, projectionHeight,
                                  minValue, maxValue, ColorMap::MAPS[i]);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

}
