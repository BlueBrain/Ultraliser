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
                               const double *projection,
                               const size_t &projectionWidth,
                               const size_t &projectionHeight,
                               const double &minValue,
                               const double &maxValue,
                               const std::string &colorMap)
{

    // Construct the colored projection array
    auto coloredProjection = std::make_unique< Vector3f[] > (projectionWidth * projectionHeight);

    // Fill the color-coded project
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (size_t index = 0; index < projectionWidth * projectionHeight; ++index)
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
                                               const double* projection,
                                               const size_t &projectionWidth,
                                               const size_t &projectionHeight,
                                               const double& minValue,
                                               const double& maxValue)
{
    // Starts the timer
    TIMER_SET;

    LOOP_STARTS("Color Coded Projections");
    for (size_t i = 0; i < ColorMap::MAPS.size(); ++i)
    {
        LOOP_PROGRESS(i, ColorMap::MAPS.size());
        saveColorMappedProjection(prefix, projection, projectionWidth, projectionHeight,
                                  minValue, maxValue, ColorMap::MAPS[i]);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void writeProjections(const std::string &prefix,
                      const double* projectionImage,
                      const size_t &projectionWidth,
                      const size_t &projectionHeight,
                      const bool &colorCodedProjections,
                      const bool &verbose)
{
    // Create normalized projection array (0 - 255)
    uint16_t* normalizedProjectionImage = new uint16_t[projectionWidth * projectionHeight]();

    // Get the maximum value
    double maxValue = 0.0;
    for (int64_t index = 0; index < projectionWidth * projectionHeight; ++index)
    {
        if (projectionImage[index] > maxValue)
            maxValue = projectionImage[index];
    }

    // Construct the normalized projection
    OMP_PARALLEL_FOR
    for (int64_t index = 0; index < projectionWidth * projectionHeight; ++index)
    {
        // Compute float pixel value
        double pixelValue = 255 * (projectionImage[index] / maxValue);

        // Convert to uint8_t to be able to write it to the image
        normalizedProjectionImage[index] = F2UI16(pixelValue);
    }

    // Save the projection into a PPM image
    Utilities::savePPMLuminanceImage(prefix, normalizedProjectionImage,
                                     projectionWidth, projectionHeight);

    // Free the normalized projections
    delete[] normalizedProjectionImage;

    // Save color coded projections with all possible color-maps
    if (colorCodedProjections)
    {
        saveColorMappedProjectionWithAllColorMaps(prefix, projectionImage,
                                                  projectionWidth, projectionHeight, 0, maxValue);
    }
}

}
