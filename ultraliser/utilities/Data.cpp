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

#include "Data.h"

namespace Ultraliser
{

std::string formatSize(const size_t &size)
{
    std::stringstream stream;

    if (size < 1024)
    {
        stream << size;
    }
    else if (size < 1024 * 1024)
    {
        stream << float(size) / 1024.f << " k";
    }
    else if (size < 1024 * 1024 * 1024)
    {
        stream << float(size) / (1024.f * 1024.f) << " M";
    }
    else
    {
        stream << float(size) / (1024.f * 1024.f * 1024.f) << " G";
    }

    return stream.str();
}

void computeHstogram(const std::vector<float> data, const size_t nbins,
                     std::vector< size_t > &histogram, std::vector< float > &binEdges)
{
    // Get the minimum value
    const auto minValue = *std::min_element(data.begin(), data.end());

    // Get the maximum value
    const auto maxValue = *std::max_element(data.begin(), data.end());

    // Compute the range
    const auto range = maxValue - minValue;

    // Compute the step
    const auto step = range / (nbins);

    // Construct the histogram
    std::vector < float > xaxis;
    std::vector < size_t > yaxis;

    // Resize the elements based on the number of bins
    xaxis.resize(nbins + 1);
    yaxis.resize(nbins);

    // Initialization
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nbins; ++i)
    {
        xaxis[i] = 0.0; yaxis[i] = 0;
    }

    // Histogram calculcations
    OMP_PARALLEL_FOR
    for (size_t i = 0; i <= nbins; ++i)
    {
        // Step range, for comparison
        const auto value0 = minValue + i * step;
        const auto value1 = minValue + (i + 1) * step;

        // Set the x-axis value
        xaxis[i] = value0;

        // Set the y-axis value
        for (size_t j = 0; j < data.size(); ++j)
        {
            if (data[j] >= value0 && data[j] < value1)
            {
                yaxis[i] += 1;

            }
        }
    }

    // Return the values
    histogram = yaxis;
    binEdges = xaxis;
}

}
