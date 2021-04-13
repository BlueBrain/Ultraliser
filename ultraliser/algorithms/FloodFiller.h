/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_ALGORITHMS_FLOOD_FILLER_H
#define ULTRALISER_ALGORITHMS_FLOOD_FILLER_H

#include <data/images/Image.h>

namespace Ultraliser
{

/**
 * @brief The FloodFiller class
 */
class FloodFiller
{
public:

    /**
     * @brief isSameColor
     * This function checks if the colors of two pixels match or not.
     *
     * @param c1
     * The color of the first pixel.
     * @param c2
     * The color of the second pixel.
     * @return
     * True if they match and false otherwise.
     */
    static bool isSameColor(const PIXEL_COLOR c1, const PIXEL_COLOR c2);

    /**
     * @brief isSameColor
     * This function checks if the filling of two pixels match or not.
     *
     * @param c1
     * The value of the pixel, either true or false.
     * @param c2
     * The value of the pixel, either true or false.
     * @return
     * True if they match and false otherwise.
     */
    static bool isSameColor(const bool c1, const bool c2);

    /**
     * @brief fill
     * Flood-fill the image.
     *
     * @param image
     * An image to be flood-filled.
     * @param nx
     * Number of pixels on the x-axis.
     * @param ny
     * Number of pixels on the y-axis.
     * @param x
     * Seed pixel on the x-axis.
     * @param y
     * Seed pixel on the y-axis.
     * @param background
     * The default background color.
     * @param fillColor
     * The filling color.
     */
    static void fill(Image* image, const int64_t &nx, const int64_t &ny,
                     const int64_t &x, const int64_t &y,
                     PIXEL_COLOR background, PIXEL_COLOR fillColor);
};

}

#endif // ULTRALISER_ALGORITHMS_FLOOD_FILLER_H
