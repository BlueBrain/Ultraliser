/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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


#include "FloodFiller.h"

namespace Ultraliser
{

bool FloodFiller::isSameColor(const PIXEL_COLOR c1, const PIXEL_COLOR c2)
{
    return (c1 == c2);
}

bool FloodFiller::isSameColor(const bool c1, const bool c2)
{
    return (c1 == c2);
}

void FloodFiller::fill(Image* image,
                       const int64_t &nx, const int64_t &ny,
                       const int64_t &x, const int64_t &y,
                       PIXEL_COLOR background, PIXEL_COLOR fillColor)
{

    if (!isSameColor(background, image->getPixelColor(x, y)))
        return;

    /// Bottom
    int64_t y1 = y;
    while (y1 < ny && isSameColor(background, image->getPixelColor(x, y1)))
        image->setPixelColor(x, y1++, fillColor);

    /// Top
    y1 = y - 1;
    while (y1 >= 0 && isSameColor(background, image->getPixelColor(x, y1)))
        image->setPixelColor(x, y1--, fillColor);

    /// Left
    y1 = y;
    while (y1 < ny && isSameColor(fillColor, image->getPixelColor(x, y1)))
    {
        if (x > 0 && isSameColor(background, image->getPixelColor(x - 1, y1)))
            fill(image, nx, ny, x - 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && isSameColor(fillColor, image->getPixelColor(x, y1)))
    {
        if (x > 0 && isSameColor(background, image->getPixelColor(x - 1, y1)))
            fill(image, nx, ny, x - 1, y1, background, fillColor);
        y1--;
    }

    /// Right
    y1 = y;
    while (y1 < ny && isSameColor(fillColor, image->getPixelColor(x, y1)))
    {
        if (x < nx - 1 && isSameColor(background, image->getPixelColor(x + 1, y1)))
            fill(image, nx, ny, x + 1, y1, background, fillColor);
        y1++;
    }

    y1 = y - 1;
    while (y1 >= 0 && isSameColor(fillColor,image->getPixelColor(x, y1)))
    {
        if (x < nx - 1 && isSameColor(background, image->getPixelColor(x + 1, y1)))
            fill(image, nx, ny, x + 1, y1, background, fillColor);
        y1--;
    }
    return;
}

}
