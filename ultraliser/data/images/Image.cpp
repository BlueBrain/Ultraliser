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

#include "Image.h"
#include <utilities/Utilities.h>

namespace Ultraliser
{

Image::Image(const Vec2i_64 &dimensions)
{
    _dimensions = dimensions;
    _numberPixels = I2UI64(dimensions.v[0] * dimensions.v[1]);

    _allocateMemory();
}

Image::Image(const int64_t &width, const int64_t &height)
{
    _dimensions.v[0] = width;
    _dimensions.v[1] = height;
    _numberPixels = I2UI64(_dimensions.v[0] * _dimensions.v[1]);

    _allocateMemory();
}

void Image::_allocateMemory(void)
{
    // Allocate the BIT ARRAY
    _data = new uint8_t[_numberPixels];

    // Fill black for the moment
    for (uint64_t i = 0; i < _numberPixels; ++i)
        _data[i] = BLACK;
}

void Image::_freeMemory(void)
{
    delete [] _data;
}

int64_t Image::dimension(const int& i) const
{
    switch (i)
    {
        case 0 :
            return getWidth();
        case 1:
            return getHeight();
    }

    return 0;
}

void Image::writePPM(const std::string &prefix) const
{
    std::stringstream stream;
    stream << prefix << ".ppm";
    FILE *image = fopen(stream.str().c_str(), "wb");
    fprintf(image, "P6\n%ld %ld\n255\n", getWidth(), getHeight());

    uint64_t index = 0;
    for (int64_t i = 0; i < getWidth(); ++i)
    {
        for (int64_t j = 0; j < getHeight(); ++j)
        {
            uint64_t index1D = I2UI64(getWidth() * getHeight()) - index;
            uint8_t value;

            if (PIXEL_COLOR(_data[index1D]) == WHITE)
                value = 255;
            else if (PIXEL_COLOR(_data[ index1D ]) == GRAY)
                value = 128;
            else
                value = 0;

            uint8_t color[3];
            color[0] = value; // R
            color[1] = value; // G
            color[2] = value; // B

            fwrite(color, 1, 3, image);

            index++;
        }
    }

    fclose(image);
}

Image::~Image()
{
    _freeMemory();
}

}
