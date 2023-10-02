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

#include <common/Common.h>
#include <math/Vector.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

enum PIXEL_COLOR
{
    BLACK = 0,
    GRAY = 128,
    WHITE = 255
};

/**
 * @brief The Image class
 */
class Image
{
public:

    /**
     * @brief Image
     * Constructor
     * @param width
     * The width of the image in pixels.
     * @param height
     * The height of the image in pixels.
     */
    Image(const size_t& width, const size_t& height);
    ~Image();

    /**
     * @brief getWidth
     * Returns the width of the image in pixels.
     * @return
     * The width of the image in pixels.
     */
    size_t getWidth() const { return _width; }

    /**
     * @brief getHeight
     * Returns the height of the image in pixels.
     * @return
     * The height of the image in pixels.
     */
    size_t getHeight() const { return _height; }

    /**
     * @brief getNumberPixels
     * @return
     */
    size_t getNumberPixels() const { return _numberPixels; }

    /**
     * @brief mapTo1DIndex
     * Maps the XY index of a given pixel into its corresponding 1D index.
     * @param x
     * The X-index of the pixel. Note that this value cannot be negative.
     * @param y
     * The Y-index of the pixel. Note that this value cannot be negative.
     * @return
     * The corresponding 1D index of the pixel in the array.
     */
    size_t mapTo1DIndex(const size_t &x, const size_t &y) const
    {
        if (x >= getWidth())
            LOG_ERROR("Index [%d, %d] is out of bound", x, y);
        if (y >= getHeight())
            LOG_ERROR("Index [%d, %d] is out of bound", x, y);
        return (x + (_width * y));
    }

    /**
     * @brief mapTo1DIndexWOBC
     * Maps the XY index of a given pixel into its corresponding 1D index. This function is
     * similar to @mapTo1DIndex, but it is faster becuase it does not have to verify the bounds.
     * The X-index of the pixel. Note that this value cannot be negative.
     * @param y
     * The Y-index of the pixel. Note that this value cannot be negative.
     * @return
     * The corresponding 1D index of the pixel in the array.
     */
    size_t mapTo1DIndexWOBC(const size_t &x, const size_t &y) const
    {
        return (x + (_width * y));
    }

    /**
     * @brief fill
     * Fills the entire image with a specific color.
     * @param color
     * The color with which the image will be filled.
     */
    void fill(const PIXEL_COLOR& color)
    {
        for(size_t i = 0; i < _width; ++i)
            for(int64_t j = 0; j < _height; ++j)
                setPixelColor(i, j, color);
    }

    /**
     * @brief dimension
     * @param i
     * @return
     */
    size_t dimension(const size_t &i) const;

    /**
     * @brief setPixelColor
     * @param index
     * @param color
     */
    void setPixelColor(const size_t &index, const PIXEL_COLOR& color)
    {
        _data[index] = color;
    }

    /**
     * @brief setPixelColor
     * Sets the color of a specific pixel specified by its X and Y indices.
     * @param x
     * The X-index of the pixel to be set.
     * @param y
     * The Y-index of the pixel to be set.
     * @param color
     * The color of the pixel.
     */
    void setPixelColor(const int64_t &x, const int64_t &y, const PIXEL_COLOR& color)
    {
        setPixelColor(mapTo1DIndex(x, y), color);
    }

    /**
     * @brief getPixelColor
     * Gets the color of a pixel specified by 1D index.
     * @param index
     * The 1D index of the pixel.
     * @return
     * The color of the pixel.
     */
    PIXEL_COLOR getPixelColor(const size_t &index) const
    {
        return PIXEL_COLOR(_data[ index ]);
    }

    /**
     * @brief getPixelColor
     * Gets the color of a pixel specified by its X and Y indices.
     * @param x
     * The X-index of the pixel.
     * @param y
     * The Y-index of the pixel.
     * @return
     * The color of the pixel.
     */
    PIXEL_COLOR getPixelColor(const size_t &x, const size_t &y)
    {
        return getPixelColor(mapTo1DIndex(x, y));
    }

    /**
     * @brief isFilled
     * Checks if a specific pixel is filled or not.
     * @param x
     * The X-index of the pixel.
     * @param y
     * The Y-index of the pixel.
     * @return
     * True if the pixel if filled, and false otherwise.
     */
    bool isFilled(const size_t& x, const size_t& y)
    {
        // If the given index is not within the extent of the pixel, return zero
        if(x >= getWidth())
            return false;
        if (y >= getHeight())
            return false;

        // Return true if the data is set, otherwise false
        if (_data[(x + (_width * y))])
            return true;
        return false;
    }

    /**
     * @brief writeEXR
     * @param prefix
     */
    void writeEXR(const std::string &prefix);

    /**
     * @brief writePPM
     * @param prefix
     */
    void writePPM(const std::string &prefix) const;

private:

    /**
     * @brief _allocateMemory
     */
    void _allocateMemory();

    /**
     * @brief _freeMemory
     */
    void _freeMemory();

private:

    /**
     * @brief _width
     * The width of the image in pixels.
     */
    size_t _width;

    /**
     * @brief height
     * The height of the image in pixels.
     */
    size_t _height;

    /**
     * @brief _numberPixels
     * The total number of pixels in the image.
     */
    size_t _numberPixels;

    /**
     * @brief _data
     * The array the contains the data of the image. The image is represented by an 8-bit unsigned
     * data to save space. It allows the image to have gray scale data only.
     */
    uint8_t* _data;
};

}
