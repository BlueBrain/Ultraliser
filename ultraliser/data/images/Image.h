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
     * @param dimensions
     */
    Image(const Vec2i_64& dimensions);

    /**
     * @brief Image
     * @param width
     * @param height
     */
    Image(const int64_t &width, const int64_t &height);
    ~Image();

    /**
     * @brief getWidth
     * @return
     */
    int64_t getWidth() const { return _dimensions.v[0]; }

    /**
     * @brief getHeight
     * @return
     */
    int64_t getHeight() const { return _dimensions.v[1]; }

    /**
     * @brief getNumberPixels
     * @return
     */
    uint64_t getNumberPixels() const { return _numberPixels; }

    /**
     * @brief mapToIndex
     * @param x
     * @param y
     * @return
     */
    uint64_t mapToIndex(const int64_t &x,
                        const int64_t &y) const
    {

        if(x >= getWidth()  || x < 0 || y >= getHeight() || y < 0)
            LOG_ERROR("Index [%d, %d] is out of bound", x, y);
        return I2UI64((x + (_dimensions.v[0] * y)));
    }

    /**
     * @brief fill
     * @param color
     */
    void fill(PIXEL_COLOR color)
    {
        for(int64_t i = 0; i < getWidth(); i++)
            for(int64_t j = 0; j < this->getHeight(); j++)
                setPixelColor(i, j, color);
    }

    /**
     * @brief dimension
     * @param i
     * @return
     */
    int64_t dimension(const int& i) const;

    /**
     * @brief setPixelColor
     * @param index
     * @param color
     */
    void setPixelColor(const uint64_t &index, const PIXEL_COLOR& color)
    {
        _data[ index ] = color;
    }

    /**
     * @brief setPixelColor
     * @param x
     * @param y
     * @param color
     */
    void setPixelColor(const int64_t &x,
                       const int64_t &y,
                       const PIXEL_COLOR& color)
    {
        setPixelColor(mapToIndex(x, y), color);
    }

    /**
     * @brief getPixelColor
     * @param index
     * @return
     */
    PIXEL_COLOR getPixelColor(const uint64_t &index) const
    {
        return PIXEL_COLOR(_data[ index ]);
    }

    /**
     * @brief getPixelColor
     * @param x
     * @param y
     * @return
     */
    PIXEL_COLOR getPixelColor(const int64_t &x, const int64_t &y)
    {
        return getPixelColor(mapToIndex(x, y));
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
    void _allocateMemory(void);

    /**
     * @brief _freeMemory
     */
    void _freeMemory(void);

private:

    /**
     * @brief _dimensions
     */
    Vec2i_64 _dimensions;

    /**
     * @brief _numberPixels
     */
    uint64_t _numberPixels;

    /**
     * @brief _data
     */
    u_int8_t* _data;
};

}
