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

#ifndef ULTRALISER_DATA_VOLUME_GRIDINDEX_H
#define ULTRALISER_DATA_VOLUME_GRIDINDEX_H

#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The GridIndex struct
 */
struct GridIndex
{
public:

    /**
     * @brief index
     */
    int64_t index[3];

public:
    /**
     * @brief GridIndex
     * @param a
     * @param b
     * @param c
     */
    GridIndex(const int64_t &a = 0, const int64_t &b = 0, const int64_t &c = 0)
    {
        index[0] = a; index[1] = b; index[2] = c;
    }

public:
    /**
     * @brief operator =
     * @param i
     * @return
     */
    GridIndex& operator=(const GridIndex & i)
    {
        for(int64_t j = 0; j < DIMENSIONS; j++)
            index[0] = i.index[0];

        return *this;
    }

    /**
     * @brief operator []
     * @param idx
     * @return
     */
    int64_t& operator[](uint64_t idx) { return index[idx]; }

    /**
     * @brief operator []
     * @param idx
     * @return
     */
    int64_t operator[](uint64_t idx) const { return index[idx]; }

    /**
     * @brief operator <
     * @param grid
     * @return
     */
    bool operator<(const GridIndex& grid) const
    {
        for(size_t i = 0; i < DIMENSIONS; i++)
        {
            if(index[i] == grid.index[i])
                continue;

            if(index[i] < grid.index[i])
                return true;

            return false;
        }

        return false;
    }

    /**
     * @brief isValid
     * @return
     */
    bool isValid() const
    {
        return _valid;
    }

private:

    /**
     * @brief _valid
     * Means that the grid index is a valid one that can be indexed.
     */
    bool _valid;
};

}

#endif // ULTRALISER_DATA_VOLUME_GRIDINDEX_H
