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

#ifndef ULTRALISER_DATA_DMC_VOXEL_H
#define ULTRALISER_DATA_DMC_VOXEL_H

#include <common/Common.h>

/**
 * @brief The DMC_EDGE_SIDE enum
 */
enum DMC_EDGE_SIDE
{
    X,
    Y,
    Z,
    UNKNOWN
};

/**
 * @brief The DMCVoxel struct
 */
struct DMCVoxel
{
public:

    /**
     * @brief DMCVoxel
     * Constructor
     * @param i
     * X-coordinates
     * @param j
     * Y-coordinates
     * @param k
     * Z-coordinates
     * @param enter
     * Entering the voxel
     * @param exit
     * Exiting the voxel
     * @param edgeSide
     * Side axis
     */
    DMCVoxel(const int64_t& i, const int64_t& j, const int64_t& k,
             const bool& enter, const bool& exit, const DMC_EDGE_SIDE& edgeSide)
    {
        // Corrdinates
        x = i; y = j;  z = k;
        entering = enter; exiting = exit;
        side = edgeSide;
    }

    /**
     * @brief x
     * X-axis of the voxel.
     */
    int64_t x = -1;

    /**
     * @brief y
     * Y-axis of the voxel.
     */
    int64_t y = -1;

    /**
     * @brief z
     * Z-axis of the voxel.
     */
    int64_t z = -1;

    /**
     * @brief entering
     * Are we entering the voxel?
     */
    bool entering = false;

    /**
     * @brief exiting
     * Are we exiting the voxel?
     */
    bool exiting = false;

    /**
     * @brief side
     * DMC edge side value
     */
    DMC_EDGE_SIDE side = DMC_EDGE_SIDE::UNKNOWN;
};

/**
 * @brief DMCVoxels
 * A list (std::vector) of DMCVoxel's.
 */
typedef std::vector< DMCVoxel* > DMCVoxels;

/**
 * @brief DMCVoxelsList
 * A list (std::vector) of lists (std::vector's) of DMCVoxel's.
 */
typedef std::vector< DMCVoxels > DMCVoxelsList;

#endif // ULTRALISER_DATA_DMC_VOXEL_H
