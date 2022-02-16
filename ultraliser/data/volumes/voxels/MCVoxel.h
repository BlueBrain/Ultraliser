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

#ifndef ULTRALISER_DATA_MC_VOXEL_H
#define ULTRALISER_DATA_MC_VOXEL_H

#include <common/Common.h>

namespace Ultraliser
{

struct MCVoxel
{
    /**
     * @brief MCVoxel
     * @param i
     * @param j
     * @param k
     * @param mcIndex
     */
    MCVoxel(const int64_t& i, const int64_t& j, const int64_t& k, uint64_t mcIndex, double* v)
    {
        x = i; y = j; z = k;
        cubeIndex = mcIndex;

        vConfig[0] = v[0];
        vConfig[1] = v[1];
        vConfig[2] = v[2];
        vConfig[3] = v[3];
        vConfig[4] = v[4];
        vConfig[5] = v[5];
        vConfig[6] = v[6];
        vConfig[7] = v[7];
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
     * @brief cubeIndex
     */
    uint64_t cubeIndex;

    /**
     * @brief v
     */
    double vConfig[8];
};

/**
 * @brief MCVoxels
 */
typedef std::vector< MCVoxel* > MCVoxels;

/**
 * @brief MCVoxelsList
 */
typedef std::vector< MCVoxels > MCVoxelsList;

}

#endif // ULTRALISER_DATA_MC_VOXEL_H
