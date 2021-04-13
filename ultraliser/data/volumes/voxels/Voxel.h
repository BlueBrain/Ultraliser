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

#ifndef ULTRALISER_DATA_VOXEL_H
#define ULTRALISER_DATA_VOXEL_H

#include <common/Common.h>

/**
 * @brief The Voxel struct
 */
struct Voxel
{
public:

    /**
     * @brief Voxel
     * @param voxelValue
     */
    Voxel(const uint8_t voxelValue = 0)
    {
        value = voxelValue;
    }

public:

    /**
     * @brief value
     */
    uint8_t value;
};

/**
 * @brief Voxels
 */
typedef std::vector< Voxel > Voxels;



#endif // ULTRALISER_DATA_VOXEL_H
