/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <common/Common.h>
#include "DirectionalMask.h"

namespace Ultraliser
{

/**
 * @brief The Thinning6Iterations class
 * An implementation of the 3D thinning algorithm proposed by Kalman Palagyi  and  Attila Kuba.
 * This approach checks if the current voxel in the volume matches a set of six directional masks
 * depending on the given direction.
 * The indices of the faces are set as follows:
 *
 *     0  1  2
 *     3  4  5
 *     6  7  8
 *
 *     9 10 11
 *     12    13
 *     14 15 16
 *
 *     17 18 19
 *     20 21 22
 *     23 24 25
 *
 */
class Thinning6Iterations
{
public:

    /**
     * @brief Thinning6Iterations
     * Constructor.
     */
    Thinning6Iterations();
    ~Thinning6Iterations();

    /**
     * @brief matches
     * @param direction
     * @param subVolume
     * @return
     */
    bool matches(const size_t& direction, const int8_t* subVolume);

private:

    /**
     * @brief _directionalMasks
     * A set of six directional masks that are used to check if the voxel should be removed or not.
     */
    DirectionalMask *_directionalMasks[6];
};

}
