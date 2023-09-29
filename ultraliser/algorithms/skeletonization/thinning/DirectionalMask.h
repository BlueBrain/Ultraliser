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

#include "Mask.h"

namespace Ultraliser
{

/**
 * @brief The DirectionalMask class
 * Store M1, M2, M3, M4, M5 and M6 masks for a given direction.
 * This directional mask contains the base mask and and the four rotations.
 * This class then includes 6 x 4 = 24 masks for the same direction.
 * NOTE: The i-th element of the base mask corresponds to the mask Mi
 */
class DirectionalMask
{
public:

    /**
     * @brief DirectionalMask
     * @param direction
     */
    DirectionalMask(char direction);
    ~DirectionalMask();

    /**
     * @brief matches
     * Does the mask match an input sub-volume?
     * @param subVolume
     * @return
     */
    bool matches(const int8_t *subVolume);

private:

    /**
     * @brief _baseMasks
     */
    Mask *_baseMasks;
};

}
