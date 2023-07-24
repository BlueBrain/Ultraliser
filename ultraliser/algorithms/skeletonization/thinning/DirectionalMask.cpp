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

#include "DirectionalMask.h"

namespace Ultraliser
{

// Templates for 27-neighborhood and  the  six masks of the up direction.
int8_t ssMu[6][26] =
{{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3 },
 { 2, 2, 2, 0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 2, 2, 2, 0, 0, 2, 0, 0, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 0, 0, 0, 3, 1, 3, 3, 0, 3, 0, 0, 0 },
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 0, 2, 2, 1, 2, 0, 0, 1, 0, 0, 2 }};

DirectionalMask::DirectionalMask(char direction)
{
    _baseMasks = new Mask[6];

    for (int i = 6; i--;)
    {
        _baseMasks[i].setDirection(direction);
        _baseMasks[i].set_mask_from_u(ssMu[i]);
        _baseMasks[i].generateRotations();
    }
}

DirectionalMask::~DirectionalMask()
{
    delete[]
    _baseMasks;
}


bool DirectionalMask::matches(const int8_t *subVolume)
{
    for (size_t i = 0; i < 6; ++i)
    {
        if (_baseMasks[i].matches(subVolume))
            return true;
    }
    return false;
}

}

