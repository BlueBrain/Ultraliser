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

#include "Thinning6Iterations.h"

namespace Ultraliser
{

Thinning6Iterations::Thinning6Iterations()
{
    // Generates the 24 masks for each direction (6x24 = 144)
    _directionalMasks[0] = new DirectionalMask('u');
    _directionalMasks[1] = new DirectionalMask('d');
    _directionalMasks[2] = new DirectionalMask('n');
    _directionalMasks[3] = new DirectionalMask('s');
    _directionalMasks[4] = new DirectionalMask('e');
    _directionalMasks[5] = new DirectionalMask('w');
}

Thinning6Iterations::~Thinning6Iterations()
{
    delete _directionalMasks[0];
    delete _directionalMasks[1];
    delete _directionalMasks[2];
    delete _directionalMasks[3];
    delete _directionalMasks[4];
    delete _directionalMasks[5];
}

bool Thinning6Iterations::matches(const size_t &direction, const int8_t *subVolume)
{
    if (_directionalMasks[direction]->matches(subVolume))
        return true;
    return false;
}

}
