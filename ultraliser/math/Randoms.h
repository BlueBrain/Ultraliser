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
#include <math/Vector4f.h>

namespace Ultraliser
{

/**
 * @brief getRandomAngle
 * Gets a random angle in 2D.
 * @return
 */
float getRandomAngle(void);

/**
 * @brief getRandomInt
 * Gets a random integer between two values.
 * @param min
 * @param max
 * @return
 */
size_t getRandomInt(const size_t &min, const size_t &max);

/**
 * @brief getRandomFloat
 * Gets a random float between two values.
 * @param min
 * @param max
 * @return
 */
float getRandomFloat(const float min = 0.f, const float max = 1.0);

/**
 * @brief getRandomRGBAColor
 * Gets a random RGBA color between two values.
 * @param min
 * @param max
 * @return
 */
Vector4f getRandomRGBAColor(const float min = 0.f, const float max = 1.0);

}
