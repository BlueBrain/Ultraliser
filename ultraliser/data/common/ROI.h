/***************************************************************************************************
 * Copyright (c) 202
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

#ifndef ULTRALISER_DATA_COMMON_ROI_H
#define ULTRALISER_DATA_COMMON_ROI_H

#include <geometry/Geometry.h>

namespace Ultraliser
{

/**
 * @brief The ROI struct
 * A region of interest defines a sphere where we need to perform local operations.
 */
struct ROI
{
    /**
     * @brief ROI
     * Constructor.
     *
     * @param center
     * The center of the ROI.
     * @param radius
     * The radius of the ROI.
     */
    ROI(const Vector3f& center, const float radius)
    {
        this->center = center;
        this->radius = radius;
    }

public:

    /**
     * @brief center
     * The center of the ROI.
     */
    Vector3f center;

    /**
     * @brief radius
     * The radius of the ROI.
     */
    float radius;
};

/**
 * @brief ROIs
 * A vector of multiple region of interests.
 */
typedef std::vector< ROI* > ROIs;

}

#endif // ULTRALISER_DATA_COMMON_ROI_H