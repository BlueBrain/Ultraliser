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

#ifndef MORPHOLOGYOPERTIONS_H
#define MORPHOLOGYOPERTIONS_H

#include <data/morphologies/Sample.h>
#include <data/morphologies/Section.h>

namespace Ultraliser
{

/**
 * @brief computeSegmentLength
 * Computes the length of a given segment.
 * @param sample1
 * Sample 1 of the segment.
 * @param sample2
 * Sample 2 of the segment.
 * @return
 * The length of the segment.
 */
float computeSegmentLength(const Sample* sample1, const Sample* sample2);

/**
 * @brief computeSegmentSurfaceArea
 * Computes the surface area of a given segment.
 * @param sample1
 * Sample 1 of the segment.
 * @param sample2
 * Sample 2 of the segment.
 * @return
 * The surface area of the segment.
 */
float computeSegmentSurfaceArea(const Sample* sample1, const Sample* sample2);

/**
 * @brief computeSegmentVolume
 * Computes the volume of a given segment.
 * @param sample1
 * Sample 1 of the segment.
 * @param sample2
 * Sample 2 of the segment.
 * @return
 * The volume of the segment.
 */
float computeSegmentVolume(const Sample* sample1, const Sample* sample2);

}

#endif // MORPHOLOGYOPERTIONS_H
