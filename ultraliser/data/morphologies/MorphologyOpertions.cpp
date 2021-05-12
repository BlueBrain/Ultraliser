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

#include "MorphologyOpertions.h"
#include <common/Common.h>

namespace Ultraliser
{

float computeSegmentLength(const Sample* sample1, const Sample* sample2)
{

    const Vector3f delta = sample1->getPosition() - sample2->getPosition();
    return delta.abs();
}

float computeSegmentSurfaceArea(const Sample* sample1, const Sample* sample2)
{
    // Radii
    const float r1 = sample1->getRadius();
    const float r2 = sample2->getRadius();
    const float radiiSum = r1 + r2;
    const float radiiDiff = r1 - r2;

    // Segment length
    const float length = computeSegmentLength(sample1, sample2);

    // Later area
    const float lateralArea = ULTRALISER_PIF * radiiSum * sqrt((radiiDiff * radiiDiff) + length);

    // The surface area of the segment
    const float surfaceArea = lateralArea + ULTRALISER_PIF * ((r1 * r1) + (r2 * r2));
    return surfaceArea;
}

float computeSegmentVolume(const Sample* sample1, const Sample* sample2)
{
    // Radii
    const float r1 = sample1->getRadius();
    const float r2 = sample2->getRadius();

    // Segment length
    const float length = computeSegmentLength(sample1, sample2);

    const float volume = (1.f / 3.f) * ULTRALISER_PIF * length * (r1 * r1 + r1 * r2 + r2 * r2);
    return volume;
}

void resampleSectionRelaxed(Section* section)
{

}

void resampleSectionPacked(Section* section)
{

}

void resampleSectionWithFixedStep(Section* section, const float& step)
{
    // Initially, compute section length

    // If the section length if less than the step, do not resample it and return

    // If the section length is less than double the step, add a new bi-centric sample

    // Otherwise, resample the section.
}

}
