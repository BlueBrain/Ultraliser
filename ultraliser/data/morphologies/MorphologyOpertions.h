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
