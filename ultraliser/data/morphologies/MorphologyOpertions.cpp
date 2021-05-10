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

}
