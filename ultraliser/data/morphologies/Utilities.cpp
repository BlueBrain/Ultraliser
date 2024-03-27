#include "Utilities.h"

namespace Ultraliser
{

void computeSectionsBoundingBox(Sections& sections,
                                Vector3f& pMin, Vector3f& pMax)
{

    // Bounding box data
    pMin = Vector3f(std::numeric_limits<float>::max());
    pMax = Vector3f(-1 * std::numeric_limits<float>::max());

    for (const auto& section: sections)
    {
        for (const auto& sample: section->getSamples())
        {
            const auto position = sample->getPosition();
            const auto radius = sample->getRadius();

            Vector3f pMaxSample = position + Vector3f(radius);
            Vector3f pMinSample = position - Vector3f(radius);

            if (pMaxSample.x() > pMax.x()) pMax.x() = pMaxSample.x();
            if (pMaxSample.y() > pMax.y()) pMax.y() = pMaxSample.y();
            if (pMaxSample.z() > pMax.z()) pMax.z() = pMaxSample.z();

            if (pMinSample.x() < pMin.x()) pMin.x() = pMinSample.x();
            if (pMinSample.y() < pMin.y()) pMin.y() = pMinSample.y();
            if (pMinSample.z() < pMin.z()) pMin.z() = pMinSample.z();
        }
    }
}

}
