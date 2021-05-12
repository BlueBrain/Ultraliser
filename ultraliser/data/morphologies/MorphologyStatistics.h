#ifndef MORPHOLOGYSTATISTICS_H
#define MORPHOLOGYSTATISTICS_H

#include <common/Common.h>
#include <data/morphologies/Morphologies.h>

namespace Ultraliser
{

/**
 * @brief The MorphologyStatistics class
 */
class MorphologyStatistics
{
public:
    MorphologyStatistics(const Morphology* morphology);

public:

    /**
     * @brief computeSamplesRadiiDistribution
     * @return
     */
    std::vector< float > computeSamplesRadiiDistribution() const;

    /**
     * @brief computeSectionAverageRadiiDistribution
     * @return
     */
    std::vector< float > computeSectionAverageRadiiDistribution() const;

    /**
     * @brief computeSegmentsLengthDistribution
     * @return
     */
    std::vector< float > computeSegmentsLengthDistribution() const;

    /**
     * @brief computeSectionsLengthDistribution
     * @return
     */
    std::vector< float > computeSectionsLengthDistribution() const;

    /**
     * @brief computeSegmentsSurfaceAreaDistribution
     * @return
     */
    std::vector< float > computeSegmentsSurfaceAreaDistribution() const;

    /**
     * @brief computeSectionsSurfaceAreaDistribution
     * @return
     */
    std::vector< float > computeSectionsSurfaceAreaDistribution() const;

    /**
     * @brief computeSegmentsVolumeDistribution
     * @return
     */
    std::vector< float > computeSegmentsVolumeDistribution() const;

    /**
     * @brief computeSectionsVolumeDistribution
     * @return
     */
    std::vector< float > computeSectionsVolumeDistribution() const;

    /**
     * @brief writeStatsDistributions
     * @param reference
     */
    void writeStatsDistributions(const std::string &prefix);

private:

    /**
     * @brief _morphology
     */
    const Morphology* _morphology;
};

}

#endif // MORPHOLOGYSTATISTICS_H
