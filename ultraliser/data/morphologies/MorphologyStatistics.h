/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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
     * @brief computeNumberSamplesPerSectionDistribution
     * @return
     */
    std::vector< float > computeNumberSamplesPerSectionDistribution() const;

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