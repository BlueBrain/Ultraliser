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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_SECTION_H
#define ULTRALISER_DATA_MORPHOLOGIES_SECTION_H

#include <data/morphologies/Sample.h>

namespace Ultraliser
{

class Section
{
public:
    Section(const uint64_t &index);

public:

    /**
     * @brief getIndex
     * @return The unique index of the section.
     */
    uint64_t getIndex() const;

    /**
     * @brief addSample
     * Adds a new sample to the samples list.
     * @param sample
     */
    void addSample(Sample* sample);

    /**
     * @brief addParentIndex
     * Adds the index of a parent section.
     * @param index
     * The index of the parent section.
     */
    void addParentIndex(const uint64_t index);

    /**
     * @brief addChildIndex
     * Adds the index of a child section.
     * @param index
     * The index of the child section.
     */
    void addChildIndex(const uint64_t index);

    /**
     * @brief getSamples
     * @return
     */
    Samples getSamples() const;

    /**
     * @brief getFirstSample
     * Gets the first sample along the section.
     * @return
     * A pointer to the first sample along the section.
     */
    Sample* getFirstSample() const;

    /**
     * @brief getLastSample
     * Gets the last sample along the section.
     * @return
     * A pointer to the last sample along the section.
     */
    Sample* getLastSample() const;

    /**
     * @brief resampleSectionUniformly
     * Resamples the section uniformy.
     * @param step
     * Resampling step.
     */
    void resampleSectionUniformly(const float& step);

    /**
     * @brief resampleSectionAdaptively
     * Resamples the section adaptively.
     */
    void resampleSectionAdaptively();

    /**
     * @brief getParentIndices
     * @return Returns a list of the indices of the parents.
     */
    std::vector< uint64_t > getParentIndices() const;

    /**
     * @brief getChildrenIndices
     * @return Returns a list of the children indices.
     */
    std::vector< uint64_t > getChildrenIndices() const;

    /**
     * @brief computeLength
     * Computes the length of the section.
     * @return
     * The length of the section.
     */
    float computeLength() const;

    /**
     * @brief computeSegmentsLengthDistribution
     * Computes the distribution of the lengths of the segments in the section.
     * @return
     */
    std::vector< float > computeSegmentsLengthDistribution() const;

    /**
     * @brief computeSurfaceArea
     * Computes the surface area of the section.
     * @return
     * The surface area of the section.
     */
    float computeSurfaceArea() const;

    /**
     * @brief computeSegmentsSurfaceAreaDistribution
     * Computes the distribution of the surface area of the segments in the section.
     * @return
     * The distribution of the surface area of the segments in the section.
     */
    std::vector< float > computeSegmentsSurfaceAreaDistribution() const;

    /**
     * @brief computeVolume
     * Computes the volume of the section.
     * @return
     * The volume of the section.
     */
    float computeVolume() const;

    /**
     * @brief computeSegmentsVolumeDistribution
     * Computes the distribution of the volume of the segments in the section.
     * @return
     * The distribution of the volume of the segments in the section.
     */
    std::vector< float > computeSegmentsVolumeDistribution() const;

    /**
     * @brief computeAverageRadius
     * Computes the average radius of the section.
     * @return
     * Computes the average radius of the section.
     */
    float computeAverageRadius() const;

private:

    /**
     * @brief _index
     * The unique index of the section.
     */
    uint64_t _index;

    /**
     * @brief _parentIndices
     * A list of the indices of the parent sections.
     */
    std::vector< uint64_t > _parentsIndices;

    /**
     * @brief _childrenIndices
     * A list of the indices of the children sections.
     */
    std::vector< uint64_t > _childrenIndices;

    /**
     * @brief _samples
     * A list of samples along a section.
     */
    Samples _samples;
};

/**
 * @brief Sections
 * A list of sections;
 */
typedef std::vector< Section* > Sections;

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_SECTION_H
