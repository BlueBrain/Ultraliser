/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#pragma once

#include <data/morphologies/ProcessType.h>
#include <data/morphologies/Sample.h>

namespace Ultraliser
{

class Section;
typedef std::vector< Section* > Sections;

class Section
{
public:

    /**
     * @brief Section
     * @param index
     */
    Section(const size_t &index, const PROCESS_TYPE& type);

public:

    /**
     * @brief getType
     * @return
     */
    PROCESS_TYPE getType() const;

    /**
     * @brief setType
     * Set the section type
     * @param sectionType
     */
    void setType(PROCESS_TYPE type);

    /**
     * @brief getIndex
     * @return The unique index of the section.
     */
    size_t getIndex() const;

    /**
     * @brief setIndex
     * Set the section index
     * @param index 
     */
    void setIndex(const size_t &index);

    /**
     * @brief addSample
     * Adds a new sample to the samples list.
     * @param sample
     * A given sample.
     */
    void addSample(Sample* sample);

    /**
     * @brief addSamples
     * Adds a list of new samples to the current samples list.
     * @param samples
     * A list of samples.
     */
    void addSamples(const Samples& samples);

    /**
     * @brief addParent
     * Adds a reference to a parent section.
     * @param section
     * A pointer to the parent section.
     */
    void addParent(Section* section);

    /**
     * @brief addChild
     * Adds a reference to a child section.
     * @param section
     * A pointer to the child section.
     */
    void addChild(Section* section);

    /**
     * @brief addParentIndex
     * Adds the index of a parent section.
     * @param index
     * The index of the parent section.
     */
    void addParentIndex(const size_t index);

    /**
     * @brief clearParentIndices
     * Clear the indices of the parent sections.
     */
    void clearParentsIndices();

    /**
     * @brief addChildIndex
     * Adds the index of a child section.
     * @param index
     * The index of the child section.
     */
    void addChildIndex(const size_t index);

    /**
     * @brief clearChildrenIndices
     * Clear the indices of the child sections.
     */
    void clearChildrenIndices();

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
    void resampleUniformly(const float& step);

    /**
     * @brief resampleSectionAdaptively
     * Resamples the section adaptively.
     * @param relaxed
     */
    void resampleAdaptively(const bool& relaxed = true);




    size_t removeInnerSamples();

    size_t removeIntermediateSamples();

    size_t interpolateLongSegments();


    /**
     * @brief verifyMinimumSampleRadius
     * Ensures that all the samples along the section have a minimum radius that is greater then
     * the given one.
     * @param radius
     */
    void verifyMinimumSampleRadius(const float& radius);

    float getMinimumSampleRadius() const;


    /**
     * @brief isRoot
     * @return
     */
    bool isRoot();

    /**
     * @brief getParentIndices
     * @return Returns a list of the indices of the parents.
     */
    std::vector< size_t > getParentIndices() const;

    /**
     * @brief getChildrenIndices
     * @return Returns a list of the children indices.
     */
    std::vector< size_t > getChildrenIndices() const;

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

    /**
     * @brief computeNumberZeroRadiusSamples
     * Computes the number of zero-radius samples that have radii below the given threshold.
     * @param threshold
     * Minimum radius value.
     * @return
     * The number of zero-radius samples.
     */
    size_t computeNumberZeroRadiusSamples(const float& threshold) const;

    /**
     * @brief computeSamplingDensityPerMicron
     * @return
     */
    float computeSamplingDensityPerMicron() const;

    /**
     * @brief updateSamplesIndices
     * Update the local indices of the samples
     */
    void updateSamplesIndices();

    void compileSWCTableRecursively(Samples& samples,
                                    size_t &currentSampleIndex,
                                    const size_t& branchingSampleIndex = 0);


    /**
     * @brief getBranchingOrder
     * @return
     */
    size_t getBranchingOrder() const;

private:

    /**
     * @brief _type
     * The type of the section, for example: axon, basal or apical dendrite, vasculature, etc ...
     */
    const PROCESS_TYPE _type;

    /**
     * @brief _index
     * The unique index of the section.
     */
    size_t _index;

    /**
     * @brief _parentIndices
     * A list of the indices of the parent sections.
     */
    std::vector< size_t > _parentsIndices;

    /**
     * @brief _parents
     * A list of all the parent sections.
     */
    Sections _parents;

    /**
     * @brief _childrenIndices
     * A list of the indices of the children sections.
     */
    std::vector< size_t > _childrenIndices;

    /**
     * @brief _children
     * A list of all the child sections.
     */
    Sections _children;

    /**
     * @brief _samples
     * A list of samples along a section.
     */
    Samples _samples;

    /**
     * @brief _branchingOrder
     * The branching order of the section in the arbor tree.
     * NOTE: This parameter is only for acylic graphs, such as neurons and astrocytes. For
     * vasculature, it will always be 0.
     * NOTE: If this parameter is set to 0, then it is not initialized. It must be initialized
     * during the construction from the morphology file. Since sections cannot be removed from the
     * morphology, it must stay the same during the lifetime of the application.
     */
    size_t _branchingOrder;
};

/**
 * @brief Sections
 * A list of sections;
 */
typedef std::vector< Section* > Sections;

}
