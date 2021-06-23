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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H

#include <data/morphologies/Sample.h>
#include <data/morphologies/Section.h>
#include <data/common/ROI.h>

namespace Ultraliser
{

class Morphology
{
public:

    /**
     * @brief Morphology
     * Constructor
     */
    Morphology();

    /**
     * @brief Morphology
     * Destructor
     */
    virtual ~Morphology();

public:

    /**
     * @brief getSections
     * @return Returns a reference to the list of sections in the morphology.
     */
    Sections getSections() const;

    /**
     * @brief getSamples
     * Return a reference to the list of samples in the morphology.
     * @return
     */
    Samples getSamples() const;

    /**
     * @brief getConnectedPathsFromParentsToChildren
     * @param section
     * @return
     */
    Paths getConnectedPathsFromParentsToChildren(const Section* section) const;

    /**
     * @brief getBoundingBox
     * Returns the bounding box of the morphology.
     * @param pMin
     * The computed pMin of the bounding box.
     * @param pMax
     * The computed pMax of the bounding box.
     * @param bounds
     * The computed bounds of the bounding box.
     * @param center
     * The center of the bounding box.
     */
    void getBoundingBox(Vector3f& pMin, Vector3f& pMax, Vector3f& bounds, Vector3f &center) const;

    /**
     * @brief computeTotalLength
     * Computes the total length of the morphology.
     * @return
     * The total length of the morphology.
     */
    float computeTotalLength() const;

    /**
     * @brief computeTotalSurfaceArea
     * Computes the total surface area of the morphology.
     * @return
     * The total surface area of the morphology.
     */
    float computeTotalSurfaceArea() const;

    /**
     * @brief computeTotalVolume
     * Computes the total volume of the morphology.
     * @return
     * The total volume of the morphology.
     */
    float computeTotalVolume() const;

    /**
     * @brief computeNumberSegments
     * Computes the numbers of segments in the morphology.
     * @return
     * The number of segments in the morphology.
     */
    uint64_t computeNumberSegments() const;

    /**
     * @brief collectRegionsWithThinStructures
     * Compiles a list of small structures in the morphology for the meshing.
     * @return
     */
    ROIs collectRegionsWithThinStructures() const;

    /**
     * @brief computeMinMaxAvgSampleRadius
     * Computes the minimum, maximum and average sample radius in the morphology.
     * @param minSampleRadius
     * Returned minimum sample radius.
     * @param maxSampleRadius
     * Returned maximum sample radius.
     * @param avgSampleRadius
     * Returned average sample radius.
     */
    void computeMinMaxAvgSampleRadius(float& minSampleRadius,
                                      float& maxSampleRadius,
                                      float& avgSampleRadius) const;

    /**
     * @brief computeMinMaxAvgSegmentLength
     * Computes the minimum, maximum and average segment length in the morphology.
     * @param minSegmentLength
     * Resulting minimum segment length.
     * @param maxSegmentLength
     * Resulting maximum segment length.
     * @param avgSegmentLength
     * Resulting average segment length.
     */
    void computeMinMaxAvgSegmentLength(float& minSegmentLength,
                                       float& maxSegmentLength,
                                       float& avgSegmentLength,
                                       const uint64_t &numberMorphologySegments) const;

    /**
     * @brief computeMinMaxAvgSegmentSurfaceArea
     * Computes the minimum, maximum and average segments surface areas.
     * @param minSegmentSurfaceArea
     * Resulting minimum segment surface area.
     * @param maxSegmentSurfaceArea
     * Resulting maximum segment surface area.
     * @param avgSegmentSurfaceArea
     * Resulting average segment surface area.
     */
    void computeMinMaxAvgSegmentSurfaceArea(float& minSegmentSurfaceArea,
                                            float& maxSegmentSurfaceArea,
                                            float& avgSegmentSurfaceArea,
                                            const uint64_t& numberMorphologySegments) const;

    /**
     * @brief computeMinMaxAvgSegmentVolume
     * Computes the minimum, maximum and average segments volumes.
     * @param minSegmentVolume
     * Resulting minimum segment volume.
     * @param maxSegmentVolume
     * Resulting maximum segment volume.
     * @param avgSegmentVolume
     * Resulting average segment volume.
     */
    void computeMinMaxAvgSegmentVolume(float& minSegmentVolume,
                                       float& maxSegmentVolume,
                                       float& avgSegmentVolume,
                                       const uint64_t &numberMorphologySegments) const;

    /**
     * @brief computeMinMaxAvgSectionLength
     * Computes the minimum, maximum and average segment length in the morphology.
     * @param minSectionLength
     * Resulting minimum section length.
     * @param maxSectionLength
     * Resulting maximum section length.
     * @param avgSectionLength
     * Resulting average section length.
     */
    void computeMinMaxAvgSectionLength(float& minSectionLength,
                                       float& maxSectionLength,
                                       float& avgSectionLength) const;

    /**
     * @brief computeMinMaxAvgSectionSurfaceArea
     * omputes the minimum, maximum and average sections surface areas.
     * @param minSectionSurfaceArea
     * Resulting minimum section surface area.
     * @param maxSectionSurfaceArea
     * Resulting maximum section surface area.
     * @param avgSectionSurfaceArea
     * Resulting average section surface area.
     */
    void computeMinMaxAvgSectionSurfaceArea(float& minSectionSurfaceArea,
                                            float& maxSectionSurfaceArea,
                                            float& avgSectionSurfaceArea) const;

    /**
     * @brief computeMinMaxAvgSectionVolume
     * Computes the minimum, maximum and average sections volumes.
     * @param minSectionVolume
     * Resulting minimum section volume.
     * @param maxSectionVolume
     * Resulting maximum  section volume.
     * @param avgSectionVolume
     * Resulting average section volume.
     */
    void computeMinMaxAvgSectionVolume(float& minSectionVolume,
                                       float& maxSectionVolume,
                                       float& avgSectionVolume) const;

    /**
     * @brief resampleSectionsUniformly
     * Resample each sections in the morphology uniformly with a fixed step.
     * @param step
     * A given step that will be used to sample the sections in the morphology.
     */
    void resampleSectionsUniformly(const float step);

    /**
     * @brief resampleSectionsAdaptively
     * Resample each section in the morphology adaptively, i.e. with respect to the radii of the
     * samples along the morphology.
     * @param relaxed
     * Relaxed resampling, if this flag is set to true, the resampling will be use the diameters
     * of the samples instead of their radii.
     */
    void resampleSectionsAdaptively(const bool& relaxed);

    /**
     * @brief printStats
     * Prints the morphology stats.
     * @param reference
     * Reference string for the morphology.
     * @param prefix
     * File prefix.
     */
    void printStats(const std::string &reference, const std::string *prefix) const;

    /**
     * @brief printDistributions
     * Prints the distributions of the morphology.
     * @param reference
     * Reference string for the morphology.
     * @param prefix
     * Output file prefix.
     */
    void printDistributions(const std::string *prefix) const;

protected:

    /**
     * @brief _samples
     * A list of all the actual samples of the morphology.
     */
    Samples _samples;

    /**
     * @brief _sections
     * A list of all the actual sections of the morphology.
     */
    Sections _sections;

    /**
     * @brief _pMin
     */
    Vector3f _pMin;

    /**
     * @brief _pMax
     */
    Vector3f _pMax;
};

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H
