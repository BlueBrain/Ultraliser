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

#include "Morphology.h"
#include "MorphologyOpertions.h"
#include "MorphologyStatistics.h"
#include <common/Common.h>
#include <utilities/Utilities.h>


namespace Ultraliser
{

Morphology::Morphology()
{
    /// EMPTY CONSTRUCTOR
}

Morphology::~Morphology()
{
    for (auto section: _sections)
    {
        delete section;
    }
    _sections.clear();
    for (auto sample: _samples)
    {
        delete sample;
    }
    _samples.clear();
}

Sections Morphology::getSections() const
{
    return _sections;
}

Samples Morphology::getSamples() const
{
    return _samples;
}

Sections Morphology::getRootSections() const
{
    // Collect all the root sections in this list
    Sections rootSections;
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        // The root section has no parents, indeed
        if (_sections[i]->getParentIndices().size() == 0)
            rootSections.push_back(_sections[i]);
    }

    // Return a list with all the root nodes
    return rootSections;
}

Sections Morphology::getSubsectionsInBoundingBox(const Section* section,
                                                 const Vector3f& center,
                                                 const float& width,
                                                 const float& height,
                                                 const float& depth) const
{
    // Collecting the list of sections that are located within the bounding box
    Sections internalSections;
    Samples internalSamples;
    size_t internalSectionIndex = 0;

    // Get all the samples in the section
    const Samples samples = section->getSamples();

    size_t sampleIndex = 0;
    while (1)
    {
        // If the samples is located within the bounding box
        if (samples[sampleIndex]->isLocatedInBoundingBox(center, width, height, depth))
        {
            // Add the sample to the list, and proceed
            internalSamples.push_back(samples[sampleIndex]);
        }

        // If the sample is not located within the bounding box
        else
        {
            // If the internalSamples has some collected samples, construct it
            if (internalSamples.size() > 1)
            {
                // Create a new internal section
                Section* internalSection = new Section(internalSectionIndex, section->getType());
                internalSectionIndex++;

                // Add the samples to the section
                for (size_t i = 0; i < internalSamples.size(); ++i)
                {
                    internalSection->addSample(internalSamples[i]);
                }

                // Append the section to the internal sections list
                internalSections.push_back(internalSection);
            }

            // Clear the internalSamples to collect the rest of the section
            internalSamples.clear();


        }

        // Increment the sample index
        sampleIndex++;

        if (sampleIndex >= samples.size())
        {
            // If the internalSamples has some collected samples, construct it
            if (internalSamples.size() > 1)
            {
                // Create a new internal section
                Section* internalSection = new Section(internalSectionIndex, section->getType());
                internalSectionIndex++;

                // Add the samples to the section
                for (size_t i = 0; i < internalSamples.size(); ++i)
                {
                    internalSection->addSample(internalSamples[i]);
                }

                // Append the section to the internal sections list
                internalSections.push_back(internalSection);
            }

            // Clear the internalSamples to collect the rest of the section
            internalSamples.clear();

            // Escape the loop
            break;
        }
    }

    // Return all the internal sections
    return internalSections;
}

Paths Morphology::getConnectedPathsFromParentsToChildren(const Section* section) const
{
    // All possible combination of paths along the section (from parents to children)
    Paths paths;

    // Get the parents data
    const std::vector< size_t > parentsIndices = section->getParentIndices();
    const size_t numberParents = parentsIndices .size();

    // Get the children data
    const std::vector< size_t > childrenIndices = section->getChildrenIndices();
    const size_t numberChildren = childrenIndices.size();

    // The samples along the section
    Samples sectionSamples = section->getSamples();

    // If the section is root and leaf at the same moment,
    if (numberParents == 0 && numberChildren == 0)
    {
        // The samples along the path from a parent section to
        Samples pathSamples;

        // Add the samples of the section
        pathSamples.insert(pathSamples.begin(), sectionSamples.begin(), sectionSamples.end());

        // Only a single path
        paths.push_back(pathSamples);

        // Done
        return paths;
    }

    // If the section is root and has children
    else if (numberParents == 0 && numberChildren > 0)
    {
        for (size_t k = 0; k < numberChildren; ++k)
        {
            // The samples along the path from a parent section to
            Samples pathSamples;

            // Get child data
            Section* childSection = _sections[childrenIndices[k]];
            Samples childSamples = childSection->getSamples();

            // Section samples
            pathSamples.insert(pathSamples.begin(),
                               sectionSamples.begin(),
                               sectionSamples.end());

            // Child section samples
            pathSamples.insert(pathSamples.end(),
                               childSamples.begin() + 1,
                               childSamples.end());

            // Add the path
            paths.push_back(pathSamples);
        }

        // Done
        return paths;
    }

    // If the section is a leaf (has no children) node but has parents
    else if (numberParents > 0 && numberChildren == 0)
    {
        for (size_t j = 0; j < numberParents; ++j)
        {
            // The samples along the path from a parent section to
            Samples pathSamples;

            // Get parent data
            Section* parentSection = _sections[parentsIndices[j]];
            Samples parentSamples = parentSection->getSamples();

            // Parent section samples
            pathSamples.insert(pathSamples.begin(),
                               parentSamples.begin(),
                               parentSamples.end());

            // Section samples
            pathSamples.insert(pathSamples.end(),
                               sectionSamples.begin() + 1,
                               sectionSamples.end());

            // Add the path
            paths.push_back(pathSamples);
        }

        // Done
        return paths;
    }

    // If the section has parents and children
    else if (numberParents > 0 && numberChildren > 0)
    {
        for (size_t j = 0; j < numberParents; ++j)
        {
            // Parents data
            Section* parentSection = _sections[parentsIndices[j]];
            Samples parentSamples = parentSection->getSamples();

            for (size_t k = 0; k < numberChildren; ++k)
            {
                // The samples along the path from a parent section to
                Samples pathSamples;

                // Children data
                Section* childSection = _sections[childrenIndices[k]];
                Samples childSamples = childSection->getSamples();

                // Parent section samples
                pathSamples.insert(pathSamples.begin(),
                                   parentSamples.begin(),
                                   parentSamples.end());

                // Ignore the first and last samples to avoid duplicated
                // samples along the poly-line. This is mandatory to be
                // able to draw a VTK poly-line.
                pathSamples.insert(pathSamples.end(),
                                   sectionSamples.begin() + 1,
                                   sectionSamples.end() - 1);

                // Child section samples
                pathSamples.insert(pathSamples.end(),
                                   childSamples.begin(),
                                   childSamples.end());

                // Add the path
                paths.push_back(pathSamples);
            }
        }

        // Done
        return paths;
    }

    // Return the paths
    return paths;
}

void Morphology::getBoundingBox(Vector3f& pMin,
                                Vector3f& pMax,
                                Vector3f& bounds,
                                Vector3f& center) const
{
    pMin = _pMin;
    pMax = _pMax;
    bounds = _pMax - _pMin;
    center = _pMin + (0.5f * bounds);
}

float Morphology::computeTotalLength() const
{
    auto length = 0.f;
    for (const auto section: _sections)
    {
        length += section->computeLength();
    }
    return length;
}

float Morphology::computeTotalSurfaceArea() const
{
    auto area = 0.f;
    for (const auto section: _sections)
    {
        area += section->computeSurfaceArea();
    }
    return area;
}

float Morphology::computeTotalVolume() const
{
    auto volume = 0.f;
    for (const auto section: _sections)
    {
        volume += section->computeVolume();
    }
    return volume;
}

void Morphology::computeMinMaxAvgSampleRadius(float& minSampleRadius,
                                              float& maxSampleRadius,
                                              float& avgSampleRadius) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Sample* sample: _samples)
    {
        const auto sampleRadius = sample->getRadius();

        if (sampleRadius < minValue)
            minValue = sampleRadius;

        if (sampleRadius > maxValue)
            maxValue = sampleRadius;

        avgValue += sampleRadius;
    }

    // Return the values
    minSampleRadius = minValue;
    maxSampleRadius = maxValue;
    avgSampleRadius = avgValue / _samples.size();
}

size_t Morphology::computeNumberSegments() const
{
    size_t numberSegments = 0;
    for (Section* section : _sections)
    {
        // Note that if the section has N samples, it has N - 1 segments
        numberSegments += section->getSamples().size() - 1;
    }
    return numberSegments;
}

void Morphology::computeMinMaxAvgSegmentLength(float& minSegmentLength,
                                               float& maxSegmentLength,
                                               float& avgSegmentLength,
                                               const size_t &numberMorphologySegments) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        for (size_t i = 0; i < section->getSamples().size(); ++i)
        {
            const auto segmentLength = computeSegmentLength(_samples[i], _samples[i + 1]);

            if (segmentLength < minValue)
                minValue = segmentLength;

            if (segmentLength > maxValue)
                maxValue = segmentLength;

            avgValue += segmentLength;
        }
    }

    // Return the values
    minSegmentLength = minValue;
    maxSegmentLength = maxValue;
    avgSegmentLength = avgValue / numberMorphologySegments;
}

void Morphology::computeMinMaxAvgSectionLength(float& minSectionLength,
                                               float& maxSectionLength,
                                               float& avgSectionLength) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        const auto sectionLength = section->computeLength();

        if (sectionLength < minValue)
            minValue = sectionLength;

        if (sectionLength > maxValue)
            maxValue = sectionLength;

        avgValue += sectionLength;
    }

    // Return the values
    minSectionLength = minValue;
    maxSectionLength = maxValue;
    avgSectionLength = avgValue / _sections.size();
}

void Morphology::computeMinMaxAvgSegmentSurfaceArea(float& minSegmentSurfaceArea,
                                                    float& maxSegmentSurfaceArea,
                                                    float& avgSegmentSurfaceArea,
                                                    const size_t &numberMorphologySegments) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        for (size_t i = 0; i < section->getSamples().size(); ++i)
        {
            const auto segmentSurfaceArea = computeSegmentSurfaceArea(_samples[i], _samples[i + 1]);

            if (segmentSurfaceArea < minValue)
                minValue = segmentSurfaceArea;

            if (segmentSurfaceArea > maxValue)
                maxValue = segmentSurfaceArea;

            avgValue += segmentSurfaceArea;
        }
    }

    // Return the values
    minSegmentSurfaceArea = minValue;
    maxSegmentSurfaceArea = maxValue;
    avgSegmentSurfaceArea = avgValue / numberMorphologySegments;
}

void Morphology::computeMinMaxAvgSectionSurfaceArea(float& minSectionSurfaceArea,
                                                    float& maxSectionSurfaceArea,
                                                    float& avgSectionSurfaceArea) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        const auto sectionSurfaceArea = section->computeSurfaceArea();

        if (sectionSurfaceArea < minValue)
            minValue = sectionSurfaceArea;

        if (sectionSurfaceArea > maxValue)
            maxValue = sectionSurfaceArea;

        avgValue += sectionSurfaceArea;
    }

    // Return the values
    minSectionSurfaceArea = minValue;
    maxSectionSurfaceArea = maxValue;
    avgSectionSurfaceArea = avgValue / _sections.size();
}


void Morphology::computeMinMaxAvgSegmentVolume(float& minSegmentVolume,
                                               float& maxSegmentVolume,
                                               float& avgSegmentVolume,
                                               const size_t& numberMorphologySegments) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        for (size_t i = 0; i < section->getSamples().size(); ++i)
        {
            const auto segmentVolume = computeSegmentVolume(_samples[i], _samples[i + 1]);

            if (segmentVolume < minValue)
                minValue = segmentVolume;

            if (segmentVolume > maxValue)
                maxValue = segmentVolume;

            avgValue += segmentVolume;
        }
    }

    // Return the values
    minSegmentVolume = minValue;
    maxSegmentVolume = maxValue;
    avgSegmentVolume = avgValue / numberMorphologySegments;
}

void Morphology::computeMinMaxAvgSectionVolume(float& minSectionVolume,
                                               float& maxSectionVolume,
                                               float& avgSectionVolume) const
{
    float minValue = std::numeric_limits<float>::max();
    float maxValue = std::numeric_limits<float>::lowest();
    float avgValue = 0.f;

    for (const Section* section : _sections)
    {
        const auto sectionVolume = section->computeVolume();

        if (sectionVolume < minValue)
            minValue = sectionVolume;

        if (sectionVolume > maxValue)
            maxValue = sectionVolume;

        avgValue += sectionVolume;
    }

    // Return the values
    minSectionVolume = minValue;
    maxSectionVolume = maxValue;
    avgSectionVolume = avgValue / _sections.size();
}

void Morphology::verifyMinimumSampleRadius(const float& radius)
{
    LOG_STATUS("Verifying Minimum Sample Radius");

    // Starting the timer
    TIMER_SET;

    LOOP_STARTS("Checking Morphology")
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        _sections[i]->verifyMinimumSampleRadius(radius);

        LOOP_PROGRESS(PROGRESS, _sections.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Morphology::resampleSectionsUniformly(const float step)
{
    LOG_STATUS("Resampling Morphology");

    // Starting the timer
    TIMER_SET;

    LOOP_STARTS("Uniform Resampling")
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        _sections[i]->resampleUniformly(step);

        LOOP_PROGRESS(PROGRESS, _sections.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Morphology::resampleSectionsAdaptively(const bool& relaxed)
{
    LOG_STATUS("Resampling Morphology");

    // Starting the timer
    TIMER_SET;

    LOOP_STARTS("Adaptive Resampling")
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        _sections[i]->resampleAdaptively(relaxed);

        LOOP_PROGRESS(PROGRESS, _sections.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Morphology::resampleSectionsSmartly()
{
    LOG_STATUS("Resampling Morphology");

    // Starting the timer
    TIMER_SET;

    std::vector< size_t > sectionCounter;
    sectionCounter.resize(_sections.size());

    LOOP_STARTS("Removing Redundant Samples")
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        sectionCounter[i] = _sections[i]->removeInnerSamples();

        LOOP_PROGRESS(PROGRESS, _sections.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    size_t totalNumberRemovedSamples = 0;
    for (size_t i = 0; i < sectionCounter.size(); ++i)
        totalNumberRemovedSamples += sectionCounter[i];
    LOG_SUCCESS("Number of Redundant Sample(s) [%d]", totalNumberRemovedSamples);

    LOOP_STARTS("Intepolating Long Sections")
    PROGRESS_RESET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        sectionCounter[i] = _sections[i]->interpolateLongSegments();

        LOOP_PROGRESS(PROGRESS, _sections.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    size_t totalNumberInterpolatedSamples = 0;
    for (size_t i = 0; i < sectionCounter.size(); ++i)
        totalNumberInterpolatedSamples += sectionCounter[i];
    LOG_SUCCESS("Number of Added Sample(s) [%d]", totalNumberInterpolatedSamples);
}

ROIs Morphology::collectRegionsWithThinStructures(const float& threshold) const
{
    // Regions of interest
    ROIs regions;

    // Starting the timer
    TIMER_SET;

    LOOP_STARTS("Collecting ROIs");
    PROGRESS_SET;
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        // Update the progress bar
        LOOP_PROGRESS_FRACTION(PROGRESS, _sections.size());
        PROGRESS_UPDATE;

        const Samples& samples = _sections[i]->getSamples();
        for (size_t j = 0; j < samples.size(); ++j)
        {
            if (samples[j]->getRadius() < threshold)
            {
                // NOTE: We scale the region by 1.15 to guarantee that the triangles will be covered
                regions.push_back(new ROI(samples[j]->getPosition(), samples[j]->getRadius() * 2.0));
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_SUCCESS("[ %d ] ROIs selected", regions.size());

    // Return the regions of interest
    return regions;
}


void Morphology::printDistributions(const std::string *prefix) const
{
    // Starting the timer
    TIMER_SET;

    LOG_TITLE("Morphology Distributions");

    LOG_STATUS("Collecting Stats.");

    std::unique_ptr< MorphologyStatistics > stats = std::make_unique< MorphologyStatistics >(this);
    stats->writeStatsDistributions(*prefix);

    LOG_STATUS_IMPORTANT("Gathering Morphology Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Morphology::printStats(const std::string &reference, const std::string *prefix) const
{
    // Starting the timer
    TIMER_SET;

    LOG_TITLE("Morphology Statistics");

    LOG_STATUS("Collecting Stats.");

    // Calculate the stats.
    const float length = computeTotalLength();
    const float area = computeTotalSurfaceArea();
    const float volume = computeTotalVolume();

    const size_t numberSegments = computeNumberSegments();

    float minSampleRadius, maxSampleRadius, avgSampleRadius;
    computeMinMaxAvgSampleRadius(minSampleRadius, maxSampleRadius, avgSampleRadius);

    float minSegmentLength, maxSegmentLength, avgSegmentLength;
    computeMinMaxAvgSegmentLength(minSegmentLength, maxSegmentLength, avgSegmentLength,
                                  numberSegments);

    float minSectionLength, maxSectionLength, avgSectionLength;
    computeMinMaxAvgSectionLength(minSectionLength, maxSectionLength, avgSectionLength);

    float minSegmentArea, maxSegmentArea, avgSegmentArea;
    computeMinMaxAvgSegmentSurfaceArea(minSegmentArea, maxSegmentArea, avgSegmentArea,
                                       numberSegments);

    float minSectionArea, maxSectionArea, avgSectionArea;
    computeMinMaxAvgSectionSurfaceArea(minSectionArea, maxSectionArea, avgSectionArea);

    float minSegmentVolume, maxSegmentVolume, avgSegmentVolume;
    computeMinMaxAvgSegmentVolume(minSegmentVolume, maxSegmentVolume, avgSegmentVolume,
                                  numberSegments);

    float minSectionVolume, maxSectionVolume, avgSectionVolume;
    computeMinMaxAvgSectionVolume(minSectionVolume, maxSectionVolume, avgSectionVolume);

    Vector3f pMin, pMax, bounds, center;
    getBoundingBox(pMin, pMax, bounds, center);

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + MORPHOLOGY_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
        fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                F2D(pMin.x()), F2D(pMin.y()), F2D(pMin.z()));
        fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                F2D(pMax.x()), F2D(pMax.y()), F2D(pMax.z()));
        fprintf(info, "\t* Morphology Length     | %f \n",
                F2D(length));
        fprintf(info, "\t* Number Samples        | %s \n",
                FORMAT(_samples.size()));
        fprintf(info, "\t* Number Sections       | %s \n",
                FORMAT(_sections.size()));
        fprintf(info, "\t* Surface Area          | %f² \n",
                F2D(area));
        fprintf(info, "\t* Volume                | %f³ \n",
                F2D(volume));

        // Close the file
        fclose(info);
    }

    LOG_STATUS_IMPORTANT("Morphology Stats. [ %s ]", reference.c_str());
    LOG_INFO("\t* Bounding Box:         | [%.5f, %.5f, %.5f]",
             F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
    LOG_INFO("\t* pMin:                 | [%.5f, %.5f, %.5f]",
             F2D(pMin.x()), F2D(pMin.y()), F2D(pMin.z()));
    LOG_INFO("\t* pMax:                 | [%.5f, %.5f, %.5f]",
             F2D(pMax.x()), F2D(pMax.y()), F2D(pMax.z()));
    LOG_INFO("\t* Number Samples        | %s",
             FORMAT(_samples.size()));
    LOG_INFO("\t* Number Sections       | %s",
             FORMAT(_sections.size()));
    LOG_INFO("\t* Number Segments       | %s",
             FORMAT(numberSegments));
    LOG_INFO("\t* Morphology Length     | %f",
             F2D(length));
    LOG_INFO("\t* Surface Area          | %f²",
             F2D(area));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volume));
    LOG_INFO("\t* Distributions:");
    LOG_INFO("\t* Samples Radii         | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSampleRadius), F2D(maxSampleRadius), F2D(avgSampleRadius));
    LOG_INFO("\t* Segments Lengths      | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSegmentLength), F2D(maxSegmentLength), F2D(avgSegmentLength));
    LOG_INFO("\t* Sections Lengths      | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSectionLength), F2D(maxSectionLength), F2D(avgSectionLength));

    LOG_INFO("\t* Segments Surf. Areas  | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSegmentArea), F2D(maxSegmentArea), F2D(avgSegmentArea));
    LOG_INFO("\t* Sections Surf. Areas  | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSectionArea), F2D(maxSectionArea), F2D(avgSectionArea));
    LOG_INFO("\t* Segments Volumes      | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSegmentVolume), F2D(maxSegmentVolume), F2D(avgSegmentVolume));
    LOG_INFO("\t* Sections Volumes      | [Min.: %.5f, Max.: %.5f, Avg.: %.5f]",
             F2D(minSectionVolume), F2D(maxSectionVolume), F2D(avgSectionVolume));
    LOG_INFO("");

    LOG_STATUS_IMPORTANT("Gathering Morphology Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}


void Morphology::exportToH5(const std::string& prefix)
{

}

void Morphology::exportToSWC(const std::string& prefix)
{
    reIndexMorphology();

    // Open the file
    std::string fileName = prefix + SWC_EXTENSION;
    std::ofstream stream(fileName.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot write morphology file [ %s ]", fileName.c_str());
    }

    LOG_STATUS("Exporting SWC Morphology : [ %s ]", fileName.c_str());

    // Start the time
    TIMER_SET;

    // Write the vertices
    LOOP_STARTS("Writing Vertices");
    for (size_t i = 0; i < _samples.size(); ++i)
    {
        // LOOP_PROGRESS_FRACTION(i, _samples.size());

        auto& sample = _samples[i];
        stream << sample->getIndex() + 1 << SPACE
               << mapNeuronProcessTypeToSWCIndex(sample->getType()) << SPACE
               << sample->getPosition().x() << SPACE
               << sample->getPosition().y() << SPACE
               << sample->getPosition().z() << SPACE
               << sample->getRadius()       << SPACE
               << sample->getParentIndex() + 1 << NEW_LINE;
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file stream
    stream.close();
}

}
