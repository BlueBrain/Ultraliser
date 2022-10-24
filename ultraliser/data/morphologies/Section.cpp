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

#include "Section.h"
#include "MorphologyOpertions.h"

namespace Ultraliser
{

Section::Section(const size_t &index, const PROCESS_TYPE &type)
    : _index(index)
    , _type(type)
{
    /// EMPTY CONSTRUCTOR
}

PROCESS_TYPE Section::getType() const
{
    return _type;
}

size_t Section::getIndex() const
{
    return _index;
}

void Section::setIndex(const size_t &index)
{
    _index = index;
}

void Section::addParentIndex(const size_t index)
{
    _parentsIndices.push_back(index);
}

void Section::addParent(Section* section)
{
    _parents.push_back(section);
}

void Section::addChild(Section* section)
{
    _children.push_back(section);
}

void Section::clearParentsIndices()
{
    _parentsIndices.clear();
}

void Section::addChildIndex(const size_t index)
{
    _childrenIndices.push_back(index);
}

void Section::clearChildrenIndices()
{
    _childrenIndices.clear();
}

bool Section::isRoot()
{
    if (_parentsIndices.size() == 0)
    {
        return true;
    }
    return false;
}

size_t Section::getBranchingOrder() const
{
    return _branchingOrder;
}

void Section::addSample(Sample* sample)
{
    _samples.push_back(sample);
}

void Section::addSamples(const Samples& samples)
{
    _samples.insert(_samples.end(), samples.begin(), samples.end());
}

Samples Section::getSamples() const
{
    return _samples;
}

Sample* Section::getFirstSample() const
{
    return _samples.front();
}

Sample* Section::getLastSample() const
{
    return _samples.back();
}

std::vector< size_t > Section::getParentIndices() const
{
    return _parentsIndices;
}

std::vector< size_t > Section::getChildrenIndices() const
{
    return _childrenIndices;
}

float Section::computeLength() const
{
    auto sectionLength = 0.f;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionLength += computeSegmentLength(_samples[i], _samples[i + 1]);
    }
    return sectionLength;
}

std::vector < float > Section:: computeSegmentsLengthDistribution() const
{
    std::vector < float > distribution;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        distribution.push_back(computeSegmentLength(_samples[i], _samples[i + 1]));
    }
    return distribution;
}

float Section::computeSurfaceArea() const
{
    auto sectionSurfaceArea = 0.f;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionSurfaceArea += computeSegmentSurfaceArea(_samples[i], _samples[i + 1]);
    }
    return sectionSurfaceArea;
}

std::vector < float > Section:: computeSegmentsSurfaceAreaDistribution() const
{
    std::vector < float > distribution;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        distribution.push_back(computeSegmentSurfaceArea(_samples[i], _samples[i + 1]));
    }
    return distribution;
}

float Section::computeVolume() const
{
    auto sectionVolume = 0.f;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionVolume += computeSegmentVolume(_samples[i], _samples[i + 1]);
    }
    return sectionVolume;
}

std::vector < float > Section:: computeSegmentsVolumeDistribution() const
{
    std::vector < float > distribution;
    for (size_t i = 0; i < _samples.size() - 1; ++i)
    {
        distribution.push_back(computeSegmentVolume(_samples[i], _samples[i + 1]));
    }
    return distribution;
}

float Section::computeAverageRadius() const
{
    auto sectionAverageRadius = 0.f;
    for (const auto sample: _samples)
    {
        sectionAverageRadius += sample->getRadius();
    }
    return sectionAverageRadius / _samples.size();
}

size_t Section::computeNumberZeroRadiusSamples(const float& threshold) const
{
    size_t count = 0;
    for (const auto sample: _samples)
    {
        if (sample->getRadius() < threshold)
        {
            count++;
        }
    }
    return count;
}

void Section::updateSamplesIndices()
{
    // Update the indices of the samples, locally
    for (size_t i = 0; i < _samples.size(); ++i)
        _samples[i]->setIndex(i);
}

void Section::resampleSegment(const float& step, Samples& samples)
{
    // The segment can be resampled, if it only has two samples, otherwise return
    if (samples.size() == 2)
    {
        const auto sample0 = samples[0];
        const auto sample1 = samples[1];

        // Compuet the distance between the two samples, i.e. the segmentLength
        const auto segmentLength = (sample0->getPosition() - sample1->getPosition()).abs();

        // If the step is greater than the segmentLength, then we cannot resample this segment
        if (step > segmentLength)
            return;

        // Compute the direction of the segment
        Vector3f direction = sample1->getPosition() - sample0->getPosition();
        direction.normalize();

        // Sampling index
        size_t index = 1;
        while (true)
        {
            // Compute the position of a new sample
            const auto position = sample0->getPosition() + direction * (index) * step;

            // Compute the delta
            const auto delta = (position - sample0->getPosition()).abs();

            // If we exceed the segmentLength, then break
            if (delta > segmentLength)
            {
                break;
            }

            // Compute the radius
            const auto radius = 0.5 *
                    (samples[index - 1]->getRadius() + samples[index]->getRadius());

            // Construct the new sample
            Sample* sample = new Sample(position, radius, sample0->getType(),
                    AUXILIARY_SAMPLE_INDEX, AUXILIARY_SAMPLE_INDEX);
            samples.insert(samples.begin() + index, sample);
            index++;
        }
    }
}

void Section::resampleUniformly(const float& step)
{
    // Get the total number of samples in the section
    const auto numberSamples = _samples.size();

    // If the section has less than two samples, return. The section cannot be resampled.
    if (numberSamples < 2)
        return;

    // If the section has two samples only, then verify its length with respect to the given step
    if (numberSamples == 2)
    {
        resampleSegment(step, _samples);
    }
    else
    {
        //
        Samples sectionSamples;
        sectionSamples.push_back(_samples[0]);
        for (size_t i = 0; i < _samples.size() - 1; ++i)
        {
            Samples segmentSamples;
            segmentSamples.resize(2);
            segmentSamples[0] = _samples[i];
            segmentSamples[1] = _samples[i + 1];
            resampleSegment(step, segmentSamples);

            for (size_t j = 1; j < segmentSamples.size(); ++j)
            {
                sectionSamples.push_back(segmentSamples[j]);
            }
        }

        // Clear the old samples list
        _samples.clear();
        _samples.shrink_to_fit();

        // Update the samples list
        _samples = sectionSamples;
    }

    // Update the indices of the samples
    updateSamplesIndices();
}

size_t Section::removeIntermediateSamples()
{
    // Get the total number of samples in the section
    const size_t numberOriginalSamples = _samples.size();

    // If the section has less than two samples, return
    if (numberOriginalSamples < 2)
        return 0;

    // If the section has two samples, return. The section cannot be resampled.
    if (numberOriginalSamples == 2)
        return 0;

    Samples newSamples;
    newSamples.push_back(new Sample(_samples.front()));
    newSamples.push_back(new Sample(_samples.back()));

    for (auto& sample: _samples) { delete sample; }
    _samples.clear();
    _samples.shrink_to_fit();

    _samples = newSamples;

    // Keep track on the removed samples
    size_t numberRemovedSamples = numberOriginalSamples - _samples.size();

    // Return the total number of removed samples.
    return numberRemovedSamples;
}

void Section::removeOverlappingSamples()
{
    if (_samples.size() <=2)
        return;

    size_t index = 0;
    while(true)
    {
        if (index >= _samples.size() - 1)
            break;

        auto sample0 = _samples[index];
        auto sample1 = _samples[index + 1];

        const auto r0 = sample0->getRadius();
        const auto r1 = sample1->getRadius();

        const auto p0 = sample0->getPosition();
        const auto p1 = sample1->getPosition();

        const auto delta = (p1 - p0).abs();
        const auto radiiSum = r0 + r1;

        // If they are overlapping, then remove the second sample
        if (delta < radiiSum)
        {
            // The second sample must not be a terminal sample
            if (index >= _samples.size() - 1)
                break;

            // Delete the second sample
            _samples.erase(_samples.begin() + index + 1);

            // Reset
            index = 0;
        }

        // Next sample
        index++;
    }
}

size_t Section::removeInnerSamples()
{
    // Get the total number of samples in the section
    const size_t numberSamples = _samples.size();

    // If the section has less than two samples, return. The section cannot be resampled.
    if (numberSamples < 2)
        return 0;

    // If the section has two samples, return. The section cannot be resampled.
    if (numberSamples == 2)
        return 0;

    // Keep track on the removed samples
    size_t numberRemovedSamples = 0;

    // Keep track on current sample, its index
    size_t i = 0;
    while (true)
    {
        // Compute the distance between samples
        const auto distance = (_samples[i + 1]->getPosition() - _samples[i]->getPosition()).abs();

        // Compute the radii sum
        const auto& r0 = _samples[i]->getRadius();
        const auto& r1 = _samples[i + 1]->getRadius();

        // If sample0 is located entirely within sample1
        if (distance < r1)
        {
            // If sample0 is the first sample in the section
            if (i == 0)
            {
                // Use sample1 as the first sample of the section, and remove sample0
                _samples[i + 1]->setPosition(_samples[i]->getPosition());
                _samples.erase(_samples.begin());
            }
            // Otherwise, remove sample0 from the section
            else
            {
                // Remove sample0 @i
                _samples.erase(_samples.begin() + i);
            }

            // Increment the counter
            numberRemovedSamples++;
        }

        // If sample1 is located entirely within sample0
        else if (distance < r0)
        {
            // If sample1 is the last sample in the section
            if (i == _samples.size() - 1)
            {
                // Use sample0 as the last sample of the section and remove sample1
                _samples[i]->setPosition(_samples[i + 1]->getPosition());
                _samples.erase(_samples.end());
            }
            // Otherwise, remove sample1 from the section
            else
            {
                // Remove sample1 @i+1
                _samples.erase(_samples.begin() + i + 1);
            }

            // Increment the counter
            numberRemovedSamples++;
        }
        else
        {
            // Next sample
            i++;
        }

        // We have reached the terminal
        if (i == _samples.size() - 2)
            break;
    }

    // Return the number of removed samples
    return numberRemovedSamples;
}

size_t Section::interpolateLongSegments()
{
    // Get the total number of samples in the section
    const size_t numberSamples = _samples.size();

    // If the section has less than two samples, return
    if (numberSamples < 2)
        return 0;

    // Keep track on the total number of samples added
    size_t numberAddedSamples = 0;

    // Keep track on current sample, its index
    size_t i = 0;
    while (true)
    {
        auto sample0 = _samples[i];
        auto sample1 = _samples[i + 1];

        // Compute the distance between samples
        const auto distance = (sample1->getPosition() - sample0->getPosition()).abs();

        // If the distance between the two samples is greater than the radii, insert a sample
        if (distance > sample0->getRadius() + sample1->getRadius())
        {
            // Compute the direction of the segment from sample0 to sample1
            const auto direction = (sample1->getPosition() - sample0->getPosition()).normalized();

            // Compute the interpolated radius
            const auto radius = (sample0->getRadius() + sample1->getRadius()) * 0.5;

            // Compute the new position of the sample
            const auto position =
                    sample0->getPosition() + (direction * radius) - Vector3f(FLT_EPSILON);

            // Create the new sample, for the moment use the index 0 until all the auxillary
            // samples are created
            Sample* interpolatedSample = new Sample(position, radius, sample0->getType(), 0);

            // Add the new sample
            _samples.insert(_samples.begin() + i + 1, interpolatedSample);

            // Increment the counter
            numberAddedSamples++;
        }

        // Next sample
        i++;

        if (i >= _samples.size() - 2)
            break;
    }

    // Return the number of added samples
    return numberAddedSamples;
}

void Section::resampleAdaptively(const bool& relaxed)
{
    // Get the total number of samples in the section
    const size_t numberSamples = _samples.size();

    // If the section has less than two samples, return as it is not valid
    if (numberSamples < 2)
        return;

    // Get the smallest radius
    const auto smallestRadius = getMinimumSampleRadius();

    // Resample uniformly based on the radius of the smallest sample
    resampleUniformly(smallestRadius);

    // Remove the ovrlapping samples
    removeOverlappingSamples();
}

float Section::getMinimumSampleRadius() const
{
    // Initially, set it to maximum float value
    auto minRadius = std::numeric_limits< float >::max();

    // Compare with all the section samples
    for (auto& sample: _samples)
    {
        // If less, then update the minRadius value
        if (sample->getRadius() < minRadius)
        {
            minRadius = sample->getRadius();
        }
    }

    // Returns the minimum radius
    return minRadius;
}

void Section::verifyMinimumSampleRadius(const float& radius)
{
    for (auto& sample: _samples)
    {
        // If the sample radius is smaller than the threshold, update the sample radius
        if (sample->getRadius() < radius)
        {
            sample->setRadius(radius);
        }
    }
}

void Section::compileSWCTableRecursively(Samples& samples,
                                         size_t &currentSampleIndex,
                                         const size_t& branchingSampleIndex)
{
    /// NOTE: Section samples are _samples, but the collecting list is samples
    if (isRoot())
    {
        // The parentIndex of the first sample is 1 if the section isRoot()
        samples.push_back(new Sample(_samples[0]->getPosition(),
                                     _samples[0]->getRadius(),
                                     _samples[0]->getType(),
                                     currentSampleIndex++, 1));

        // Start directly at the second sample to complete the section
        for (size_t i = 1; i < _samples.size(); ++i)
        {
            samples.push_back(new Sample(_samples[i]->getPosition(),
                                         _samples[i]->getRadius(),
                                         _samples[i]->getType(),
                                         currentSampleIndex,
                                         currentSampleIndex - 1));
            currentSampleIndex++;
        }
    }
    else
    {
        // The parentIndex of the sample is the branchingSampleIndex if not isRoot()
        samples.push_back(new Sample(_samples[1]->getPosition(),
                                     _samples[1]->getRadius(),
                                     _samples[1]->getType(),
                                     currentSampleIndex++, branchingSampleIndex));

        // Start at the third sample
        for (size_t i = 2; i < _samples.size(); ++i)
        {
            samples.push_back(new Sample(_samples[i]->getPosition(),
                                         _samples[i]->getRadius(),
                                         _samples[i]->getType(),
                                         currentSampleIndex,
                                         currentSampleIndex - 1));
            currentSampleIndex++;
        }
    }

    // The branchingSampleIndex will be the same for all the children sections
    const size_t childBranchingSampleIndex = samples[samples.size() - 1]->getIndex();

    // Re-index the children in a recursive manner
    for (size_t i = 0; i < _children.size(); ++i)
    {
        // Re-index the children sections
        _children[i]->compileSWCTableRecursively(samples,
                                                 currentSampleIndex,
                                                 childBranchingSampleIndex);
    }
}
}
