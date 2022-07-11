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

Section::Section(const uint64_t &index)
    : _index(index)
    , _type(UNKNOWN)
{
    /// EMPTY CONSTRUCTOR
}

SECTION_TYPE Section::getType() const
{
    return _type;
}

void Section::setType(SECTION_TYPE type)
{
    _type = type;
}

uint64_t Section::getIndex() const
{
    return _index;
}

void Section::setIndex(const uint64_t &index)
{
    _index = index;
}

void Section::addParentIndex(const uint64_t index)
{
    _parentsIndices.push_back(index);
}

void Section::clearParentsIndices()
{
    _parentsIndices.clear();
}

void Section::addChildIndex(const uint64_t index)
{
    _childrenIndices.push_back(index);
}

void Section::clearChildrenIndices()
{
    _childrenIndices.clear();
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

std::vector< uint64_t > Section::getParentIndices() const
{
    return _parentsIndices;
}

std::vector< uint64_t > Section::getChildrenIndices() const
{
    return _childrenIndices;
}

float Section::computeLength() const
{
    auto sectionLength = 0.f;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionLength += computeSegmentLength(_samples[i], _samples[i + 1]);
    }
    return sectionLength;
}

std::vector < float > Section:: computeSegmentsLengthDistribution() const
{
    std::vector < float > distribution;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
    {
        distribution.push_back(computeSegmentLength(_samples[i], _samples[i + 1]));
    }
    return distribution;
}

float Section::computeSurfaceArea() const
{
    auto sectionSurfaceArea = 0.f;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionSurfaceArea += computeSegmentSurfaceArea(_samples[i], _samples[i + 1]);
    }
    return sectionSurfaceArea;
}

std::vector < float > Section:: computeSegmentsSurfaceAreaDistribution() const
{
    std::vector < float > distribution;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
    {
        distribution.push_back(computeSegmentSurfaceArea(_samples[i], _samples[i + 1]));
    }
    return distribution;
}

float Section::computeVolume() const
{
    auto sectionVolume = 0.f;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
    {
        sectionVolume += computeSegmentVolume(_samples[i], _samples[i + 1]);
    }
    return sectionVolume;
}

std::vector < float > Section:: computeSegmentsVolumeDistribution() const
{
    std::vector < float > distribution;
    for (uint64_t i = 0; i < _samples.size() - 1; ++i)
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

void Section::resampleUniformly(const float& step)
{
    // Get the total number of samples in the section
    const uint64_t numberSamples = _samples.size();

    // If the section has less than two samples, return as it is not valid
    if (numberSamples < 2)
        return;

    // If the section has two samples only, then verify its length with respect to the given step
    if (numberSamples == 2)
    {
        const auto sample0 = _samples[0];
        const auto sample1 = _samples[1];

        // Compuet the distance between the two samples
        const auto distance = (sample0->getPosition() - sample1->getPosition()).abs();

        // If the distance is less than the step, then we cannot resample this section
        if (distance < step)
            return;

        // Otherwise, the section could be resampled, so RESAMPLE it
        else
        {
            // Create a new samples list
            Samples newSamples;

            // Add the first sample to the section
            newSamples.push_back(_samples[0]);

            // Compute the direction of the section
            Vector3f direction = _samples[1]->getPosition() - _samples[0]->getPosition();
            direction.normalize();

            // This index will keep track on the current sample along the section
            uint64_t index = 1;
            while (true)
            {
                // Compute the position of a new sample
                Vector3f position = newSamples[0]->getPosition() + direction * index * step;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->getPosition()).abs() > distance)
                    break;

                // Otherwise, interpolate the new sample and add it
                const float radius = 0.5 *
                        (newSamples[index - 1]->getRadius() + _samples[1]->getRadius());

                // Add the new sample to the list of new sample
                Sample* sample = new Sample(position, radius, index);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
            }

            // Add the first sample to the section
            newSamples.push_back(_samples.back());

            // Clear the old samples list
            _samples.clear();
            _samples.shrink_to_fit();

            // Update the samples list
            _samples = newSamples;
        }
    }

    // If the section has more than two samples
    else
    {
        // Create a new samples list
        Samples newSamples;

        // Add the first sample to the new list
        newSamples.push_back(_samples[0]);

        // This index will keep track on the current sample along the section
        uint64_t index = 1;

        for (uint64_t i = 0; i < numberSamples - 2; ++i)
        {
            const auto sample0 = _samples[i];
            const auto sample1 = _samples[i + 1];

            // Compute the direction of the section
            Vector3f direction = sample1->getPosition() - sample0->getPosition();
            direction.normalize();

            // Compuet the distance between the two samples
            const auto distance = (sample1->getPosition() - sample0->getPosition()).abs();

            // Proceed wiht the resampling
            uint64_t perSegmentIndex = 1;
            while (true)
            {
                // Compute the position of a new sample
                Vector3f position = sample0->getPosition() + direction * perSegmentIndex * step;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->getPosition()).abs() > distance)
                {
                    // Add the last sample of the segment
                    newSamples.push_back(sample1);
                    break;
                }
                // Otherwise, interpolate the new sample and add it
                const float radius = 0.5 *
                        (newSamples[index - 1]->getRadius() + sample1->getRadius());

                // Add the new sample to the list of new sample
                Sample* sample = new Sample(position, radius, index);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
                ++perSegmentIndex;
            }
        }

        // Add the last sample to ensure to completensee of the section
        newSamples.push_back(_samples.at(_samples.size() - 1));

        // Clear the old samples list
        _samples.clear();
        _samples.shrink_to_fit();

        // Update the samples list
        _samples = newSamples;
    }
}

void Section::resampleAdaptively(const bool& relaxed)
{
    // Get the total number of samples in the section
    const uint64_t numberSamples = _samples.size();

    // If the section has less than two samples, return as it is not valid
    if (numberSamples < 2)
        return;

    // If the sampling is relaxed, then use a distance factor of 2 to account for diameters
    // instead of radii
    auto factor = 1.f;
    if (relaxed)
        factor = 2.f;

    // If the section has two samples only, then verify if it can be resampled or not
    if (numberSamples == 2)
    {
        const auto sample0 = _samples[0];
        const auto sample1 = _samples[1];

        // Compute the distance between samples
        const auto distance = (sample1->getPosition() - sample0->getPosition()).abs();

        // Compute the radii sum
        const auto radii = sample0->getRadius() + sample1->getRadius();

        // If the distance is less than the radii sum, then we cannot resample this section
        if (distance < radii)
            return;

        // Otherwise, the section could be resampled, so RESAMPLE it
        else
        {
            // Create a new samples list
            Samples newSamples;

            // Add the first sample to the section
            newSamples.push_back(_samples[0]);

            // Compute the direction of the section
            Vector3f direction = _samples[1]->getPosition() - _samples[0]->getPosition();
            direction.normalize();

            // This index will keep track on the current sample along the section
            uint64_t index = 1;
            while (true)
            {
                // Get the current sample
                const auto currentSample = newSamples[index - 1];

                // Geth the last sample
                const auto lastSample = _samples[1];

                // Compute the radius of the new sample radius by interpolation with the last sample
                const auto newSampleRadius =
                        0.5f * (currentSample->getRadius() + lastSample->getRadius());

                // Compute the position of the new sample
                const auto position =
                        currentSample->getPosition() + direction * newSampleRadius * factor;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->getPosition()).abs() > distance)
                    break;

                // Add the new sample to the list of new sample
                auto sample = new Sample(position, newSampleRadius, index);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
            }

            // Add the first sample to the section
            newSamples.push_back(_samples.back());

            // Clear the old samples list
            _samples.clear();
            _samples.shrink_to_fit();

            // Update the samples list
            _samples = newSamples;
        }
    }

    // If the section has more than two samples, resample the segments
    else
    {
        // Create a new samples list
        Samples newSamples;

        // Add the first sample to the new list
        newSamples.push_back(_samples[0]);

        // This index will keep track on the current sample along the section
        uint64_t index = 1;

        // On a per-segment basis
        for (uint64_t i = 0; i < numberSamples - 2; ++i)
        {
            const auto sample0 = _samples[i];
            const auto sample1 = _samples[i + 1];

            // Compute the direction of the section
            Vector3f direction = sample1->getPosition() - sample0->getPosition();
            direction.normalize();

            // Compuet the distance between the two samples
            const auto distance = (sample1->getPosition() - sample0->getPosition()).abs();

            // Proceed wiht the resampling
            uint64_t perSegmentIndex = 1;
            while (true)
            {
                // Get the current sample
                const auto currentSample = newSamples[index - 1];

                // Geth the last sample
                const auto lastSample = _samples[1];

                // Compute the radius of the new sample radius by interpolation with the last sample
                const auto newSampleRadius =
                        0.5f * (currentSample->getRadius() + lastSample->getRadius());

                // Compute the position of the new sample
                const auto position =
                        currentSample->getPosition() + direction * newSampleRadius * factor;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->getPosition()).abs() > distance)
                { 
                    // Add the last sample of the segment
                    newSamples.push_back(sample1);
                    break;
                }

                // Add the new sample to the list of new sample
                auto sample = new Sample(position, newSampleRadius, index);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
                ++perSegmentIndex;
            }
        }

        // Clear the old samples list
        _samples.clear();
        _samples.shrink_to_fit();

        // Update the samples list
        _samples = newSamples;
    }
}

}
