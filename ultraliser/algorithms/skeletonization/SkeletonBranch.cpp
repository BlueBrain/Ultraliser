/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#define VALIDAITY_BIT_INDEX     0
#define VISITED_BIT_INDEX       1
#define TERMINAL_BIT_INDEX      2
#define ROOT_BIT_INDEX          3
#define DUPLICATE_BIT_INDEX     4
#define SPINE_BIT_INDEX         5
#define INSIDE_SOMA_BIT_INDEX   6

#include "SkeletonBranch.h"

namespace Ultraliser
{

SkeletonBranch::SkeletonBranch()
{
    _flags = new BitArray(8);
    _flags->clearAll();
}

SkeletonBranch::~SkeletonBranch()
{
    _flags->~BitArray();
}

void SkeletonBranch::printTree(size_t order)
{
    std::cout << std::string(order * 4, '-') << index << "\n";
    for (size_t i = 0; i < children.size(); ++i)
    {
        children[i]->printTree(order + 1);
    }
}

bool SkeletonBranch::hasTerminalNodes(const size_t& node1Index, const size_t& node2Index)
{
    const auto& frontNode = nodes.front();
    const auto& backNode = nodes.back();

    if (node1Index == backNode->index && node2Index == frontNode->index)
    {
        return true;
    }

    if (node1Index == frontNode->index && node2Index == backNode->index)
    {
        return true;
    }

    return false;
}

void SkeletonBranch::adjustDirection(const size_t& frontNodeIndex, const size_t& backNodeIndex)
{
    if (nodes.front()->index == backNodeIndex && nodes.back()->index == frontNodeIndex)
    {
        std::reverse(this->nodes.begin(), this->nodes.end());
    }
}

bool SkeletonBranch::isLoop() const
{
    return this->nodes.front()->index == this->nodes.back()->index;
}

void SkeletonBranch::setValid()
{
    _flags->setBit(VALIDAITY_BIT_INDEX);
}

void SkeletonBranch::setInvalid()
{
    _flags->clearBit(VALIDAITY_BIT_INDEX);
}

bool SkeletonBranch::isValid() const
{
    return _flags->bit(VALIDAITY_BIT_INDEX);
}

void SkeletonBranch::setVisited()
{
    _flags->setBit(VISITED_BIT_INDEX);
}

void SkeletonBranch::setUnvisited()
{
     _flags->clearBit(VISITED_BIT_INDEX);
}

bool SkeletonBranch::visited() const
{
    return _flags->bit(VISITED_BIT_INDEX);
}

void SkeletonBranch::setTerminal()
{
    _flags->setBit(TERMINAL_BIT_INDEX);
}

void SkeletonBranch::unsetTermainal()
{
    _flags->clearBit(TERMINAL_BIT_INDEX);
}

bool SkeletonBranch::isTerminal() const
{
    return _flags->bit(TERMINAL_BIT_INDEX);
}

void SkeletonBranch::setRoot()
{
    _flags->setBit(ROOT_BIT_INDEX);
}

void SkeletonBranch::unsetRoot()
{
    _flags->clearBit(ROOT_BIT_INDEX);
}

bool SkeletonBranch::isRoot() const
{
    return _flags->bit(ROOT_BIT_INDEX);
}

void SkeletonBranch::setDuplicate()
{
    _flags->setBit(DUPLICATE_BIT_INDEX);
}

void SkeletonBranch::unsetDuplicate()
{
    _flags->clearBit(DUPLICATE_BIT_INDEX);
}

bool SkeletonBranch::isDuplicate() const
{
    return _flags->bit(DUPLICATE_BIT_INDEX);
}

float SkeletonBranch::getMinimumRadius() const
{
    float minRadius = 1e32;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i]->radius < minRadius)
        {
            minRadius = nodes[i]->radius;
        }
    }
    return minRadius;
}

float SkeletonBranch::getMaximumRadius() const
{
    float maxRadius = 0.f;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        if (nodes[i]->radius > maxRadius)
        {
            maxRadius = nodes[i]->radius;
        }
    }
    return maxRadius;
}

float SkeletonBranch::computeLength() const
{
    float length = 0;
    for (size_t i = 0; i < nodes.size() - 1; ++i)
    {
        auto p1 = nodes[i]->point;
        auto p2 = nodes[i + 1]->point;
        length += p1.distance(p2);
    }

    return length;
}

float SkeletonBranch::computeAverageRadius() const
{
    float averageRadius = 0;
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        averageRadius += nodes[i]->radius;
    }
    averageRadius /= nodes.size();

    return averageRadius;
}

void SkeletonBranch::setSpine()
{
    _flags->setBit(SPINE_BIT_INDEX);
}

bool SkeletonBranch::isSpine()
{
    return _flags->bit(SPINE_BIT_INDEX);
}

void SkeletonBranch::setInsideSoma()
{
    _flags->setBit(INSIDE_SOMA_BIT_INDEX);
}

bool SkeletonBranch::isInsideSoma()
{
    return _flags->bit(INSIDE_SOMA_BIT_INDEX);
}

void SkeletonBranch::resampleAdaptively(const bool& relaxed)
{
    // Get the total number of samples in the section
    const size_t numberSamples = nodes.size();

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
        const auto sample0 = nodes[0];
        const auto sample1 = nodes[1];

        // Compute the distance between samples
        const auto distance = (sample1->point - sample0->point).abs();

        // Compute the radii sum
        const auto radii = sample0->radius + sample1->radius;

        // If the distance is less than the radii sum, then we cannot resample this section
        if (distance < radii)
            return;

        // Otherwise, the section could be resampled, so RESAMPLE it
        else
        {
            // Create a new samples list
            SkeletonNodes newSamples;

            // Add the first sample to the section
            newSamples.push_back(nodes[0]);

            // Compute the direction of the section
            Vector3f direction = nodes[1]->point - nodes[0]->point;
            direction.normalize();

            // This index will keep track on the current sample along the section
            size_t index = 1;
            while (true)
            {
                // Get the current sample
                const auto currentSample = newSamples[index - 1];

                // Geth the last sample
                const auto lastSample = nodes[1];

                // Compute the radius of the new sample radius by interpolation with the last sample
                const auto newSampleRadius =
                        0.5f * (currentSample->radius + lastSample->radius);

                // Compute the position of the new sample
                const auto position =
                        currentSample->point + direction * newSampleRadius * factor;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->point).abs() > distance)
                    break;

                // Add the new sample to the list of new sample
                auto sample = new SkeletonNode(index, position, newSampleRadius);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
            }

            // Add the first sample to the section
            newSamples.push_back(nodes.back());

            // Clear the old samples list
            nodes.clear();
            nodes.shrink_to_fit();

            // Update the samples list
            nodes = newSamples;
        }
    }

    // If the section has more than two samples, resample the segments
    else
    {
        // Create a new samples list
        SkeletonNodes newSamples;

        // Add the first sample to the new list
        newSamples.push_back(nodes[0]);

        // This index will keep track on the current sample along the section
        size_t index = 1;

        // On a per-segment basis
        for (size_t i = 0; i < numberSamples - 2; ++i)
        {
            const auto sample0 = nodes[i];
            const auto sample1 = nodes[i + 1];

            // Compute the direction of the section
            Vector3f direction = sample1->point - sample0->point;
            direction.normalize();

            // Compuet the distance between the two samples
            const auto distance = (sample1->point - sample0->point).abs();

            // Proceed wiht the resampling
            size_t perSegmentIndex = 1;
            while (true)
            {
                // Get the current sample
                const auto currentSample = newSamples[index - 1];

                // Geth the last sample
                const auto lastSample = nodes[1];

                // Compute the radius of the new sample radius by interpolation with the last sample
                const auto newSampleRadius =
                        0.5f * (currentSample->radius + lastSample->radius);

                // Compute the position of the new sample
                const auto position = currentSample->point + direction * newSampleRadius * factor;

                // If the position of the new sample goes beyond the segment length, break
                if ((position - sample0->point).abs() > distance)
                {
                    // Add the last sample of the segment
                    newSamples.push_back(sample1);
                    break;
                }

                // Add the new sample to the list of new sample
                auto sample = new SkeletonNode(index, position, newSampleRadius);
                newSamples.push_back(sample);

                // Increase the sample index
                ++index;
                ++perSegmentIndex;
            }
        }

        // Clear the old samples list
        nodes.clear();
        nodes.shrink_to_fit();

        // Update the samples list
        nodes = newSamples;
    }
}

}
