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

bool SkeletonBranch::isSpine() const
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

void SkeletonBranch::resampleAdaptively()
{
    // Get the total number of samples in the section
    const size_t numberSamples = nodes.size();

    // If the section has less than two samples, return as it is not valid
    if (numberSamples < 2)
        return;

    // TODO: If the section has only two samples, it is better not to resample it UFN
    if (numberSamples == 2)
        return;

    /// Initialization
    // Make a copy of the skeleton nodes
    SkeletonNodes newNodes;
    for (size_t i = 0; i < nodes.size(); ++i) { newNodes.push_back(nodes[i]); }

    /// Resample the section
    size_t nodeIndex = isRoot()? 1 : 0;
    // Resample until there is nothing to resample
    while (true)
    {
        // Break if the last sample is reached
        if (nodeIndex >= newNodes.size() - 2) { break; }

        // Compute the distance between the current sample and the next one
        auto& sample0 = newNodes[nodeIndex];
        auto& sample1 = newNodes[nodeIndex + 1];

        const auto distance = sample1->point.distance(sample0->point);

        // Get the extent of the sample, where no other samples should be located
        const auto extent = sample1->radius + sample0->radius;

        // If the next sample is located within the extent of the current sample, then remove it
        if (distance < extent)
        {
            newNodes.erase(newNodes.begin() + nodeIndex + 1);
            continue;
        }

        // If the sample is at a greater step, then add a new sample
        else
        {
            // Compute the radius of the new sample
            const auto radius = (sample0->radius + sample1->radius) * 0.5f;

            // Compute the direction of the next sample
            const auto direction = (sample1->point - sample0->point).normalized();

            // Compute the location of the new sample
            const auto point = sample0->point + (direction * extent);

            // Construct the new sample
            auto newNode = new SkeletonNode(nodeIndex + 1, point, radius);

            // Insert the new node to the list
            newNodes.insert(newNodes.begin() + nodeIndex + 1, newNode);

            // Update the index
            nodeIndex++;

            // If we reach the last sample, then break
            if (nodeIndex >= newNodes.size() - 1) { break; }
        }
    }

    /// Clean-up
    nodes.clear();
    nodes.shrink_to_fit();
    nodes = newNodes;
}

void SkeletonBranch::getBoundingBox(Vector3f& pMin, Vector3f& pMax,
                                    Vector3f& bounds, Vector3f &center)
{
    Vector3f _pMin(FLT_MAX);
    Vector3f _pMax(-1 * FLT_MAX);

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        const auto& point = nodes[i]->point;
        if (point.x() < _pMin.x()) _pMin.x() = point.x();
        if (point.y() < _pMin.y()) _pMin.y() = point.y();
        if (point.z() < _pMin.z()) _pMin.z() = point.z();

        if (point.x() > _pMax.x()) _pMax.x() = point.x();
        if (point.y() > _pMax.y()) _pMax.y() = point.y();
        if (point.z() > _pMax.z()) _pMax.z() = point.z();
    }

    pMin = _pMin;
    pMax = _pMax;
    bounds = _pMax - _pMin;
    center = _pMin + bounds * 0.5f;
}

}
