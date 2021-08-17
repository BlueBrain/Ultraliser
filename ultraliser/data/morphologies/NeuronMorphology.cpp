/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero <juanjose.garcia@epfl.ch>
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

#include "NeuronMorphology.h"
#include <common/Common.h>
#include <data/morphologies/NeuronSWCReader.h>
#include <utilities/Utilities.h>
#include <bits/stdc++.h>
#include <stack>

namespace Ultraliser
{

NeuronMorphology::NeuronMorphology(const NeuronSWCSamples& swcSamples)
{
    // Construct the sections
    _constructMorphology(swcSamples);
}

Samples NeuronMorphology::getSomaSamples() const
{
    return _somaSamples;
}

Sections NeuronMorphology::getFirstSections() const
{
    return _firstSections;
}

Vector3f NeuronMorphology::getSomaCenter() const
{
    return _somaCenter;
}

float NeuronMorphology::getSomaMeanRadius() const
{
    return _somaMeanRadius;
}

float NeuronMorphology::getSomaMinRadius() const
{
    return _somaMinRadius;
}

float NeuronMorphology::getSomaMaxRadius() const
{
    return _somaMaxRadius;
}

void NeuronMorphology::_constructMorphology(const NeuronSWCSamples& swcSamples)
{
    // The _pMin and _pMax will be used to compute the bounding box of the neuron morphology
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples and the sections list from the tree structure
    std::stack<NeuronSWCSample*> samplesStack;
    std::stack<Section*> sectionStack;
    uint64_t sectionId = 0;

    // The first sample is always the SOMA
    samplesStack.push(swcSamples[0]);

    // Initially, nothing is in the section stack
    sectionStack.push(nullptr);

    // Process the samples stack until empty where by then the sections are constructed
    while (!samplesStack.empty())
    {
        // Gets the first element in the samples stack, and then removes it
        auto swcSample = samplesStack.top();
        samplesStack.pop();

        // ?
        auto section = sectionStack.top();
        sectionStack.pop();

        Ultraliser::Vector3f position(swcSample->x, swcSample->y, swcSample->z);
        float radius = swcSample->r;

        Vector3f pMaxSample = position + Vector3f(radius);
        Vector3f pMinSample = position - Vector3f(radius);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();

        // Construct a new sample
        Sample* sample = new Sample(position, radius);

        // If this sample is somatic, either the soma average sample or a profile point
        if (swcSample->type == SWCSampleType::SOMA)
        {
            for (auto child : swcSample->childrenSamples)
            {
                samplesStack.push(child);
                if (child->type == SWCSampleType::SOMA)
                {
                    sectionStack.push(nullptr);
                }
                else
                {
                    auto newSection = new Section(sectionId);
                    _sections.push_back(newSection);
                    _firstSections.push_back(newSection);
                    ++sectionId;
                    sectionStack.push(newSection);
                }
            }
        }

        // Skeleton samples
        else
        {
            // Add the sample to the skeleton samples
            _samples.push_back(sample);

            // Add the sample to the section
            section->addSample(sample);

            // Get the children samples
            auto children = swcSample->childrenSamples;

            // If there is only a single child
            if (children.size() == 1)
            {
                samplesStack.push(children[0]);
                sectionStack.push(section);
            }

            // Multiple children
            else
            {
                for (auto child : children)
                {
                    samplesStack.push(child);
                    section->addChildIndex(sectionId);

                    // Construct a new section corresponding to the child
                    auto newSection = new Section(sectionId);

                    // Add it to the list of sections
                    _sections.push_back(newSection);

                    // Update the counters to keep track
                    ++sectionId;

                    // Parent-child link
                    newSection->addParentIndex(section->getIndex());

                    // Add the sample to the new section
                    newSection->addSample(sample);

                    // Add the new section to the stack
                    sectionStack.push(newSection);
                }
            }
        }
    }

    // Get the somatic samples from the initial sections of all the neurites
    for (const auto& section: _firstSections)
    {
        const auto sample = section->getFirstSample();
        _somaSamples.push_back(sample);
    }

    // Initially, the soma center is set to Zero.
    _somaCenter = Vector3f(0.0f);
    _somaMeanRadius = 0.0f;
    _somaMinRadius = 1e10;
    _somaMaxRadius = -1e10;

    if (_somaSamples.size() > 1)
    {
        // Compute the soma center
        for (auto sample : _somaSamples)
        {
            _somaCenter += sample->getPosition();
        }
        _somaCenter /= _somaSamples.size();

        // Compute the minimum, average and maximum radii
        for (auto sample : _somaSamples)
        {
            const auto radius = (sample->getPosition() - _somaCenter).abs();

            if (radius < _somaMinRadius)
                _somaMinRadius = radius;

            if (radius > _somaMaxRadius)
                _somaMaxRadius = radius;

            _somaMeanRadius += radius;
        }

        // Compute the minimum, average and maximum radii
        _somaMeanRadius /= _somaSamples.size();
    }
    else
    {
        _somaCenter = _somaSamples[0]->getPosition();
        _somaMeanRadius = _somaSamples[0]->getRadius();
        _somaMinRadius = _somaMeanRadius;
        _somaMaxRadius = _somaMeanRadius;
    }
}

NeuronMorphology* readNeuronMorphology(std::string& morphologyPath)
{
    if (String::subStringFound(morphologyPath, std::string(".swc")) ||
        String::subStringFound(morphologyPath, std::string(".SWC")) ||
        String::subStringFound(morphologyPath, std::string(".Swc")))
    {
        // Read the file
        auto reader = std::make_unique< NeuronSWCReader >(morphologyPath);

        // Get a pointer to the morphology to start using it
        return reader->getMorphology();
    }
    else
    {
        LOG_ERROR("Unrecognized morphology file format [ %s ]", morphologyPath.c_str());
    }

    // To avoid any warning issues.
    return nullptr;
}

}  // namespace Ultraliser
