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

#include "AstrocyteMorphology.h"
#include <data/morphologies/h5/AstrocyteH5Reader.h>
#include <utilities/String.h>

namespace Ultraliser
{

AstrocyteMorphology::AstrocyteMorphology(Samples& samples, EndfeetPatches& endfeetPatches)
{
    _samples = samples;
    _endfeetPatches = endfeetPatches;

    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());
    for (auto sample: _samples)
    {
        Vector3f pMaxSample = sample->getPosition() + Vector3f(sample->getRadius());
        Vector3f pMinSample = sample->getPosition() - Vector3f(sample->getRadius());

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }    
}


AstrocyteMorphology::AstrocyteMorphology(const H5Samples& h5Samples,
                                         const H5Sections& h5Sections,
                                         const EndfeetPatches &endfeetPatches,
                                         const Vector3f &center)
{
    _endfeetPatches = endfeetPatches;

    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples list that will be used later to build the sections
    for (auto h5Sample : h5Samples)
    {
        Ultraliser::Vector3f position(h5Sample.x, h5Sample.y, h5Sample.z);
        Sample* sample = new Sample(
            Ultraliser::Vector3f(h5Sample.x, h5Sample.y, h5Sample.z),
                                    h5Sample.r * 0.5f, 0);
        _samples.push_back(sample);

        Vector3f pMaxSample = position + Vector3f(h5Sample.r * 0.5f);
        Vector3f pMinSample = position - Vector3f(h5Sample.r * 0.5f);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }

    // Create the soma section, for indexing only
    Section* somaSection = new Section(0);
    _sections.push_back(somaSection);

    // The soma is located at index 0, but we will consider it a section for indexing!
    for (uint64_t i = 1; i < h5Sections.size() - 1; ++i)
    {
        // Current section index
        const auto sectionIndex = i;

        // First point index
        const auto firstPointIdx = h5Sections[i].offsetIndex;

        // Last point index
        const auto lastPointIdx = h5Sections[i + 1].offsetIndex - 1;

        // Create a new morphology section
        Section* section = new Section(sectionIndex);

        // Construct a list of samples
        Samples samples;
        for (uint64_t s = firstPointIdx; s <= lastPointIdx; s++)
        {
            // Create the sample
            Sample* sample =
                    new Sample(Ultraliser::Vector3f(h5Samples[s].x, h5Samples[s].y, h5Samples[s].z),
                               h5Samples[s].r * 0.5f, 0);

            // Update the samples list
            section->addSample(sample);
        }

        // Update the section parent
        section->addParentIndex(h5Sections[i].parentIndex);

        // Section type
        const auto sectionType = h5Sections[i].sectionType;

        switch (h5Sections[i].sectionType)
        {
        case 2:
        {
            section->setType(NEURON_AXON);
        }
        case 3:
        {
            section->setType(NEURON_BASAL_DENDRITE);
        }
        case 4:
        {
            section->setType(NEURON_APICAL_DENDRITE);
        }
        default:
        {
            section->setType(NEURON_BASAL_DENDRITE);
        }
        }

        // Add the section to the list
        _sections.push_back(section);

        // Determine if this is a first section or not!
        if (h5Sections[i].parentIndex == 0)
        {
            _firstSections.push_back(section);
        }
    }

    // Build the tree (add the children indices) from the linear list
    for (uint64_t i = 1; i < _sections.size(); ++i)
    {
        // Get parent index
        uint64_t parentSectionIndex = _sections[i]->getParentIndices()[0];

        // Update the choldren list
        _sections[parentSectionIndex]->addChildIndex(_sections[i]->getIndex());
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

Samples AstrocyteMorphology::getSomaSamples() const
{
    return _somaSamples;
}

Sections AstrocyteMorphology::getFirstSections() const
{
    return _firstSections;
}

Vector3f AstrocyteMorphology::getSomaCenter() const
{
    return _somaCenter;
}

float AstrocyteMorphology::getSomaMeanRadius() const
{
    return _somaMeanRadius;
}

float AstrocyteMorphology::getSomaMinRadius() const
{
    return _somaMinRadius;
}

float AstrocyteMorphology::getSomaMaxRadius() const
{
    return _somaMaxRadius;
}

EndfeetPatches AstrocyteMorphology::getEndfeetPatches() const
{
    return _endfeetPatches;
}

AstrocyteMorphology* readAstrocyteMorphology(std::string& morphologyPath)
{
     if (String::subStringFound(morphologyPath, std::string(".h5")) ||
         String::subStringFound(morphologyPath, std::string(".H5")))
     {
         // Read the file
         auto reader = std::make_unique< AstrocyteH5Reader >(morphologyPath);

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

}
