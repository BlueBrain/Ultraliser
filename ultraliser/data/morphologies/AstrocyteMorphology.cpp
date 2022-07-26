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
#include <data/morphologies/swc/NeuronSWCReader.h>
#include <utilities/String.h>
#include <stack>

namespace Ultraliser
{

AstrocyteMorphology::AstrocyteMorphology(const NeuronSWCSamples& swcSamples)
{
    // Construct the morphology from the loaded data
    _constructMorphologyFromSWC(swcSamples);
}

AstrocyteMorphology::AstrocyteMorphology(Samples& samples, EndfeetPatches& endfeetPatches)
{
    _samples = samples;
    _endfeetPatches = endfeetPatches;

    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    _smallestRadiusInMorphology = std::numeric_limits<float>::max();
    _largestRadiusInMorphology = std::numeric_limits<float>::lowest();

    for (auto sample: _samples)
    {
        const float radius = sample->getRadius();
        if (radius  < _smallestRadiusInMorphology && radius > 0)
            _smallestRadiusInMorphology = radius ;

        if (radius  > _largestRadiusInMorphology)
            _largestRadiusInMorphology = radius;

        Vector3f pMaxSample = sample->getPosition() + Vector3f(sample->getRadius());
        Vector3f pMinSample = sample->getPosition() - Vector3f(sample->getRadius());

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }    

    // Extend the bounding box based on endfeet patches as well
    for (size_t i = 0; i < _endfeetPatches.size(); ++i)
    {
        const EndfootPatches patches = _endfeetPatches[i];

        for (size_t j = 0; j < patches.size(); ++j)
        {
            const EndfootPatch* patch = patches[j];

            if (patch->sample0->getRadius() < _smallestRadiusInMorphology)
                _smallestRadiusInMorphology = patch->sample0->getRadius();

            if (patch->sample0->getRadius() > _largestRadiusInMorphology)
                _largestRadiusInMorphology = patch->sample0->getRadius();

            // First sample
            const Vector3f sample0Max =
                    patch->sample0->getPosition() + Vector3f(patch->sample0->getRadius()* 0.5f);
            const Vector3f sample0Min =
                    patch->sample0->getPosition() - Vector3f(patch->sample0->getRadius()* 0.5f);

            if (sample0Max.x() > _pMax.x()) _pMax.x() = sample0Max.x();
            if (sample0Max.y() > _pMax.y()) _pMax.y() = sample0Max.y();
            if (sample0Max.z() > _pMax.z()) _pMax.z() = sample0Max.z();

            if (sample0Min.x() < _pMin.x()) _pMin.x() = sample0Min.x();
            if (sample0Min.y() < _pMin.y()) _pMin.y() = sample0Min.y();
            if (sample0Min.z() < _pMin.z()) _pMin.z() = sample0Min.z();

            // Second sample
            const Vector3f sample1Max =
                    patch->sample1->getPosition() + Vector3f(patch->sample1->getRadius()* 0.5f);
            const Vector3f sample1Min =
                    patch->sample1->getPosition() - Vector3f(patch->sample1->getRadius()* 0.5f);

            if (sample1Max.x() > _pMax.x()) _pMax.x() = sample1Max.x();
            if (sample1Max.y() > _pMax.y()) _pMax.y() = sample1Max.y();
            if (sample1Max.z() > _pMax.z()) _pMax.z() = sample1Max.z();

            if (sample1Min.x() < _pMin.x()) _pMin.x() = sample1Min.x();
            if (sample1Min.y() < _pMin.y()) _pMin.y() = sample1Min.y();
            if (sample1Min.z() < _pMin.z()) _pMin.z() = sample1Min.z();

            // Third sample
            const Vector3f sample2Max =
                    patch->sample2->getPosition() + Vector3f(patch->sample2->getRadius()* 0.5f);
            const Vector3f sample2Min =
                    patch->sample2->getPosition() - Vector3f(patch->sample2->getRadius()* 0.5f);

            if (sample2Max.x() > _pMax.x()) _pMax.x() = sample2Max.x();
            if (sample2Max.y() > _pMax.y()) _pMax.y() = sample2Max.y();
            if (sample2Max.z() > _pMax.z()) _pMax.z() = sample2Max.z();

            if (sample2Min.x() < _pMin.x()) _pMin.x() = sample2Min.x();
            if (sample2Min.y() < _pMin.y()) _pMin.y() = sample2Min.y();
            if (sample2Min.z() < _pMin.z()) _pMin.z() = sample2Min.z();
        }
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

    _smallestRadiusInMorphology = std::numeric_limits<float>::max();
    _largestRadiusInMorphology = std::numeric_limits<float>::lowest();

    // Construct the samples list that will be used later to build the sections
    // Compute the bounding box as well
    for (auto h5Sample : h5Samples)
    {
        Ultraliser::Vector3f position(h5Sample.x, h5Sample.y, h5Sample.z);
        Sample* sample = new Sample(
            Ultraliser::Vector3f(h5Sample.x, h5Sample.y, h5Sample.z),
                                    h5Sample.r * 0.5f, 0);
        _samples.push_back(sample);

        if (h5Sample.r < _smallestRadiusInMorphology && h5Sample.r > 0)
            _smallestRadiusInMorphology = h5Sample.r;

        if (h5Sample.r > _largestRadiusInMorphology)
            _largestRadiusInMorphology = h5Sample.r;

        Vector3f pMaxSample = position + Vector3f(h5Sample.r * 0.5f);
        Vector3f pMinSample = position - Vector3f(h5Sample.r * 0.5f);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }

    // Extend the bounding box based on endfeet patches as well
    for (size_t i = 0; i < _endfeetPatches.size(); ++i)
    {
        const EndfootPatches patches = _endfeetPatches[i];

        for (size_t j = 0; j < patches.size(); ++j)
        {
            const EndfootPatch* patch = patches[j];

            if (patch->sample0->getRadius() < _smallestRadiusInMorphology)
                _smallestRadiusInMorphology = patch->sample0->getRadius();

            if (patch->sample0->getRadius() > _largestRadiusInMorphology)
                _largestRadiusInMorphology = patch->sample0->getRadius();

            // First sample
            const Vector3f sample0Max =
                    patch->sample0->getPosition() + Vector3f(patch->sample0->getRadius()* 0.5f);
            const Vector3f sample0Min =
                    patch->sample0->getPosition() - Vector3f(patch->sample0->getRadius()* 0.5f);

            if (sample0Max.x() > _pMax.x()) _pMax.x() = sample0Max.x();
            if (sample0Max.y() > _pMax.y()) _pMax.y() = sample0Max.y();
            if (sample0Max.z() > _pMax.z()) _pMax.z() = sample0Max.z();

            if (sample0Min.x() < _pMin.x()) _pMin.x() = sample0Min.x();
            if (sample0Min.y() < _pMin.y()) _pMin.y() = sample0Min.y();
            if (sample0Min.z() < _pMin.z()) _pMin.z() = sample0Min.z();

            // Second sample
            const Vector3f sample1Max =
                    patch->sample1->getPosition() + Vector3f(patch->sample1->getRadius()* 0.5f);
            const Vector3f sample1Min =
                    patch->sample1->getPosition() - Vector3f(patch->sample1->getRadius()* 0.5f);

            if (sample1Max.x() > _pMax.x()) _pMax.x() = sample1Max.x();
            if (sample1Max.y() > _pMax.y()) _pMax.y() = sample1Max.y();
            if (sample1Max.z() > _pMax.z()) _pMax.z() = sample1Max.z();

            if (sample1Min.x() < _pMin.x()) _pMin.x() = sample1Min.x();
            if (sample1Min.y() < _pMin.y()) _pMin.y() = sample1Min.y();
            if (sample1Min.z() < _pMin.z()) _pMin.z() = sample1Min.z();

            // Third sample
            const Vector3f sample2Max =
                    patch->sample2->getPosition() + Vector3f(patch->sample2->getRadius()* 0.5f);
            const Vector3f sample2Min =
                    patch->sample2->getPosition() - Vector3f(patch->sample2->getRadius()* 0.5f);

            if (sample2Max.x() > _pMax.x()) _pMax.x() = sample2Max.x();
            if (sample2Max.y() > _pMax.y()) _pMax.y() = sample2Max.y();
            if (sample2Max.z() > _pMax.z()) _pMax.z() = sample2Max.z();

            if (sample2Min.x() < _pMin.x()) _pMin.x() = sample2Min.x();
            if (sample2Min.y() < _pMin.y()) _pMin.y() = sample2Min.y();
            if (sample2Min.z() < _pMin.z()) _pMin.z() = sample2Min.z();
        }
    }

    // Create the soma section, for indexing only
    Section* somaSection = new Section(0);
    _sections.push_back(somaSection);

    // The soma is located at index 0, but we will consider it a section for indexing!
    for (size_t i = 1; i < h5Sections.size() - 1; ++i)
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
        for (size_t s = firstPointIdx; s <= lastPointIdx; s++)
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
    for (size_t i = 1; i < _sections.size(); ++i)
    {
        // Get parent index
        size_t parentSectionIndex = _sections[i]->getParentIndices()[0];

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

float AstrocyteMorphology::getSmallestRadiusInMorphology() const
{
    return _smallestRadiusInMorphology;
}

float AstrocyteMorphology::getLargestRadiusInMorphology() const
{
    return _largestRadiusInMorphology;
}

void AstrocyteMorphology::_constructMorphologyFromSWC(const NeuronSWCSamples& swcSamples)
{
    // The _pMin and _pMax will be used to compute the bounding box of the neuron morphology
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples and the sections list from the tree structure
    std::stack< NeuronSWCSample* > samplesStack;
    std::stack< Section* > sectionStack;
    size_t sectionId = 0;

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
        Sample* sample = new Sample(position, radius, swcSample->id);

        // If this sample is somatic, either the soma average sample or a profile point
        if (swcSample->type == SWCSampleType::SOMA)
        {
            // _somaSamples.push_back(sample);
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
                    switch (child->type)
                    {
                    case AXON:
                        newSection->setType(NEURON_AXON);
                        break;
                    case BASAL:
                        newSection->setType(NEURON_BASAL_DENDRITE);
                        break;
                    case APICAL:
                        newSection->setType(NEURON_APICAL_DENDRITE);
                        break;

                    case SOMA:
                    case UNKNOWN_SAMPLE:
                        break;
                    }
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

    Vector3f pMaxSoma = _somaCenter + Vector3f(_somaMaxRadius);
    Vector3f pMinSoma = _somaCenter - Vector3f(_somaMaxRadius);

    if (pMaxSoma.x() > _pMax.x()) _pMax.x() = pMaxSoma.x();
    if (pMaxSoma.y() > _pMax.y()) _pMax.y() = pMaxSoma.y();
    if (pMaxSoma.z() > _pMax.z()) _pMax.z() = pMaxSoma.z();
    if (pMinSoma.x() < _pMin.x()) _pMin.x() = pMinSoma.x();
    if (pMinSoma.y() < _pMin.y()) _pMin.y() = pMinSoma.y();
    if (pMinSoma.z() < _pMin.z()) _pMin.z() = pMinSoma.z();
}

AstrocyteMorphology* readAstrocyteMorphology(std::string& morphologyPath)
{
    if (String::subStringFound(morphologyPath, std::string(".swc")) ||
        String::subStringFound(morphologyPath, std::string(".SWC")) ||
        String::subStringFound(morphologyPath, std::string(".Swc")))
    {
        // Read the file
        auto reader = std::make_unique< NeuronSWCReader >(morphologyPath);

        // Get an astrocyte morphology from the reader
        return reader->getAstrocyteMorphology();
    }

     else if (String::subStringFound(morphologyPath, std::string(".h5")) ||
         String::subStringFound(morphologyPath, std::string(".H5")))
     {
#ifdef ULTRALISER_USE_H5
         // Read the file
         auto reader = std::make_unique< AstrocyteH5Reader >(morphologyPath);

         // Get a pointer to the morphology to start using it
         return reader->getMorphology();
#else
         LOG_ERROR("Ultraliser compiled with NO support to read .H5 morphologies!");
         return nullptr;
#endif
     }
     else
     {
         LOG_ERROR("Unrecognized morphology file format [ %s ]", morphologyPath.c_str());
     }

    // To avoid any warning issues.
    return nullptr;
}

}
