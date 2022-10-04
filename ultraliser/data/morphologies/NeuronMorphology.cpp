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

#include "NeuronMorphology.h"
#include <data/morphologies/swc/NeuronSWCReader.h>
#include <data/morphologies/h5/NeuronH5Reader.h>
#include <data/morphologies/swc/NeuronSWCSection.hh>
#include <utilities/String.h>
#include <stack>

namespace Ultraliser
{

NeuronMorphology::NeuronMorphology(const NeuronSWCSamples& swcSamples)
{
    // Construct the morphology from the loaded data
    _constructMorphologyFromSWC(swcSamples);
}

NeuronMorphology::NeuronMorphology(const H5Samples& h5Samples, const H5Sections& h5Sections)
{
    // Construct the morphology from the loaded data
    _constructMorphologyFromH5(h5Samples, h5Sections);
}

void NeuronMorphology::trim(size_t axonBranchOrder,
                            size_t basalBranchOrder,
                            size_t apicalBranchOrder)
{
    Sections newSections;
    Sections newFirstSections;
    _samples.clear();
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());
 
    for (uint8_t i = 0; i < _rootSections.size(); ++i)
    {
        size_t maxDepth = 0;
        switch (_rootSections[i]->getType())
        {
        case PROCESS_TYPE::NEURON_AXON:
            maxDepth = axonBranchOrder;
            break;
        case PROCESS_TYPE::NEURON_BASAL_DENDRITE:
            maxDepth = basalBranchOrder;
            break;
        case PROCESS_TYPE::NEURON_APICAL_DENDRITE:
            maxDepth = apicalBranchOrder;
            break;
        case PROCESS_TYPE::SOMA:
        case PROCESS_TYPE::GLIA_PERIVASCULAR_PROCESS:
        case PROCESS_TYPE::GLIA_PROCESS:
        case PROCESS_TYPE::VASCULATURE:
        case PROCESS_TYPE::UNKNOWN_PROCESS:
        default:
            break;
        }
        
        std::stack< std::pair< Section*, size_t > > sectionStack;
        if (maxDepth > 0)
        {        
            sectionStack.push(std::make_pair(_rootSections[i], 0));
            newFirstSections.push_back(_rootSections[i]);
        }

        while(!sectionStack.empty())
        {
            auto section = sectionStack.top().first;
            auto depth = sectionStack.top().second;
            sectionStack.pop();
            
            auto childrenIndices = section->getChildrenIndices(); 
            auto parentsIndices = section->getParentIndices();

            for (auto index: childrenIndices)
                sectionStack.push(std::make_pair(_sections[index], depth + 1));

            if (depth < maxDepth)
            {
                size_t index = newSections.size();
                section->setIndex(index);
                section->clearChildrenIndices();
                newSections.push_back(section);
                for (auto childIndex: childrenIndices)
                {
                    auto childSection = _sections[childIndex];
                    childSection->clearParentsIndices();
                    childSection->addParentIndex(index);
                }   
                for (auto parentIndex: parentsIndices)
                {
                    auto parentSection = newSections[parentIndex];
                    parentSection->addChildIndex(index);
                }
               
                auto sectionSamples = section->getSamples();
                _samples.insert(_samples.end(), sectionSamples.begin(), sectionSamples.end());
                for (auto sample: sectionSamples)
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
            else
                delete section;
        }
    }

    Vector3f pMaxSoma = _somaCenter + Vector3f(_somaMaxRadius);
    Vector3f pMinSoma = _somaCenter - Vector3f(_somaMaxRadius);

    if (pMaxSoma.x() > _pMax.x()) _pMax.x() = pMaxSoma.x();
    if (pMaxSoma.y() > _pMax.y()) _pMax.y() = pMaxSoma.y();
    if (pMaxSoma.z() > _pMax.z()) _pMax.z() = pMaxSoma.z();
    if (pMinSoma.x() < _pMin.x()) _pMin.x() = pMinSoma.x();
    if (pMinSoma.y() < _pMin.y()) _pMin.y() = pMinSoma.y();
    if (pMinSoma.z() < _pMin.z()) _pMin.z() = pMinSoma.z();

    _sections.clear();
    _sections = newSections;
    _rootSections.clear();
    _rootSections = newFirstSections;
}

Samples NeuronMorphology::getSomaSamples() const
{
    return _somaSamples;
}

Sections NeuronMorphology::getFirstSections() const
{
    return _rootSections;
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

// Reindexing of acyclic morphologies, with root nodes
void NeuronMorphology::reIndexMorphology()
{
    // Global indexing for the samples
    size_t sampleIndex = 0;

    // A collector for all the samples in the morphology including the ones that have been added
    // We will clear the current samples list, and then we reconstruct a new list of samples
    _samples.clear();
    _samples.shrink_to_fit();

    // Add a zeroSample, for indexing
    Sample* zeroSample = new Sample(Vector3f(0.f), 0.f, PROCESS_TYPE::UNKNOWN_PROCESS, sampleIndex++, 0);
    _samples.push_back(zeroSample);

    // Add the somaSample
    Sample* somaSample = new Sample(_somaCenter, 0.5, PROCESS_TYPE::SOMA, sampleIndex++, -1);
    _samples.push_back(somaSample);

    // Get root sections, for neurons and astrocytes
    auto rootSections = getRootSections();

    // For each root section
    for (size_t i = 0; i < rootSections.size(); ++i)
    {
        // Get a reference to the root section
        auto& rootSection = rootSections[i];

        // Set the indices of the samples of the section
        auto rootSamples = rootSection->getSamples();

        // The first sample of the rootSection, set index, parent index and append to the list
        rootSamples.front()->setIndex(sampleIndex++);
        rootSamples.front()->setParentIndex(1);
        _samples.push_back(rootSamples.front());

        // Intermediate samples
        for (size_t j = 1; j < rootSamples.size(); ++j)
        {
            // Set index, parent index and append to the list
            rootSamples[j]->setIndex(sampleIndex);
            rootSamples[j]->setParentIndex(sampleIndex - 1);
            sampleIndex++;
            _samples.push_back(rootSamples[j]);
        }        

        // The branchingSampleIndex will be the same for all the children sections
        const size_t branchingSampleIndex = _samples[_samples.size() - 1]->getIndex();

        auto childrenIndices = rootSection->getChildrenIndices();
        for (size_t j = 0; j < childrenIndices.size(); ++j)
        {
            // Get a reference to the child section
            auto child = _sections[childrenIndices[j]];

            // Apply recursively
            child->reIndexSectionTree(_sections, sampleIndex, branchingSampleIndex, _samples);
        }
    }
}













void NeuronMorphology::_constructMorphologyFromSWC(const NeuronSWCSamples& swcSamples)
{
    // Construct the _samples list, in order
    for (size_t i = 0; i < swcSamples.size(); ++i)
    {
        const auto swcSample = swcSamples[i];
        Sample* sample = new Sample(Vector3f(swcSample->x, swcSample->y, swcSample->z),
                                    swcSample->r, swcSample->type,
                                    swcSample->id, swcSample->parentId);
        _samples.push_back(sample);
    }

    size_t numberSomaSamples = 0;
    size_t numberSomaProfilePoints = 0;

    typedef struct
    {
        // The index of the parent sample
        size_t parentSample;

        // The index of the first sample of the arbor
        size_t firstSample;

        // The index of the last sample of the arbor
        size_t lastSample;

        std::vector< size_t > branchingIndices;
    } Indices;

    // The first element represents the first sample, and the second represents the last sample
    std::vector< Indices > arborsIndices;

    // Pre-process the samples, to get the construction
    for (size_t i = 0; i < swcSamples.size(); ++i)
    {
        if (swcSamples[i]->type == PROCESS_TYPE::SOMA &&
            swcSamples[i]->id == 1 && swcSamples[i]->parentId == -1)
        {
            // This is a soma sample
            numberSomaSamples++;
        }
        else if (swcSamples[i]->type == PROCESS_TYPE::SOMA && swcSamples[i]->parentId == 1)
        {
            // This is a soma profile point
            numberSomaProfilePoints++;
        }
        else if (swcSamples[i]->type != PROCESS_TYPE::SOMA && swcSamples[i]->parentId == 1)
        {
            // Add the first sample of the arbor
            Indices arbor;
            arbor.firstSample = swcSamples[i]->id;
            arborsIndices.push_back(arbor);
        }
    }

    // Get the last sampels of the arbor
    for (size_t i = 0; i < arborsIndices.size(); ++i)
    {
        // The index of the last sample of the current arbor can be either obtained from the index
        // of the first sample of the next sample, or from the swcSamples, only for the last arbor
        arborsIndices[i].lastSample = (i == arborsIndices.size() - 1) ?
                    swcSamples.back()->id : arborsIndices[i + 1].firstSample - 1;
    }


    for (size_t i = 0; i < arborsIndices.size(); ++i)
    {
        LOG_WARNING("Arbor %d, S0 %d, S1 %d", i, arborsIndices[i].firstSample, arborsIndices[i].lastSample);
    }
    LOG_WARNING("**********************************");

    // Construct the sections
    size_t sectionIndex = 0;
    for (size_t i = 0; i < arborsIndices.size(); ++i)
    {
         // LOG_WARNING("Arbor %d, S0 %d, S1 %d", i, arborsIndices[i].firstSample, arborsIndices[i].lastSample);

        // Get the indices of the terminal samples
        const size_t& firstSampleIndex = arborsIndices[i].firstSample;
        const size_t& lastSampleIndex = arborsIndices[i].lastSample;

        // Collect the indices of the samples representing new paths
        std::vector< size_t > newPathsIndices;
        for (size_t j = firstSampleIndex + 1; j <= lastSampleIndex; ++j)
        {
            if (swcSamples[j]->id != (swcSamples[j]->parentId + 1))
            {
                newPathsIndices.push_back(j);
                // LOG_SUCCESS("A new path starts at [%d]", j);
            }
        }

        // Collect the indices of the paths, parent and along the path
        // std::vector< Indices > pathsIndices;
        std::vector< size_t > firstSamplesIndices;



        std::vector< Indices > sectionsIndices;

        // In case of no segments, then the section emanating from the soma has no children
        if (newPathsIndices.size() == 0)
        {
            // Therefore, add the indices from the first and last samples directly
            Indices path;
            path.parentSample = 1;
            path.firstSample = firstSampleIndex;
            path.lastSample = lastSampleIndex;

            // The path itself is a section
            sectionsIndices.push_back(path);
            std::cout << path.parentSample << ": "
                      << path.firstSample << " | " << path.lastSample << "\n";







            firstSamplesIndices.push_back(path.firstSample);

        }
        else if (newPathsIndices.size() == 1)
        {
            // The first path
            Indices firstPath;
            firstPath.parentSample = 1;
            firstPath.firstSample = arborsIndices[i].firstSample;
            firstPath.lastSample = newPathsIndices[0] - 1;
            firstSamplesIndices.push_back(firstPath.firstSample);

            // The second path
            Indices lastPath;
            lastPath.parentSample = swcSamples[newPathsIndices[0]]->parentId;
            lastPath.firstSample = newPathsIndices[0];
            lastPath.lastSample = arborsIndices[i].lastSample;
            firstSamplesIndices.push_back(lastPath.firstSample);

            // Split the first path in two sections and that's it

            // The first path should be split into two sections
            Indices section1, section2;
            section1.parentSample = 1;
            section1.firstSample = firstPath.firstSample;
            section1.lastSample = lastPath.parentSample;

            section2.parentSample = lastPath.parentSample;
            section2.firstSample = lastPath.parentSample + 1;
            section2.lastSample = firstPath.lastSample;

            sectionsIndices.push_back(section1);
            sectionsIndices.push_back(section2);
            sectionsIndices.push_back(lastPath);

            std::cout << section1.parentSample << ": "
                      << section1.firstSample << " | " << section1.lastSample << "\n";

            std::cout << section2.parentSample << ": "
                      << section2.firstSample << " | " << section2.lastSample << "\n";

            std::cout << lastPath.parentSample << ": "
                      << lastPath.firstSample  << " | " << lastPath.lastSample << "\n";
        }
        else
        {
            std::vector< Indices > pathsIndices;

            std::vector< size_t > branchingPointsIndices;

            // The first path
            Indices initialPath;
            initialPath.parentSample = 1;
            initialPath.firstSample = arborsIndices[i].firstSample;
            initialPath.lastSample = newPathsIndices[0] - 1;
            pathsIndices.push_back(initialPath);
            firstSamplesIndices.push_back(initialPath.firstSample);

            // The rest of the paths
            for (size_t j = 0; j < newPathsIndices.size(); ++j)
            {
                if (j < newPathsIndices.size() - 1)
                {
                    // Get the current index
                    size_t currentIndex = newPathsIndices[j];

                    // Get the parentId @swcSample[currentIndex]
                    size_t parentIndex = swcSamples[currentIndex]->parentId;

                    // Get the last sample
                    size_t lastSampleIndex = newPathsIndices[j + 1] - 1;

                    Indices inbetweenPath;
                    inbetweenPath.parentSample = parentIndex;
                    inbetweenPath.firstSample = currentIndex;
                    inbetweenPath.lastSample = lastSampleIndex;
                    pathsIndices.push_back(inbetweenPath);
                    firstSamplesIndices.push_back(inbetweenPath.firstSample);
                    branchingPointsIndices.push_back(parentIndex);
                }
                else
                {
                    // Get the current index
                    size_t currentIndex = newPathsIndices[j];

                    // Get the parentId @swcSample[currentIndex]
                    size_t parentIndex = swcSamples[currentIndex]->parentId;

                    // Get the last sample
                    size_t lastSampleIndex = arborsIndices[i].lastSample;

                    Indices inbetweenPath;
                    inbetweenPath.parentSample = parentIndex;
                    inbetweenPath.firstSample = currentIndex;
                    inbetweenPath.lastSample = lastSampleIndex;
                    pathsIndices.push_back(inbetweenPath);
                    firstSamplesIndices.push_back(inbetweenPath.firstSample);
                    branchingPointsIndices.push_back(parentIndex);
                }
            }

            // Sort
            std::sort(branchingPointsIndices.begin(), branchingPointsIndices.end());

            /// After getting the list of paths, simply make the sections

            for (size_t j = 0; j < pathsIndices.size(); ++j)
            {
                auto path = pathsIndices[j];
//                LOG_SUCCESS("* Arbor: [%d], Path: [%d], [%d], [%d]", i,
//                            path.parentSample, path.firstSample, path.lastSample);

                // Collect the branching indices
                for (size_t k = 0; k < branchingPointsIndices.size(); ++k)
                {
                    if (branchingPointsIndices[k] >= path.firstSample &&
                            branchingPointsIndices[k] <= path.lastSample)
                    {
                        path.branchingIndices.push_back(branchingPointsIndices[k]);
                    }
                }

                // Sort the branching indices
                std::sort(path.branchingIndices.begin(), path.branchingIndices.end());

                if (path.branchingIndices.size() == 0)
                {
                    sectionsIndices.push_back(path);
                    std::cout << path.parentSample << ": "
                              << path.firstSample  << " | " << path.lastSample << " *\n";

                }
                else if (path.branchingIndices.size() == 1)
                {
                    Indices section1, section2;
                    section1.parentSample = path.parentSample;
                    section1.firstSample = path.firstSample;
                    section1.lastSample = path.branchingIndices[0];

                    section2.parentSample = path.branchingIndices[0];
                    section2.firstSample = path.branchingIndices[0] + 1;
                    section2.lastSample = path.lastSample;

                    sectionsIndices.push_back(section1);
                    sectionsIndices.push_back(section2);

                    std::cout << section1.parentSample << ": "
                              << section1.firstSample << " | " << section1.lastSample << " **\n";

                    std::cout << section2.parentSample << ": "
                              << section2.firstSample << " | " << section2.lastSample << "*** \n";
                }
                else
                {
                    // Construct the sections
                    for (size_t k = 0; k <= path.branchingIndices.size(); ++k)
                    {
                        if (k == 0)
                        {
                            if (path.parentSample == 1)
                            {
                                Indices section;
                                section.parentSample = path.parentSample;
                                section.firstSample = path.firstSample;
                                section.lastSample = path.branchingIndices[k];
                                sectionsIndices.push_back(section);

                                std::cout << section.parentSample << ": "
                                          << section.firstSample << " | " << section.lastSample << " ****\n";
                            }
                            else
                            {
                                Indices section;
                                section.parentSample = path.parentSample;
                                section.firstSample = path.firstSample;
                                section.lastSample = path.branchingIndices[k];
                                sectionsIndices.push_back(section);

                                std::cout << section.parentSample << ": "
                                          << section.firstSample << " | " << section.lastSample << "***** \n";
                            }
                        }
                        else if (k == path.branchingIndices.size())
                        {

                            Indices section;
                            section.parentSample = path.branchingIndices[k - 1];
                            section.firstSample = path.branchingIndices[k - 1] + 1;
                            section.lastSample = path.lastSample;
                            sectionsIndices.push_back(section);

                            std::cout << section.parentSample << ": "
                                      << section.firstSample << " | " << section.lastSample << "****** \n";
                        }
                        else
                        {
                            Indices section;
                            section.parentSample = path.branchingIndices[k - 1];
                            section.firstSample = path.branchingIndices[k - 1] + 1;
                            section.lastSample = path.branchingIndices[k];
                            sectionsIndices.push_back(section);

                            std::cout << section.parentSample << ": "
                                      << section.firstSample << " | " << section.lastSample << "*******\n";
                        }
                    }
                }
            }
        }

        // Construct the sections

        for (size_t j = 0; j < sectionsIndices.size(); ++j)
        {
            auto& swcSection =  sectionsIndices[j];

            if (swcSection.parentSample == 1)
                LOG_SUCCESS("swcSection: [%d], Parent Sample: [%d], First Sample %d", sectionIndex, swcSection.parentSample, swcSection.firstSample);

            // Construct the new section, assign ID and type
            Section* section = new Section(sectionIndex, swcSamples[swcSection.firstSample]->type);
            sectionIndex++;

            // Construct samples
            Samples samples;

            // Soma samples are not included in the samples list
            if (swcSection.parentSample == 1) { }
            else
            {
                samples.push_back(new Sample(swcSamples[swcSection.parentSample]));
            }

            // Add the rest of the samples
            for (size_t k = swcSection.firstSample; k <= swcSection.lastSample; ++k)
            {
                samples.push_back(new Sample(swcSamples[k]));
            }

            section->addSamples(samples);

            // Add the section to the sections list
            _sections.push_back(section);
        }
    }


    LOG_SUCCESS("Sections %d", _sections.size());

//    for (size_t i = 0; i < _sections.size(); ++i)
//    {
//        LOG_WARNING("Section [%d]: %d | %d", _sections[i]->getIndex(),
//                    _sections[i]->getFirstSample()->getIndex(),
//                    _sections[i]->getLastSample()->getIndex());
//    }




    // Construct the graph from the individual sections
    /// Note that the first sample is AUXILIARY, and will NOT be taken into consideration
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        auto& iSection = _sections[i];

        for (size_t j = 0; j < _sections.size(); ++j)
        {
            auto& jSection = _sections[j];

            // Same section, ignore
            if (i == j)
                continue;

            if (iSection->getFirstSample()->getIndex() == jSection->getLastSample()->getIndex())
            {
                iSection->addParentIndex(jSection->getIndex());
                jSection->addChildIndex(iSection->getIndex());
            }

            if (iSection->getLastSample()->getIndex() == jSection->getFirstSample()->getIndex())
            {

                iSection->addChildIndex(jSection->getIndex());
                jSection->addParentIndex(iSection->getIndex());
            }
        }
    }


    // Get the root sections, for quick access
    /// Note that the first sample is AUXILIARY, and will NOT be taken into consideration
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        if (_sections[i]->isRoot())
        {
            _rootSections.push_back(_sections[i]);
            LOG_WARNING("Root %d, %d", _sections[i]->getIndex(), _sections[i]->getFirstSample()->getIndex());
        }
    }

    LOG_WARNING("Soma: %d, Profiles: %d, Sections: %d",
                numberSomaSamples, numberSomaProfilePoints, _rootSections.size());

    // The _pMin and _pMax will be used to compute the bounding box of the neuron morphology
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    for (size_t i = 1; i < _samples.size(); ++i)
    {
        Vector3f position = _samples[i]->getPosition();
        float radius = _samples[i]->getRadius();

        Vector3f pMaxSample = position + Vector3f(radius);
        Vector3f pMinSample = position - Vector3f(radius);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }

    // Get the somatic samples from the initial sections of all the neurites
    for (const auto& section: _rootSections)
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

void NeuronMorphology::_constructMorphologyFromH5(const H5Samples& h5Samples,
                                                  const H5Sections& h5Sections)
{
    // The _pMin and _pMax will be used to compute the bounding box of the neuron morphology
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Compute the boundinx box from the H5 samples
    for (auto h5Sample : h5Samples)
    {
        Ultraliser::Vector3f position(h5Sample.x, h5Sample.y, h5Sample.z);

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
    Section* somaSection = new Section(0, PROCESS_TYPE::SOMA);
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

        // Get the process type
        auto processType = mapNeuronH5IndexToType(h5Sections[i].sectionType);

        // Create a new morphology section
        Section* section = new Section(sectionIndex, processType);

        // Construct a list of samples
        Samples samples;
        for (size_t s = firstPointIdx; s <= lastPointIdx; s++)
        {
            const auto sampleIndex = s;
            size_t parentIndex;

            // Create the sample
            Sample* sample =
                    new Sample(Ultraliser::Vector3f(h5Samples[s].x, h5Samples[s].y, h5Samples[s].z),
                               h5Samples[s].r * 0.5f, processType, s, parentIndex);

            // Update the samples list
            section->addSample(sample);

            _samples.push_back(sample);
        }

        // Update the section parent
        section->addParentIndex(h5Sections[i].parentIndex);

        // Add the section to the list
        _sections.push_back(section);

        // Determine if this is a first section or not!
        if (h5Sections[i].parentIndex == 0)
        {
            _rootSections.push_back(section);
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
    for (const auto& section: _rootSections)
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
    else if (String::subStringFound(morphologyPath, std::string(".h5")) ||
             String::subStringFound(morphologyPath, std::string(".H5")))
    {
#ifdef ULTRALISER_USE_H5
        // Read the file
        auto reader = std::make_unique< NeuronH5Reader >(morphologyPath);

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
