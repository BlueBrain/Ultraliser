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

#include "NeuronMorphology.h"

#include <stack>

namespace Ultraliser
{

NeuronMorphology::NeuronMorphology(const NeuronSWCSamples& swcSamples)
{
    // Construct the sections
    _constructMorphology(swcSamples);
}

Sections NeuronMorphology::getFirstSections() const
{
    return _firstSections;
}

void NeuronMorphology::_constructMorphology(const NeuronSWCSamples& swcSamples)
{
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples and the sections list from the tree structure
    std::stack<NeuronSWCSample*> samplesStack;
    std::stack<Section*> sectionStack;
    uint64_t sectionId = 0;

    samplesStack.push(swcSamples[0]);
    sectionStack.push(nullptr);
    while(!samplesStack.empty())
    {
        auto swcSample = samplesStack.top();
        samplesStack.pop();
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

        Sample* sample = new Sample( position, radius);
        
        if (swcSample->type == SWCSampleType::SOMA)
        {
            _somaSamples.push_back(sample);
            for (auto child: swcSample->childrenSamples)
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
        else
        {
            _samples.push_back(sample);
            section->addSample(sample);
            auto children = swcSample->childrenSamples;
            if (children.size() == 1)
            {
                samplesStack.push(children[0]);
                sectionStack.push(section);
            }
            else
            {
                for (auto child: children)
                {
                    samplesStack.push(child);
                    section->addChildIndex(sectionId);
                    auto newSection = new Section(sectionId);
                    _sections.push_back(newSection);
                    ++sectionId;
                    newSection->addParentIndex(section->getIndex());
                    newSection->addSample(sample);
                    sectionStack.push(newSection);
                }
            }
        }
    }
}

}