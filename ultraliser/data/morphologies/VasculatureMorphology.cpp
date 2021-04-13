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

#include "VasculatureMorphology.h"

namespace Ultraliser
{

VasculatureMorphology::VasculatureMorphology
(const VasculatureH5Samples& h5Samples, const VasculatureH5Sections& h5Sections,
 const VasculatureH5ConnectivityList& h5Connectivity)
    : _h5Samples(h5Samples)
    , _h5Sections(h5Sections)
    , _h5Connectivity(h5Connectivity)
{
    // Construct the sections
    _constructSections();

    // Connect the sections
    _connectSections();
}

Sections VasculatureMorphology::getSections() const
{
    return _sections;
}

void VasculatureMorphology::_constructSections()
{
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::min());

    // Construct the samples list that will be used later to build the sections
    for (auto h5Sample : _h5Samples)
    {
        Ultraliser::Vector3f position(h5Sample.x, h5Sample.y, h5Sample.z);
        Sample* sample = new Sample(Ultraliser::Vector3f(h5Sample.x,
                                                         h5Sample.y,
                                                         h5Sample.z),
                                    h5Sample.r * 0.5f);
        _samples.push_back(sample);

        Vector3f pMaxSample = position + Vector3f(h5Sample.r * 0.5);
        Vector3f pMinSample = position - Vector3f(h5Sample.r * 0.5);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }

    // Construct the sections list
    for (uint64_t i = 0; i < _h5Sections.size(); ++i)
    {
        Section* section = new Section(i);

        // The initial sample always stays the same
        const uint64_t startingSampleIndex =
                I2UI64(_h5Sections[i].firstSampleIndex);

        // The last sample is different
        uint64_t endSampleIndex;
        if (i < _h5Sections.size() - 1)
            endSampleIndex = I2UI64(_h5Sections[i + 1].firstSampleIndex);
        else
            endSampleIndex = _h5Samples.size();

        // Add the samples to the section
        for (uint64_t j = startingSampleIndex; j < endSampleIndex; ++j)
            section->addSample(_samples[j]);

        // Append the reconstructed section to the list
        _sections.push_back(section);
    }
}

void VasculatureMorphology::_connectSections()
{
    for (uint64_t i = 0; i < _h5Connectivity.size(); ++i)
    {
        // Get the parent section and add to it the index of a child section
        Section* parentSection = _sections[
                I2UI64(_h5Connectivity[i].parentSectionIndex)];
        parentSection->addChildIndex(
                    I2UI64(_h5Connectivity[i].childSectionIndex));

        // Get the child section and add to it the index of a parent section
        Section* childSection = _sections[
                I2UI64(_h5Connectivity[i].childSectionIndex)];
        childSection->addParentIndex(
                    I2UI64(_h5Connectivity[i].parentSectionIndex));
    }
}

std::vector< Samples >
VasculatureMorphology::_getConnectedSamplesFromParentToChild() const
{
    // A list of all the connected samples from the parents to the children
    // to avoid any gaps in the reconstructed morphology before converting it
    // to a mesh
    std::vector< Samples > connectedSamples;

    for (uint64_t i = 0; i < _sections.size(); ++i)
    {
        // Get the section
        Section* section = _sections[i];
        Samples sectionSamples = section->getSamples();

        // Get the number of parents
        const std::vector< uint64_t > parentsIndices =
                section->getParentIndices();
        const uint64_t numberParents = parentsIndices .size();

        // Get the number of children
        const std::vector< uint64_t > childrenIndices =
                section->getChildrenIndices();
        const uint64_t numberChildren = childrenIndices.size();

        // If the section has parents and children
        if (numberParents > 0 && numberChildren > 0)
        {
            for (uint64_t j = 0; j < numberParents; ++j)
            {
                Section* parentSection = _sections[parentsIndices[j]];
                Samples parentSamples = parentSection->getSamples();

                for (uint64_t k = 0; k < numberChildren; ++k)
                {
                    Section* childSection = _sections[childrenIndices[k]];
                    Samples childSamples = childSection->getSamples();

                    // Construct the final list
                    Samples samples;

                    // Parent section samples
                    samples.insert(samples.begin(),
                                   parentSamples.begin(),
                                   parentSamples.end());

                    // Ignore the first and last samples to avoid duplicated
                    // samples along the poly-line. This is mandatory to be
                    // able to draw a VTK poly-line.
                    samples.insert(samples.end(),
                                   sectionSamples.begin() + 1,
                                   sectionSamples.end() - 1);

                    // Child section samples
                    samples.insert(samples.end(),
                                   childSamples.begin(),
                                   childSamples.end());


                    // Append the list to the final output
                    connectedSamples.push_back(samples);
                }
            }
        }

        // If the section is root and has children
        else if (numberParents == 0 && numberChildren > 0)
        {
            for (uint64_t k = 0; k < numberChildren; ++k)
            {
                Section* childSection = _sections[childrenIndices[k]];
                Samples childSamples = childSection->getSamples();

                // Construct the final list
                Samples samples;

                // Section samples
                samples.insert(samples.begin(),
                               sectionSamples.begin(),
                               sectionSamples.end());

                // Child section samples
                samples.insert(samples.end(),
                               childSamples.begin() + 1,
                               childSamples.end());

                // Append the list to the final output
                connectedSamples.push_back(samples);
            }
        }

        // If the section is a leaf node but had parents
        else if (numberParents > 0 && numberChildren == 0)
        {
            for (uint64_t j = 0; j < numberParents; ++j)
            {
                Section* parentSection = _sections[parentsIndices[j]];
                Samples parentSamples = parentSection->getSamples();

                // Construct the final list
                Samples samples;

                // Parent section samples
                samples.insert(samples.begin(),
                               parentSamples.begin(),
                               parentSamples.end());

                // Section samples
                samples.insert(samples.end(),
                               sectionSamples.begin() + 1,
                               sectionSamples.end());

                // Append the list to the final output
                connectedSamples.push_back(samples);
            }
        }

        // If the section is root and lead at the same moment
        else if (numberParents == 0 && numberChildren == 0)
        {
            connectedSamples.push_back(sectionSamples);
        }
        else
        {
            /// TODO: HANDLE THIS
        }
    }

    return connectedSamples;
}

Paths VasculatureMorphology::getConnectedPathsFromParentsToChildren(const Section* section) const
{
    // All possible combination of paths along the section (from parents to children)
    Paths paths;

    // Get the parents data
    const std::vector< uint64_t > parentsIndices =
            section->getParentIndices();
    const uint64_t numberParents = parentsIndices .size();

    // Get the children data
    const std::vector< uint64_t > childrenIndices =
            section->getChildrenIndices();
    const uint64_t numberChildren = childrenIndices.size();

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
        for (uint64_t k = 0; k < numberChildren; ++k)
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
        for (uint64_t j = 0; j < numberParents; ++j)
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
        for (uint64_t j = 0; j < numberParents; ++j)
        {
            // Parents data
            Section* parentSection = _sections[parentsIndices[j]];
            Samples parentSamples = parentSection->getSamples();

            for (uint64_t k = 0; k < numberChildren; ++k)
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

void VasculatureMorphology::getBoundingBox(Vector3f& pMin,
                                           Vector3f& pMax,
                                           Vector3f& bounds,
                                           Vector3f& center)
{
    pMin = _pMin;
    pMax = _pMax;
    bounds = _pMax - _pMin;
    center = _pMin + (0.5f * bounds);
}

}
