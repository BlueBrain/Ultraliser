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

#include "Morphology.h"

namespace Ultraliser
{

Morphology::Morphology()
{

}

Morphology::~Morphology()
{
    for (auto section: _sections)
    {
        delete section;
    }
    _sections.clear();
    for (auto sample: _samples)
    {
        delete sample;
    }
    _samples.clear();
}

Sections Morphology::getSections() const
{
    return _sections;
}

Paths Morphology::getConnectedPathsFromParentsToChildren(const Section* section) const
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

void Morphology::getBoundingBox(Vector3f& pMin,
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
