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

#include "VasculatureMorphology.h"
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <utilities/TypeConversion.h>
#include <data/morphologies/h5/VasculatureH5Reader.h>
#include <data/morphologies/vmv/VasculatureVMVReader.h>

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

    // Resamplg the sections, relaxed
    resampleSectionsUniformly(0.1);
}

void VasculatureMorphology::_constructSections()
{
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples list that will be used later to build the sections
    for (auto h5Sample : _h5Samples)
    {
        Ultraliser::Vector3f position(h5Sample.x, h5Sample.y, h5Sample.z);
        Sample* sample = new Sample(Ultraliser::Vector3f(h5Sample.x, h5Sample.y, h5Sample.z),
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

    // Construct the sections list
    for (uint64_t i = 0; i < _h5Sections.size(); ++i)
    {
        Section* section = new Section(i);

        // The initial sample always stays the same
        const uint64_t startingSampleIndex = I2UI64(_h5Sections[i].firstSampleIndex);

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
        Section* parentSection = _sections[I2UI64(_h5Connectivity[i].parentSectionIndex)];
        parentSection->addChildIndex(I2UI64(_h5Connectivity[i].childSectionIndex));

        // Get the child section and add to it the index of a parent section
        Section* childSection = _sections[I2UI64(_h5Connectivity[i].childSectionIndex)];
        childSection->addParentIndex(I2UI64(_h5Connectivity[i].parentSectionIndex));
    }
}

VasculatureMorphology::VasculatureMorphology(Samples samples, Sections sections)
{
    _samples = samples;
    _sections = sections;

    // Bounding box data
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());

    // Construct the samples list that will be used later to build the sections
    for (auto sample : _samples)
    {
        const Ultraliser::Vector3f position = sample->getPosition();
        const float radius = sample->getRadius();

        Vector3f pMaxSample = position + Vector3f(radius);
        Vector3f pMinSample = position - Vector3f(radius);

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }
}

VasculatureMorphology* readVascularMorphology(std::string& morphologyPath)
{
    LOG_TITLE("Loading Morphology");

    if (String::subStringFound(morphologyPath, std::string(".h5")) ||
            String::subStringFound(morphologyPath, std::string(".H5")))
    {
        // Read the file
        auto reader = std::make_unique<VasculatureH5Reader>(morphologyPath);

        // Get a pointer to the morphology to start using it
        return reader->getMorphology();
    }
    else if (String::subStringFound(morphologyPath, std::string(".vmv")) ||
             String::subStringFound(morphologyPath, std::string(".VMV")))
    {
        // Read the file
        auto reader = std::make_unique<VasculatureVMVReader>(morphologyPath);

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
