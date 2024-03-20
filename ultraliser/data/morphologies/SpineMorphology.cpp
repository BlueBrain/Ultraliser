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

#include "SpineMorphology.h"
#include <common/Common.h>
#include <utilities/TypeConversion.h>
#include <utilities/Utilities.h>
#include <algorithms/mcs/DualMarchingCubes.h>

namespace Ultraliser
{

void SpineMorphology::_constructTreeFromLogicalBranches(SkeletonBranch* root,
                                                        size_t& sectionIndex)
{
    for (size_t i = 0; i < root->logicalChildren.size(); ++i)
    {
        Section* section = new Section(sectionIndex++);
        for (size_t j = 0; j < root->logicalChildren[i]->nodes.size(); ++j)
        {
            auto node = root->logicalChildren[i]->nodes[j];
            section->addSample(new Sample(node->point, node->radius * _radfiusScaleFactor, j));
        }
        _sections.push_back(section);

        _constructTreeFromLogicalBranches(root->logicalChildren[i], sectionIndex);
    }
}

SpineMorphology::SpineMorphology(SkeletonBranch* root, const bool includeDendriticSample)
{
    // Update the base point
    _basePoint = root->nodes[0]->point;

    _rootSample = new Sample(root->nodes[0]->point, root->nodes[0]->radius, 0);

    auto dendriteCenter = root->nodes[0]->point;
    auto dendriteExtent = root->nodes[0]->radius;

    // Set the spine index
    _spineIndex = root->index;

    // The section index should be the same as the index of the section in the list
    size_t sectionIndex = 0;

    // Add the root to the list of the sections
    Section* section = new Section(sectionIndex++);
    for (size_t i = 0; i < root->nodes.size(); ++i)
    {
        auto node = root->nodes[i];
        section->addSample(new Sample(node->point, node->radius * _radfiusScaleFactor, i));

    }
    _sections.push_back(section);

    _constructTreeFromLogicalBranches(root, sectionIndex);

    // Compute the bounding box of the entire morphology
    _computeBoundingBox();
}

Sections SpineMorphology::_getNonDendrticSections() const
{
    Sections nonDendriticSections;

    auto dendriticCenter = _rootSample[0].getPosition();
    auto dendriticExtent = _rootSample[0].getRadius();

    size_t validSectionIndex = 0;

    for (size_t i = 0; i < _sections.size(); ++i)
    {
        auto section = _sections[i];
        size_t numberValidSamples = 0;
        for (size_t j = 0; j < section->getSamples().size(); ++j)
        {
            auto sample = section->getSamples()[j];
            auto sampleCenter = sample->getPosition();
            if (!Utilities::isPointInsideSphere(sampleCenter, dendriticCenter, dendriticExtent))
            {
                numberValidSamples++;
            }
        }

        // If more than a sample, then it is a valid section
        if (numberValidSamples > 1)
        {
            // Construct the section
            Section* validSection = new Section(validSectionIndex);
            validSectionIndex++;

            for (size_t j = 0; j < section->getSamples().size(); ++j)
            {
                auto sample = section->getSamples()[j];
                auto sampleCenter = sample->getPosition();
                if (!Utilities::isPointInsideSphere(sampleCenter, dendriticCenter, dendriticExtent))
                {
                    Sample* newSample = new Sample(sample->getPosition(), sample->getRadius() * 2, j);
                    validSection->addSample(sample);
                }
            }
            nonDendriticSections.push_back(validSection);
        }
    }

    return nonDendriticSections;
}

SpineMorphology::SpineMorphology(SkeletonBranches branches, const size_t &index)
{
    _spineIndex = index;

    for (size_t i = 0; i < branches.size(); ++i)
    {
        const auto branch = branches[i];
        Section* section = new Section(i);
        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            auto node = branch->nodes[j];
            section->addSample(new Sample(node->point, node->radius, j));
        }
        _sections.push_back(section);
    }

    // Compute the bounding box of the entire morphology
    _computeBoundingBox();

    /// TODO: Find an algorithm to compute the base.
    /// For the moment, use the center of the bounding box.
    auto bounds = _pMax - _pMin;
    _basePoint = _pMin + 0.5 * bounds;
}

void SpineMorphology::_computeBoundingBox()
{
    // Bounding box data
    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(-1 * std::numeric_limits<float>::max());

    for (const auto& section: _sections)
    {
        for (const auto& sample: section->getSamples())
        {
            const auto position = sample->getPosition();
            const auto radius = sample->getRadius();

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
}

Volume* SpineMorphology::reconstructNonDendriticVolume(const float& voxelsPerMicron,
                                                       const float& edgeGap,
                                                       const bool & verbose)
{
    auto nonDendriticSections = _getNonDendrticSections();

    std::cout << nonDendriticSections.size() << " Non dendritic sections\n";

    if (nonDendriticSections.size() == 0)
        return nullptr;

    // Bounding box data
    Vector3f pMinInput = Vector3f(std::numeric_limits<float>::max());
    Vector3f pMaxInput = Vector3f(-1 * std::numeric_limits<float>::max());

    for (const auto& section: nonDendriticSections)
    {
        for (const auto& sample: section->getSamples())
        {
            const auto position = sample->getPosition();
            const auto radius = sample->getRadius();

            Vector3f pMaxSample = position + Vector3f(radius);
            Vector3f pMinSample = position - Vector3f(radius);

            if (pMaxSample.x() > pMaxInput.x()) pMaxInput.x() = pMaxSample.x();
            if (pMaxSample.y() > pMaxInput.y()) pMaxInput.y() = pMaxSample.y();
            if (pMaxSample.z() > pMaxInput.z()) pMaxInput.z() = pMaxSample.z();

            if (pMinSample.x() < pMinInput.x()) pMinInput.x() = pMinSample.x();
            if (pMinSample.y() < pMinInput.y()) pMinInput.y() = pMinSample.y();
            if (pMinSample.z() < pMinInput.z()) pMinInput.z() = pMinSample.z();
        }
    }

    Vector3f inputBB = pMaxInput - pMinInput;
    Vector3f inputCenter = pMinInput + (0.5f * inputBB);

    // Expand the bounding box to be able to caprture all the details missed
    pMinInput -= edgeGap * inputBB;
    pMaxInput += edgeGap * inputBB;
    inputBB = pMaxInput - pMinInput;
    inputCenter = pMinInput + (0.5f * inputBB);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, edgeGap, VOLUME_TYPE::BIT, verbose);

    // Rasterize the morphologies into the volume
    volume->surfaceVoxelizeSections(nonDendriticSections, verbose);

    // Return the volume
    return volume;
}
Volume* SpineMorphology::reconstructVolume(const float& voxelsPerMicron,
                                           const float& edgeGap,
                                           const bool & verbose)
{
    // Get the bounding box of the morphology
    Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
    getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

    // Expand the bounding box to be able to caprture all the details missed
    _pMin -= edgeGap * inputBB;
    _pMax += edgeGap * inputBB;
    inputBB = _pMax - _pMin;
    inputCenter = _pMin + (0.5f * inputBB);

    // Get the largest dimension
    float largestDimension = inputBB.getLargestDimension();
    size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

    // Construct the volume
    Volume* volume = new Volume(pMinInput, pMaxInput, resolution, edgeGap, VOLUME_TYPE::BIT, verbose);

    // Rasterize the morphologies into the volume
    volume->surfaceVoxelizeSpineMorphology(this, POLYLINE_SPHERE_PACKING);

    // Return the volume
    return volume;
}

Mesh* SpineMorphology::reconstructMesh(const float &voxelsPerMicron,
                                       const float& edgeGap,
                                       const bool &verbose)
{
    // Reconstruct the volume
    auto volume = reconstructVolume(voxelsPerMicron, edgeGap, verbose);

    // Use the DMC algorithm to reconstruct a mesh
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume, verbose);

    // Smooth the mesh to be able to have correct mapping
    mesh->smoothSurface(10, verbose);

    // Return the mesh
    return mesh;
}

Mesh* SpineMorphology::reconstructNonDendriticMesh(const float &voxelsPerMicron,
                                                   const float& edgeGap,
                                                   const bool &verbose)
{
    // Reconstruct the volume
    auto volume = reconstructNonDendriticVolume(voxelsPerMicron, edgeGap, verbose);

    if (volume == nullptr)
    {
        std::cout << "return null ptr\n";
        return nullptr;
    }

    // Use the DMC algorithm to reconstruct a mesh
    auto mesh = DualMarchingCubes::generateMeshFromVolume(volume, verbose);

    // Smooth the mesh to be able to have correct mapping
    mesh->smoothSurface(10, verbose);

    // Return the mesh
    return mesh;
}

void SpineMorphology::exportExtents(const std::string& prefix) const
{
    // Construct the file path
    std::stringstream sstream;
    sstream << prefix << "_spine_" <<_spineIndex << ".extents";
    std::fstream stream;
    stream.open(sstream.str(), std::ios::out);

    // Compute the BB
    auto bounds = _pMax - _pMin;
    auto center = _pMin + bounds * 0.5f;

    // Export the data
    stream << center.x() << " "
           << center.y() << " "
           << center.z() << " "
           << bounds.x() << " "
           << bounds.y() << " "
           << bounds.z() << "\n";

    // Close the file
    stream.close();
}

void SpineMorphology::exportBranches(const std::string &prefix,
                                     const bool& verbose)
{
    // Start the timer
    TIMER_SET;

    std::stringstream sstream;
    sstream << prefix << _spineIndex << ".branches";
    std::string filePath = sstream.str(); // prefix + Utilities::number2string(_spineIndex);
    if (verbose) LOG_STATUS("Exporting Spine Branches: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    if (verbose) LOOP_STARTS("Writing Spine Branches");
    for (size_t i = 0; i < _sections.size(); ++i)
    {
        auto section = _sections[i];

        // The @start marks a new branch in the file
        stream << "start " << section->getIndex() << "\n";

        for (auto& sample: section->getSamples())
        {
            stream << sample->getPosition().x() << " "
                   << sample->getPosition().y() << " "
                   << sample->getPosition().z() << " "
                   << sample->getRadius() << "\n";
        }
        // The @end marks the terminal sample of a branch
        stream << "end\n";

        if (verbose) LOOP_PROGRESS(i, _sections.size());
    }
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

}
