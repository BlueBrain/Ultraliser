/***************************************************************************************************
 * Copyright (c) 2016 - 2024
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

#pragma once

#include <data/morphologies/Morphology.h>
#include <data/meshes/Meshes.h>
#include <algorithms/skeletonization/SkeletonBranch.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/// Forward declaration
class Volume;

/**
 * @brief The SpineMorphology class
 * Spine morphology
 */
class SpineMorphology: public Morphology
{
public:

    /**
     * @brief SpineMorphology
     * @param root
     */
    SpineMorphology(SkeletonBranch* root);

    /**
     * @brief SpineMorphology
     * @param branches
     */
    SpineMorphology(SkeletonBranches branches, const size_t& index);

public:

    /**
     * @brief reconstructVolume
     * @param voxelsPerMicron
     * @param verbose
     * @return
     */
    Volume* reconstructVolume(const float& voxelsPerMicron,
                              const float &edgeGap = 0.1,
                              const bool &verbose = false);

    /**
     * @brief reconstructMesh
     * @param voxelsPerMicron
     * @param edgeGap
     * @param verbose
     * @return
     */
    Mesh* reconstructMesh(const float& voxelsPerMicron,
                          const float& edgeGap = 0.1,
                          const bool &verbose = false);

    /**
     * @brief exportBranches
     * @param prefix
     * @param verbose
     */
    void exportBranches(const std::string &prefix, const bool &verbose = false);

    /**
     * @brief exportExtents
     * @param prefix
     */
    void exportExtents(const std::string& prefix) const;

private:

    /**
     * @brief spineIndex
     */
    size_t _spineIndex;

private:

    /**
     * @brief _constructTreeFromLogicalBranches
     * @param root
     * @param sectionIndex
     */
    void _constructTreeFromLogicalBranches(SkeletonBranch* root, size_t &sectionIndex);

    /**
     * @brief _computeBoundingBox
     */
    void _computeBoundingBox();
};

}
