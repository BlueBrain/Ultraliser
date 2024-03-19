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
     * @param includeDendriticSample
     */
    SpineMorphology(SkeletonBranch* root, const bool includeDendriticSample = false);

    /**
     * @brief SpineMorphology
     * @param branches
     * @note It's recommended to avoid using this constructor and use the other one to construct
     * the tree.
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
     * @brief getBasePoint
     * @return
     */
    Vector3f getBasePoint() const { return _basePoint; }

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
     * @param showProgrress
     */
    void exportExtents(const std::string& prefix) const;

private:

    /**
     * @brief spineIndex
     */
    size_t _spineIndex;

    float _radfiusScaleFactor = 1.0;

    Sample* _rootSample;



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

private:

    /**
     * @brief _basePoint
     */
    Vector3f _basePoint;
};

/**
 * @brief SpineMorphologies
 */
typedef std::vector< SpineMorphology* > SpineMorphologies;

}
