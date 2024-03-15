/***************************************************************************************************
 * Copyright (c) 2016 - 2024
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <algorithms/skeletonization/Skeletonizer.h>
#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{

/**
 * @brief The NeuronSkeletonizer class
 */
class SpineSkeletonizer : public Skeletonizer
{
public:

    /**
     * @brief SpineSkeletonizer
     * @param spineVolume
     * @param useAcceleration
     * @param debugSkeleton
     * @param debuggingPrefix
     */
    SpineSkeletonizer(Volume* spineVolume,
                      const bool useAcceleration,
                      const bool debugSkeleton,
                      const std::string debuggingPrefix = NONE);

    /**
     * @brief SpineSkeletonizer
     * @param spineMesh
     * @param basePoint
     * @param options
     * @param useAcceleration
     * @param debugSkeleton
     * @param debuggingPrefix
     */
    SpineSkeletonizer(Mesh* spineMesh,
                      const Vector3f basePoint,
                      const VoxelizationOptions& options,
                      const bool useAcceleration,
                      const bool debugSkeleton,
                      const std::string debuggingPrefix = NONE);
    ~SpineSkeletonizer() { }

public:

    /**
     * @brief exportBranches
     * @param prefix
     */
    void exportBranches(const std::string& prefix, const bool = VERBOSE);

    /**
     * @brief run
     * @param verbose
     */
    void run(const bool verbose = VERBOSE);

    /**
     * @brief exportSWCFile
     * @param prefix
     * @param resampleSkeleton
     */
    void exportSWCFile(const std::string& prefix, const bool& resampleSkeleton,
                       const bool verbose = VERBOSE);

    /**
     * @brief segmentComponents
     */
    void segmentComponents(const bool verbose = VERBOSE) override;

private:

    /**
     * @brief _constructSWCTable
     * @param resampleSkeleton
     * @return
     */
    SkeletonNodes _constructSWCTable(const bool& resampleSkeleton, const bool = VERBOSE);

    /**
     * @brief _collectSWCNodes
     * @param branch
     * @param swcNodes
     * @param swcIndex
     * @param branchingNodeSWCIndex
     */
    void _collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                          int64_t& swcIndex, int64_t branchingNodeSWCIndex,
                          const bool = VERBOSE);

private:

    /**
     * @brief _basePoint
     * The base point of the spine. This value is used later to compute the orientation of the spine.
     */
    Vector3f _basePoint;
};

}
