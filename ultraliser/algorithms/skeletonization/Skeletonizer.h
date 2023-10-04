/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

#include <data/meshes/simple/Mesh.h>
#include <data/volumes/Volume.h>
#include <algorithms/skeletonization/SkeletonNode.hh>
#include <algorithms/skeletonization/SkeletonEdge.hh>
#include <algorithms/skeletonization/SkeletonBranch.h>
#include <common/Headers.hh>
#include <math/Vector3f.h>

namespace Ultraliser
{
/**
 * @brief The Skeletonizer class
 */
class Skeletonizer
{
public:

    /**
     * @brief Skeletonizer
     * @param volume
     */
    Skeletonizer(Volume *volume, const bool &useAcceleration = true);
    ~Skeletonizer() { }

    /**
     * @brief initialize
     * Initialize the skeletonizer.
     */
    void initialize();

    /**
     * @brief skeletonizeVolumeToCenterLines
     */
    virtual void skeletonizeVolumeToCenterLines() = 0;

    /**
     * @brief skeletonizeVolumeBlockByBlock
     * @param blockSize
     * @param numberOverlappingVoxels
     * @param numberZeroVoxels
     */
    virtual void skeletonizeVolumeBlockByBlock(const size_t& blockSize = 512,
                                       const size_t& numberOverlappingVoxels = 25,
                                       const size_t& numberZeroVoxels = 5) = 0;



    std::vector< Vector3f > getShellPoints();




    virtual void constructGraph();



    virtual void segmentComponents() = 0;


    void applyVolumeThinningWithDomainDecomposition();

    void thinVolumeBlockByBlock(const size_t& blockSize = 512,
                                const size_t& numberOverlappingVoxels = 25,
                                const size_t& numberZeroVoxels = 5);

    void applyVolumeThinningToVolume(Volume* volume, const bool &displayProgress = true);

protected:

    /**
     * @brief _computeShellPoints
     * Computes the shell points of the volume, i.e. the points that are representing the surface
     * shell of the volume, or those on its boundary and not internal.
     */
    void _computeShellPoints();

    /**
     * @brief _computeShellPointsUsingAcceleration
     * Computes the shell points of the volume, i.e. the points that are representing the surface
     * shell of the volume, or those on its boundary and not internal. This method uses the
     * acceleration data structures to improve the performance, but indeed requires more memory.
     * @param thinningVoxels
     * The ThinningVoxels list that contains all the filled voxels that will be used directly
     * to skeletonize the volume.
     */
    void _computeShellPointsUsingAcceleration(ThinningVoxelsUI16List &thinningVoxels);

    /**
     * @brief _scaleShellPoints
     * Scales the extracted shell points to ensure their mapping to the input mesh, if any.
     */
    void _scaleShellPoints();

    /**
     * @brief _applyVolumeThinning
     * Apply the VolumeThinning kernels to extract the center-lines of the volume.
     */
    void _applyVolumeThinning();

    /**
     * @brief _applyVolumeThinningUsingAcceleration
     * Apply the VolumeThinning kernels to extract the center-lines of the volume, using the
     * acceleration data structures that were built during the initialization stage.
     */
    void _applyVolumeThinningUsingAcceleration();

    /**
     * @brief _extractNodesFromVoxels
     * Extract the graph nodes from the centerline voxels remaining in the volume.
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxels();

    /**
     * @brief _extractNodesFromVoxelsNaive
     * Extract the graph nodes from the centerline voxels remaining in the volume using the
     * naive loops. Note that this function is just keep for benchmarking.
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsNaive();

    /**
     * @brief _extractNodesFromVoxelsUsingSlicing
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsUsingSlicing();

    /**
     * @brief _extractNodesFromVoxelsUsingAcceleration
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsUsingAcceleration();

    SkeletonBranch* _buildBranch(SkeletonNode* firstNode, SkeletonNode* edgeNode);

    void _buildAcyclicTree(SkeletonBranch* branch, SkeletonBranches &branches);





    void _inflateNodes();

    void _connectNodes(const std::map<size_t, size_t> &indicesMapper);
    void _removeTriangleLoops();

    void _buildBranchesFromNodes(const SkeletonNodes& nodes);

protected:

    /**
     * @brief _volume
     * Input volume that will be used to extract the skeleton.
     */
    Volume* _volume;

    /**
     * @brief _useAcceleration
     * Use acceleration data structure to improve the performance of the skeletonizer.
     */
    const bool _useAcceleration;

    /**
     * @brief _shellPoints
     * A list of the poins (or voxels) that represent the shell of the given volume.
     */
    std::vector< Vector3f > _shellPoints;

    /**
     * @brief _thinningVoxels
     */
    ThinningVoxelsUI16List _thinningVoxels;

    /**
     * @brief _pMinMesh
     */
    Vector3f _pMinMesh;

    /**
     * @brief _pMaxMesh
     */
    Vector3f _pMaxMesh;

    /**
     * @brief _centerMesh
     */
    Vector3f _centerMesh;

    /**
     * @brief _boundsMesh
     */
    Vector3f _boundsMesh;

    /**
     * @brief _pMinVolume
     */
    Vector3f _pMinVolume;

    /**
     * @brief _pMaxVolume
     */
    Vector3f _pMaxVolume;

    /**
     * @brief _centerVolume
     */
    Vector3f _centerVolume;

    /**
     * @brief _boundsVolume
     */
    Vector3f _boundsVolume;

    /**
     * @brief _scaleFactor
     * Mesh to volume scale factor
     */
    Vector3f _scaleFactor;

    /**
     * @brief _nodes
     * A list of all the nodes in the graph.
     */
    SkeletonNodes _nodes;

    /**
     * @brief _branches
     * The branches of the skeleton representing the reconstructed graph.
     */
    SkeletonBranches _branches;
};

}

