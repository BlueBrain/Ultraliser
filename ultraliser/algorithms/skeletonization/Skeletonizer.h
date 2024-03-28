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

#include <data/volumes/Volume.h>
#include <data/meshes/simple/Mesh.h>
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
     * @brief The VoxelizationOptions class
     */
    struct VoxelizationOptions
    {
        /**
         * @brief volumeResolution
         * The resolution used to voxelize the mesh given to the Skeletonizer.
         */
        size_t volumeResolution = 1024;

        /**
         * @brief voxelizationAxis
         * The solid voxelization axis of the Skeletonizer.
         */
        Volume::SOLID_VOXELIZATION_AXIS voxelizationAxis = Volume::SOLID_VOXELIZATION_AXIS::X;

        /**
         * @brief volumeType
         * The type of the volume, by default BIT.
         */
        VOLUME_TYPE volumeType = VOLUME_TYPE::BIT;

        /**
         * @brief edgeGapPrecentage
         * The percentage of the volume that is used to make a gap between the core volume
         * and its edges.
         */
        float edgeGapPrecentage = 0.1f;

        /**
         * @brief verbose
         * If the operations will be executed silently or not.
         */
        bool verbose = false;
    };

public:

    /**
     * @brief Skeletonizer
     * @param volume
     * @param useAcceleration
     * @param debugSkeleton
     * @param debuggingPrefix
     */
    Skeletonizer(Volume *volume,
                 const bool &useAcceleration = true,
                 const bool &debugSkeleton = false,
                 const std::string debuggingPrefix = NONE);

    /**
     * @brief Skeletonizer
     * @param mesh
     * @param options
     * @param useAcceleration
     * @param debugSkeleton
     * @param debuggingPrefix
     */
    Skeletonizer(Mesh *mesh,
                 const VoxelizationOptions& options,
                 const bool &useAcceleration = true,
                 const bool &debugSkeleton = false,
                 const std::string debuggingPrefix = NONE);

    ~Skeletonizer() { }

    /**
     * @brief initialize
     * @param verbose
     */
    void initialize(const bool verbose = VERBOSE);

    /**
     * @brief skeletonizeVolumeToCenterLines
     * @param verbose
     */
    virtual void skeletonizeVolumeToCenterLines(const bool verbose = VERBOSE);

    /**
     * @brief skeletonizeVolumeBlockByBlock
     * @param blockSize
     * @param numberOverlappingVoxels
     * @param numberZeroVoxels
     */
    virtual void skeletonizeVolumeBlockByBlock(const size_t& blockSize = 512,
                                               const size_t& numberOverlappingVoxels = 25,
                                               const size_t& numberZeroVoxels = 5,
                                               const bool verbose = VERBOSE);



    std::vector< Vector3f > getShellPoints();

    Sections getValidSections() const;





    virtual void constructGraph(const bool verbose = VERBOSE);



    virtual void segmentComponents(const bool verbose = VERBOSE) = 0;


    void applyVolumeThinningWithDomainDecomposition();

    void thinVolumeBlockByBlock(const size_t& blockSize = 512,
                                const size_t& numberOverlappingVoxels = 25,
                                const size_t& numberZeroVoxels = 5);

    void applyVolumeThinningToVolume(Volume* volume, const bool &displayProgress = true);

protected:

    /**
     * @brief _computeVolumeFromMesh
     * If the input datasets to the constructor is a mesh, compute the volume such that we can
     * proceed with the Skeletonization.
     */
    void _computeVolumeFromMesh();

    /**
     * @brief _computeShellPoints
     * Computes the shell points of the volume, i.e. the points that are representing the surface
     * shell of the volume, or those on its boundary and not internal.
     * @param verbose
     */
    void _computeShellPoints(const bool verbose = VERBOSE);

    /**
     * @brief _computeShellPointsUsingAcceleration
     * Computes the shell points of the volume, i.e. the points that are representing the surface
     * shell of the volume, or those on its boundary and not internal. This method uses the
     * acceleration data structures to improve the performance, but indeed requires more memory.
     * @param thinningVoxels
     * The ThinningVoxels list that contains all the filled voxels that will be used directly
     * to skeletonize the volume.
     * @param verbose
     */
    void _computeShellPointsUsingAcceleration(ThinningVoxelsUI16List &thinningVoxels,
                                              const bool verbose = VERBOSE);

    /**
     * @brief _scaleShellPoints
     * Scales the extracted shell points to ensure their mapping to the input mesh, if any.
     * @param verbose
     */
    void _scaleShellPoints(const bool verbose = VERBOSE);

    /**
     * @brief _applyVolumeThinning
     * Apply the VolumeThinning kernels to extract the center-lines of the volume.
     * @param verbose
     */
    void _applyVolumeThinning(const bool verbose = VERBOSE);

    /**
     * @brief _applyVolumeThinningUsingAcceleration
     * Apply the VolumeThinning kernels to extract the center-lines of the volume, using the
     * acceleration data structures that were built during the initialization stage.
     */

    /**
     * @brief _applyVolumeThinningUsingAcceleration
     * Apply the VolumeThinning kernels to extract the center-lines of the volume, using the
     * acceleration data structures that were built during the initialization stage.
     * @param verbose
     */
    void _applyVolumeThinningUsingAcceleration(const bool verbose = VERBOSE);

    /**
     * @brief _extractNodesFromVoxels
     * Extract the graph nodes from the centerline voxels remaining in the volume.
     * @param verbose
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxels(const bool verbose = VERBOSE);

    /**
     * @brief _extractNodesFromVoxelsNaive
     * Extract the graph nodes from the centerline voxels remaining in the volume using the
     * naive loops. Note that this function is just keep for benchmarking.
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsNaive(const bool verbose = VERBOSE);

    /**
     * @brief _extractNodesFromVoxelsUsingSlicing
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsUsingSlicing(const bool verbose = VERBOSE);

    /**
     * @brief _extractNodesFromVoxelsUsingAcceleration
     * @return
     */
    std::map< size_t, size_t > _extractNodesFromVoxelsUsingAcceleration(const bool verbose = VERBOSE);

    SkeletonBranch* _buildBranch(SkeletonNode* firstNode, SkeletonNode* edgeNode);

    void _buildAcyclicTree(SkeletonBranch* branch, SkeletonBranches &branches);


    /**
     * @brief _inflateNodes
     */
    void _inflateNodes(const bool verbose = VERBOSE);

    /**
     * @brief _inflateNodesUsingAcceleration
     */
    void _inflateNodesUsingAcceleration(const bool verbose = VERBOSE);

    /**
     * @brief _inflateNodesUsingAcceleration
     */
    void _inflateNodesNatively(const bool verbose = VERBOSE);


    void _connectNodesToBuildEdges(const std::map< size_t, size_t > &indicesMapper,
                                   const bool verbose = VERBOSE);

    void _removeTriangleLoops(const bool verbose = VERBOSE);

    void _buildBranchesFromNodes(const SkeletonNodes& nodes);

    /**
     * @brief _exportGraphNodes
     * Debugging function.
     * @param prefix
     * @param verbose
     */
    void _exportGraphNodes(const std::string prefix, const bool verbose = VERBOSE);

    /**
     * @brief _verifySkeletonNodes
     */
    void _verifySkeletonNodes();

    /**
     * @brief _verifySkeletonEdges
     */
    void _verifySkeletonEdges();


protected:

    /**
     * @brief _volume
     * Input volume that will be used to extract the skeleton.
     */
    Volume* _volume = nullptr;

    /**
     * @brief _mesh
     */
    Mesh* _mesh = nullptr;

    /**
     * @brief _voxelizationOptions
     */
    const VoxelizationOptions _voxelizationOptions;

    /**
     * @brief _useAcceleration
     * Use acceleration data structure to improve the performance of the skeletonizer.
     */
    const bool _useAcceleration;

    /**
     * @brief _debugSkeleton
     * Generate as many debugging files as possible to be able to debug the skeleton construction
     * in case of failure.
     */
    const bool _debugSkeleton;

    /**
     * @brief _debuggingPrefix
     * The preix, or path, where the debugging files will be writte.
     */
    const std::string _debuggingPrefix;

    /**
     * @brief _debug
     */
    const bool _debug;

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
     * @brief _edges
     */
    SkeletonEdges _edges;

    /**
     * @brief _branches
     * The branches of the skeleton representing the reconstructed graph.
     */
    SkeletonBranches _branches;

    /**
     * @brief _totalNumberNodes
     * Total number of nodes in the detected skeleton.
     */
    size_t _totalNumberNodes;

    /**
     * @brief _totalNumberEdges
     * Total number of edges in the detected skeleton.
     */
    size_t _totalNumberEdges;
};

}

