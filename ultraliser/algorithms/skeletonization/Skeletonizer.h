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
     * @param mesh
     */
    Skeletonizer(Volume *volume, const Mesh *mesh);

    virtual void skeletonizeVolume() = 0;

    virtual void skeletonizeVolumeBlockByBlock(const size_t& blockSize = 512,
                                       const size_t& numberOverlappingVoxels = 25,
                                       const size_t& numberZeroVoxels = 5) = 0;


    /**
     * @brief initialize
     * Initialize the skeletonizer.
     */
    void initialize();

    std::vector< Vector3f > getShellPoints();

    void applyVolumeThinning();
    void applyVolumeThinningUsingThinningVoxels();


    virtual void constructGraph();



    virtual void segmentComponents() = 0;


    void applyVolumeThinningWithDomainDecomposition();

    void thinVolumeBlockByBlock(const size_t& blockSize = 512,
                                const size_t& numberOverlappingVoxels = 25,
                                const size_t& numberZeroVoxels = 5);

    void applyVolumeThinningToVolume(Volume* volume, const bool &displayProgress = true);

private:


protected:


    void _computeShellPoints();
    void _computeShellPointsWithThinningVoxels(ThinningVoxelsUI16 &thinningVoxels);


    void _voxelizeMesh();




    SkeletonBranch* _buildBranch(SkeletonNode* firstNode, SkeletonNode* edgeNode);

    void _buildAcyclicTree(SkeletonBranch* branch, SkeletonBranches &branches);


    std::map< size_t, size_t > _extractNodesFromVoxels();

    std::map< size_t, size_t > _extractNodesFromVoxelsParallel();

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
     * @brief _mesh
     * Optional input mesh that will be used to compute the actual coordinates of the resulting
     * skeleton. If this input mesh is empty, i.e. a nullptr, the dimensions will be computed based
     * on the dimensions of the volume.
     */
    const Mesh* _mesh;

    /// Mesh bounding box
    Vector3f _pMinMesh, _pMaxMesh, _centerMesh, _boundsMesh;

    /// Volume bounding box, will be used if no input mesh is used.
    Vector3f _pMinVolume, _pMaxVolume, _centerVolume, _boundsVolume;

    /// Mesh to volume scale factor
    Vector3f _scaleFactor;

    /**
     * @brief _shellPoints
     * A list of the poins (or voxels) that represent the shell of the given volume.
     */
    std::vector< Vector3f > _shellPoints;

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


    bool _useThinningVoxels;

};
}

