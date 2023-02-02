#pragma once

#include <data/meshes/simple/Mesh.h>
#include <data/volumes/Volume.h>
#include <algorithms/skeletonization/SkeletonNode.hh>
#include <algorithms/skeletonization/SkeletonEdge.hh>
#include <algorithms/skeletonization/SkeletonBranch.hh>
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
    Skeletonizer(const Mesh* mesh, Volume *volume);

    std::vector< Vector3f > getShellPoints();

    void applyVolumeThinning();

    SkeletonNodes constructGraph();

    Mesh* getSomaMesh() const
    {
        return _somaMesh;
    }

    void segmentComponents(SkeletonNodes& nodes);

private:

    void _computeShellPoints();

    void _voxelizeMesh();


    SkeletonBranches _buildBranchesFromNodes(const SkeletonNodes& nodes);

    SkeletonBranch* _buildBranch(SkeletonNode* firstNode, SkeletonNode* edgeNode);

    Mesh* _reconstructSoma(const SkeletonBranches &branches);

    void _buildAcyclicTree(SkeletonBranch* branch, SkeletonBranches &branches);



private:

    /**
     * @brief _mesh
     * The input mesh to the Skeletonizer.
     */
    const Mesh* _mesh;

    /**
     * @brief _volume
     * The volume that contains a solid voxel grid of the input mesh.
     */
    Volume* _volume;

    /// Mesh bounding box
    Vector3f _pMinMesh, _pMaxMesh, _centerMesh, _boundsMesh;

    /// Volume bounding box
    Vector3f _pMinVolume, _pMaxVolume, _centerVolume, _boundsVolume;

    /// Mesh to volume scale factor
    Vector3f _scaleFactor;

    std::vector< Vector3f > _shellPoints;

    Mesh* _somaMesh = nullptr;


};
}

