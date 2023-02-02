#pragma once

#include <data/meshes/simple/Mesh.h>
#include <data/volumes/Volume.h>

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

private:

    void _computeShellPoints();

    void _voxelizeMesh();

    void _constructGraph();


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


};
}

