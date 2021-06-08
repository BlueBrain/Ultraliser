#ifndef MARCHINGCUBES_H
#define MARCHINGCUBES_H

#include <algorithms/DualMarchingCubes.hh>
#include <data/volumes/volumes/Volume.h>
#include <data/volumes/volumes/TaggedVolume.h>
#include <data/meshes/advanced/AdvancedMesh.h>

namespace Ultraliser
{

class MarchingCubes
{
public:
    MarchingCubes(Volume* volume,
                  const uint8_t isoValue = 127);

    Mesh* generateMesh();

    Mesh* marching_cubes();

    Mesh* polygonise();

    /**
     * @brief generateMeshFromVolume
     * Generate a mesh from the DMC algorithm given an input volume.
     * @param volume
     * An input volume that will be used to create the mesh.
     * @return
     * A pointer to the mesh.
     */
    static Mesh* generateMeshFromVolume(Volume *volume);
private:

    /**
     * @brief _volume
     */
    const Volume* _volume;

    /**
     * @brief _isoValue
     */
    const uint8_t _isoValue;
};

}

#endif // MARCHINGCUBES_H
