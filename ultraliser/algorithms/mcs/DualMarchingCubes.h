/***************************************************************************************************
 * Copyright (c) 2018 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Dominik Wodniok
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
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
 *
 * This file has been adapted from the original implementation of Dominik Wodniok who implemented a
 * serial version of the Dual Marching Cubes algorithm from Gregory M. Nielson. In this file, the
 * original implementation is extended by Marawn Abdellah to make a parallel implementation of the
 * algorithm. The original code of dualmc is available under the BSD 3-Clause License as indicated
 * at < https://github.com/dominikwodniok/dualmc >.
 **************************************************************************************************/

#pragma once

#include <algorithms/mcs/DualMarchingCubes.hh>
#include <data/volumes/Volume.h>
#include <data/volumes/TaggedVolume.h>
#include <data/meshes/advanced/AdvancedMesh.h>

// Default iso surface value
#define DEFAULT_ISO_VALUE 127

namespace Ultraliser
{

class DualMarchingCubes
{
public:

    /**
     * @brief DualMarchingCubes
     * @param volume
     * @param isoValue
     */
    DualMarchingCubes(Volume* volume,
                      const size_t isoValue = DEFAULT_ISO_VALUE,
                      const bool &generateManifold = true);

    /**
     * @brief DualMarchingCubes
     * @param volume
     * @param isoValue
     * @param generateManifold
     */
    DualMarchingCubes(TaggedVolume* volume,
                      const uint8_t isoValue = DEFAULT_ISO_VALUE,
                      const bool &generateManifold = true);

    /**
     * @brief generateMesh
     * @param paralle
     * @return
     */
    Mesh* generateMesh();

    /**
     * @brief generateAdvancedMesh
     * @param paralle
     * @return
     */
    AdvancedMesh* generateAdvancedMesh();

public:

    /**
     * @brief generateMeshFromVolume
     * Generate a mesh from the DMC algorithm given an input volume.
     * @param volume
     * An input volume that will be used to create the mesh.
     * @return
     * A pointer to the mesh.
     */
    static Mesh* generateMeshFromVolume(Volume *volume);

    /**
     * @brief generateMeshFromVolume
     * Generate an avanced mesh from the DMC algorithm given an input volume.
     * @param volume
     * An input volume that will be used to create the mesh.
     * @return
     * A pointer to the resulting mesh.
     */
    static AdvancedMesh* generateAdvancedMeshFromVolume(Volume* volume);

private:

    /**
     * @brief _index
     * @param x
     * @param y
     * @param z
     * @return
     */
    int64_t _index(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief _buildSharedVertices
     * Extract quad mesh with shared vertex indices, but in parallel using all
     * the CPUs available.
     * @param mesh
     */
    void _buildSharedVertices(Vertices& vertices, Triangles &triangles);

    /**
     * @brief _getCellCode
     * Get the 8-bit in-out mask for the voxel corners of the cell cube at
     * a specific (x, y, z).
     * @param x
     * @param y
     * @param z
     * @return
     */
    int _getCellCode(const int64_t &x, const int64_t &y, const int64_t &z) const;

    /**
     * @brief _getDualPointCode
     * Get the 12-bit dual point code mask, which encodes the traditional
     * marching cube vertices of the traditional marching cubes face which
     * corresponds to the dual point.
     * This is also where the advanced dual marching cubes algorithm is
     * implemented.
     * @param x
     * @param y
     * @param z
     * @param edge
     * @return
     */
    int _getDualPointCode(const int64_t x, const int64_t y, const int64_t z,
                          const DMC_EDGE_CODE edge) const;

    /**
     * @brief _calculateDualPoint
     * Given a dual point code and iso value, compute the dual point.
     * @param x
     * @param y
     * @param z
     * @param pointCode
     * @param v
     */
    void _calculateDualPoint(const int64_t &x, const int64_t &y, const int64_t &z,
                             const int32_t &pointCode, Vector3f &v) const;

    /**
     * @brief _getSharedDualPointIndex
     * Get the shared index of a dual point which is uniquly identified by its
     * cell cube index and a cube edge. The dual point is computed, if it has
     * not been computed before.
     * @param x
     * @param y
     * @param z
     * @param edge
     * @param vertices
     * @return
     */
    int64_t _getSharedDualPointIndex(const int64_t &x, const int64_t &y, const int64_t &z,
                                     const DMC_EDGE_CODE &edge, std::vector<Vector3f> &vertices);

private:

    /**
     * @brief pointToIndex
     * Hash map for shared vertex index computations
     */
    std::unordered_map< DualPointKey, int64_t, DualPointKeyHash > pointToIndex;

private:

    /**
     * @brief _volume
     */
    const Volume* _volume;

    /**
     * @brief _isoValue
     */
    const size_t _isoValue;

    /**
     * @brief _generateManifold
     */
    const bool _generateManifold;

    /**
     * @brief _dmcGenerationTime
     */
    double _meshExtractionTime;
};

}
