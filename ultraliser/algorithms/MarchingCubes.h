/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
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
 **************************************************************************************************/

#pragma once

#include <algorithms/DualMarchingCubes.hh>
#include <data/volumes/Volume.h>
#include <data/volumes/TaggedVolume.h>
#include <data/meshes/advanced/AdvancedMesh.h>

namespace Ultraliser
{

/**
 * @brief The MarchingCubes class
 */
class MarchingCubes
{
public:

    /**
     * @brief MarchingCubes
     *
     * @param volume
     * @param isoValue
     */
    MarchingCubes(Volume* volume, const size_t isoValue = 127);

    /**
     * @brief generateMesh
     * Generate the mesh.
     * @return
     * The reconstructed mesh.
     */
    Mesh* generateMesh();

    /**
     * @brief generateAdvancedMesh
     * Generate the mesh.
     * @return
     * The reconstructed mesh.
     */
    AdvancedMesh* generateAdvancedMesh();

    /**
     * @brief generateMeshFromVolume
     * Generate a mesh from the MC algorithm given an input volume.
     *
     * @param volume
     * An input volume that will be used to create the mesh.
     * @return
     * A pointer to the mesh.
     */
    static Mesh* generateMeshFromVolume(Volume *volume);

    /**
     * @brief generateAdvancedMeshFromVolume
     * Generate an advanced mesh from the MC algorithm given an input volume.
     *
     * @param volume
     * An input volume that will be used to create the mesh.
     * @return
     * A pointer to the mesh.
     */
    static AdvancedMesh* generateAdvancedMeshFromVolume(Volume *volume);

private:

    /**
     * @brief _buildSharedVertices
     * Extract triangular mesh with shared vertex indices, in parallel using OpenMP.
     *
     * @param vertices
     * A list to collect the vertices of the mesh.
     * @param triangles
     * A list to collect the triangles of the mesh.
     */
    void _buildSharedVertices(Vertices& vertices, Triangles &triangles);

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
     * @brief _meshExtractionTime
     */
    double _meshExtractionTime;
};

}
