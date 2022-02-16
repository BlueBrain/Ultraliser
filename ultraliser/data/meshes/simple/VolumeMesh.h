/***************************************************************************************************
 * Copyright (c) 2021
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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_VOLUME_MESH_H
#define ULTRALISER_DATA_MESH_SIMPLE_VOLUME_MESH_H

#include <common/Common.h>
#include <math/Math.h>
#include <data/meshes/simple/primitives/Primitives.h>

namespace Ultraliser
{

class VolumeMesh
{
public:
    /**
     * @brief VolumeMesh
     */
    VolumeMesh() { /* EMPTY CONSTRUCTOR */ }

    /**
     * @brief VolumeMesh
     * @param numberVertices
     * @param numberTriangles
     */
    VolumeMesh(const uint64_t numberVertices, const uint64_t numberTriangles);

    /**
     * @brief append
     * @param inputMesh
     */
    void append(const VolumeMesh *inputMesh);

    /**
     * @brief constructUnitCube
     * @param scale
     * @return
     */
    static VolumeMesh* constructUnitCube(const float scale = 1.f);

    /**
     * @brief constructVoxelCube
     * @param pMin
     * @param pMax
     * @return
     */
    static VolumeMesh* constructVoxelCube(const Vector3f &pMin, const Vector3f &pMax);

    /**
     * @brief translate
     * @param to
     */
    void translate(const Vector3f& to);

    /**
     * @brief rotate
     * @param matrix
     */
    void rotate(const Matrix4f &matrix);

    /**
     * @brief transform
     * @param matrix
     */
    void transform(const Matrix4f &matrix);

    /**
     * @brief uniformScale
     * @param factor
     */
    void uniformScale(const float factor);

    /**
     * @brief scale
     * @param x
     * @param y
     * @param z
     */
    void scale(const float x = 1.0, const float y = 1.0, const float z = 1.0);

    /**
     * @brief scale
     * @param factor
     */
    void scale(const Vector3f& factor);

    /**
     * @brief computeBoundingBox
     * @param pMin
     * @param pMax
     */
    void computeBoundingBox(Vector3f& pMin, Vector3f& pMax);

    /**
     * @brief centerAtOrigin
     */
    void centerAtOrigin(void);

    ~VolumeMesh();

public:

    /**
     * @brief _vertices
     */
    std::vector< Vector3f > vertices;

    /**
     * @brief _triangles
     */
    std::vector< Triangle > triangles;
};

}

#endif // ULTRALISER_DATA_MESH_SIMPLE_VOLUME_MESH_H
