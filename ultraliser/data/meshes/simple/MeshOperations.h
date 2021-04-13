/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_MESH_OPERATIONS_H
#define ULTRALISER_DATA_MESH_SIMPLE_MESH_OPERATIONS_H

#include <common/Common.h>
#include <data/common/CommonData.h>
#include <data/meshes/simple/primitives/Primitives.h>

namespace Ultraliser
{

/**
 * @brief computeNormal
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
Vector3f computeNormal(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeMeshBoundingBox
 * @param vertices
 * @param numberVertices
 * @param pMin
 * @param pMax
 * @return
 */
void computeMeshBoundingBox(const Vertex* vertices,
                            const int64_t &numberVertices,
                            Vector3f& pMin, Vector3f& pMax,
                            const bool &verbose = true);

/**
 * @brief computeMeshVolume
 * @param vertices
 * @param triangles
 * @param numberTriangles
 * @return
 */
float computeMeshVolume(const Vertex* vertices,
                        const Triangle* triangles,
                        const int64_t &numberTriangles,
                        const bool& verbose = true);

/**
 * @brief computeMeshSurfaceArea
 * @param vertices
 * @param triangles
 * @param numberTriangles
 * @return
 */
float computeMeshSurfaceArea(const Vertex* vertices,
                             const Triangle* triangles,
                             const int64_t &numberTriangles,
                             const bool& verbose = true);

/**
 * @brief importOBJ
 * @param file
 * @param vertices
 * @param triangles
 */
void importOBJ(const std::string &filePath, Vertices& vertices, Triangles& triangles);

/**
 * @brief importOBJ
 * @param file
 * @param vertices
 * @param triangles
 */
void importPLY(const std::string &filePath, Vertices& vertices, Triangles& triangles);

/**
 * @brief importSTL
 * @param filePath
 * @param vertices
 * @param triangles
 * @param verbose
 */
void importSTL(const std::string &filePath, Vertices& vertices, Triangles& triangles);

/**
 * @brief importOFF
 * @param filePath
 * @param vertices
 * @param triangles
 */
void importOFF(const std::string &filePath, Vertices& vertices, Triangles& triangles);

/**
 * @brief exportToOBJ
 * @param prefix
 * @param vertices
 * @param numberVertices
 * @param triangles
 * @param numberTriangles
 */
void exportOBJ(const std::string &prefix,
                 const Vertex* vertices,
                 const uint64_t &numberVertices,
                 const Triangle* triangles,
                 const uint64_t &numberTriangles);

/**
 * @brief exportToOFF
 * @param prefix
 * @param vertices
 * @param numberVertices
 * @param triangles
 * @param numberTriangles
 */
void exportOFF(const std::string &prefix,
                 const Vertex* vertices,
                 const uint64_t &numberVertices,
                 const Triangle* triangles,
                 const uint64_t &numberTriangles);

/**
 * @brief exportToSTL
 * @param prefix
 * @param vertices
 * @param numberVertices
 * @param triangles
 * @param numberTriangles
 */
void exportSTL(const std::string &prefix,
                 const Vertex *vertices,
                 const uint64_t &numberVertices,
                 const Triangle* triangles,
                 const uint64_t &numberTriangles);

/**
 * @brief exportPLY
 * @param prefix
 * @param vertices
 * @param numberVertices
 * @param triangles
 * @param numberTriangles
 * @param writeASCII
 */
void exportPLY(const std::string &prefix,
               const Vertex *vertices,
               const uint64_t &numberVertices,
               const Triangle* triangles,
               const uint64_t &numberTriangles,
               bool writeASCII = false);


}

#endif // ULTRALISER_DATA_MESH_SIMPLE_MESH_OPERATIONS_H
