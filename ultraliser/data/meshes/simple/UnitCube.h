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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_UNIT_CUBE_H
#define ULTRALISER_DATA_MESH_SIMPLE_UNIT_CUBE_H

#include <data/meshes/simple/primitives/Primitives.h>


namespace Ultraliser
{

/**
 * @brief The UnitCube class
 */
class UnitCube
{
public:

    /**
     * @brief UnitCube
     * @param scale
     */
    UnitCube(const float scale = 1.f);

    /**
     * @brief getVertices
     * @return
     */
    Vertex* getVertices();

    /**
     * @brief getTriangles
     * @return
     */
    Triangle* getTriangles();

private:

    /**
     * @brief _scale
     */
    const float _scale;

    /**
     * @brief _vertices
     */
    Vertex _vertices[8];

    /**
     * @brief _triangles
     */
    Triangle _triangles[12];
};

}

#endif // ULTRALISER_DATA_MESH_SIMPLE_UNIT_CUBE_H