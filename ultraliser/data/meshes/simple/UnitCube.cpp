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

#include "UnitCube.h"

namespace Ultraliser
{

UnitCube::UnitCube(const float scale)
    : _scale(scale)
{
    // Build the vertices
    _vertices[0] = Vertex(0, 0, 0) * _scale;
    _vertices[1] = Vertex(1, 0, 0) * _scale;
    _vertices[2] = Vertex(1, 1, 0) * _scale;
    _vertices[3] = Vertex(0, 1, 0) * _scale;
    _vertices[4] = Vertex(0, 0, 1) * _scale;
    _vertices[5] = Vertex(1, 0, 1) * _scale;
    _vertices[6] = Vertex(1, 1, 1) * _scale;
    _vertices[7] = Vertex(0, 1, 1) * _scale;

    // Build the triangles
    _triangles[0] = Triangle(0, 3, 1);
    _triangles[0] = Triangle(1, 3, 2);
    _triangles[0] = Triangle(5, 4, 0);
    _triangles[0] = Triangle(5, 0, 1);
    _triangles[0] = Triangle(6, 5, 1);
    _triangles[0] = Triangle(1, 2, 6);
    _triangles[0] = Triangle(3, 6, 2);
    _triangles[0] = Triangle(3, 7, 6);
    _triangles[0] = Triangle(4, 3, 0);
    _triangles[0] = Triangle(4, 7, 3);
    _triangles[0] = Triangle(7, 4, 5);
    _triangles[0] = Triangle(7, 5, 6);
}

Vertex* UnitCube::getVertices()
{
    return _vertices;
}

Triangle* UnitCube::getTriangles()
{
    return _triangles;
}

}
