/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Nadir Roman Guerrero < nadir.romanguerrero@epfl.ch >
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

#include <data/sparsevoxeloctree/SparseOctree.h>

#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{
class OctreeTriangleVoxelizer
{
public:
    explicit OctreeTriangleVoxelizer(SparseOctree &octree);

    void voxelize(const Vector3f &a, const Vector3f &b, const Vector3f &c);

private:
    SparseOctreeNode &_root;
    const Ultraliser::Bounds &_rootBounds;
    uint8_t _maxDepth;
};

class OctreeMeshVoxelizer
{
public:
    static void voxelize(SparseOctree &octree, const Mesh &mesh);
};
}