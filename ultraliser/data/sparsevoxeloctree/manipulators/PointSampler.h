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
#include <math/Vector.h>

namespace Ultraliser
{
class PointSampler
{
public:
    explicit PointSampler(const SparseOctree &octree);
    uint32_t getResolution() const;
    bool sample(const Vector3f &point) const;
    bool sample(const Vec3ui_32 &gridCoordinates) const;

private:
    const SparseOctreeNode &_root;
    Vector3f _boundMin;
    Vector3f _oneOverBoundsDims;
    uint8_t _maxDepth;
    uint32_t _gridSize;
};
}