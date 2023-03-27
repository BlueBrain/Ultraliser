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

#include "SparseOctreeGridTraverser.h"

namespace Ultraliser
{
SparseOctreeGridTraverser::SparseOctreeGridTraverser(const Vec3ui_32 &gridPoint, uint32_t gridResolution)
    : _gridPoint(gridPoint)
    , _min(0)
    , _max(gridResolution)
{
}

SparseOctreeGridTraverser::VoxelMask SparseOctreeGridTraverser::next()
{
    auto gridDims = _max - _min;
    auto nodeCenter = _min + Vec3ui_32(gridDims.x() >> 1, gridDims.y() >> 1, gridDims.z() >> 1);

    uint8_t xShift = _gridPoint.x() < nodeCenter.x() ? 0 : 1;
    uint8_t yShift = _gridPoint.y() < nodeCenter.y() ? 0 : 2;
    uint8_t zShift = _gridPoint.z() < nodeCenter.z() ? 0 : 4;

    _min[0] = xShift == 0 ? _min.x() : nodeCenter.x();
    _max[0] = xShift == 0 ? nodeCenter.x() : _max.x();
    _min[1] = yShift == 0 ? _min.y() : nodeCenter.y();
    _max[1] = yShift == 0 ? nodeCenter.y() : _max.y();
    _min[2] = zShift == 0 ? _min.z() : nodeCenter.z();
    _max[2] = zShift == 0 ? nodeCenter.z() : _max.z();

    return VoxelMask(xShift, yShift, zShift);
}
}