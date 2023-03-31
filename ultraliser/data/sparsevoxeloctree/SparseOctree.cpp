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

#include "SparseOctree.h"

#include "common/SparseOctreeNodeSlot.h"

namespace Ultraliser
{
SparseOctree::SparseOctree(const Bounds &bounds, uint8_t maxDepth)
    : _bounds(bounds)
    , _maxDepth(maxDepth)
    , _root(SparseOctreeNodeSlot::root)
{
}

SparseOctreeNode &SparseOctree::getRoot()
{
    return _root;
}

const SparseOctreeNode &SparseOctree::getRoot() const
{
    return _root;
}

uint8_t SparseOctree::getMaxDepth() const
{
    return _maxDepth;
}

const Bounds &SparseOctree::getBounds() const
{
    return _bounds;
}
/*
SparseOctree::VolumeView::VolumeView(SparseOctree &octree)
    : _octree(octree)
    , _sampler(octree)
    , _voxelSize(octree.getBounds().getDimensions() / static_cast<float>(1 << octree.getMaxDepth()))
{
}

int64_t SparseOctree::VolumeView::getWidth() const
{
    return static_cast<int64_t>(1 << _octree.getMaxDepth());
}

int64_t SparseOctree::VolumeView::getHeight() const
{
    return getWidth();
}

int64_t SparseOctree::VolumeView::getDepth() const
{
    return getWidth();
}

uint64_t SparseOctree::VolumeView::getValueUI64(int64_t x, int64_t y, int64_t z) const
{
    assert(x > -1 && x < getWidth());
    assert(y > -1 && y < getHeight());
    assert(z > -1 && z < getDepth());

    auto point = Vec3ui_32(I2UI32(x), I2UI32(y), I2UI32(z));
    return _sampler.sample(point) ? 1 : 0;
}

void SparseOctree::VolumeView::getVoxelBoundingBox(int64_t x, int64_t y, int64_t z, Vector3f &min, Vector3f &max) const
{
    assert(x > -1 && x < getWidth());
    assert(y > -1 && y < getHeight());
    assert(z > -1 && z < getDepth());

    auto oneOverResolution = 1.f / getWidth();
    auto delta = Vector3f(x * oneOverResolution, y * oneOverResolution, z * oneOverResolution);

    auto &bounds = _octree.getBounds();
    auto dimensions = bounds.getDimensions();
    auto offset = dimensions * delta;

    min = bounds.getMin() + offset;
    max = min + _voxelSize;
}
*/
}
