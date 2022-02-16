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

#include "DualMarchingCubes.h"
#include <data/volumes/voxels/DMCVoxel.h>
#include <data/meshes/simple/primitives/Primitives.h>

namespace Ultraliser
{

DualMarchingCubes::DualMarchingCubes(Volume *volume,
                                     const uint8_t isoValue,
                                     const bool& generateManifold)
    : _volume(volume)
    , _isoValue(isoValue)
    , _generateManifold(generateManifold)
{
    /// EMPTY CONSTRUCTOR
}

Mesh* DualMarchingCubes::generateMesh(const bool& paralle)
{
    LOG_TITLE("Mesh Reconstruction with DMC");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Strat the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    if (paralle)
        _buildSharedVerticesParallel(vertices, triangles);
    else
        _buildSharedVertices(vertices, triangles);

    Mesh* mesh = new Mesh(vertices, triangles);

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Dual Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

AdvancedMesh* DualMarchingCubes::generateAdvancedMesh(const bool& paralle)
{
    LOG_TITLE("Mesh Reconstruction with DMC");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Strat the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    if (paralle)
        _buildSharedVerticesParallel(vertices, triangles);
    else
        _buildSharedVertices(vertices, triangles);

    AdvancedMesh* mesh = new AdvancedMesh(vertices, triangles);

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Dual Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

int DualMarchingCubes::_getCellCode(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    // Determine for each cube corner if it is outside or inside
    int code = 0;

    if (_volume->getValue(x, y, z) >= _isoValue)
        code |= 1;
    if (_volume->getValue(x + 1, y, z) >= _isoValue)
        code |= 2;
    if (_volume->getValue(x, y + 1, z) >= _isoValue)
        code |= 4;
    if (_volume->getValue(x + 1, y + 1, z) >= _isoValue)
        code |= 8;
    if (_volume->getValue(x, y, z + 1) >= _isoValue)
        code |= 16;
    if (_volume->getValue(x + 1, y, z + 1) >= _isoValue)
        code |= 32;
    if (_volume->getValue(x, y + 1, z + 1) >= _isoValue)
        code |= 64;
    if (_volume->getValue(x + 1, y + 1, z + 1) >= _isoValue)
        code |= 128;

    return code;
}

int64_t DualMarchingCubes::_index(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    return x + _volume->getWidth() * (y + _volume->getHeight() * z);
}

int DualMarchingCubes::_getDualPointCode(const int64_t x, const int64_t y, const int64_t z,
                                         const DMC_EDGE_CODE edge) const
{
    // Get the code of the cube that corresponds to the given XYZ voxel
    int cubeCode = _getCellCode(x, y, z);

    // Is advanced dual marching cubes desired?
    if (_generateManifold)
    {
        // The Manifold Dual Marching Cubes approach from Rephael Wenger as
        // described in chapter 3.3.5 of his book "Isosurfaces: Geometry,
        // Topology, and Algorithms" is implemente here.
        // If a problematic C16 or C19 configuration shares the ambiguous face
        // with another C16 or C19 configuration we simply invert the cube code
        // before looking up dual points.
        // Doing this for these pairs ensures advanced meshes.
        // But this removes the dualism to marching cubes.

        // Check if we have a potentially problematic configuration
        const uint8_t direction = PROBLEMATIC_CONFIGS[I2UI16(cubeCode)];

        // If the direction code is in {0,...,5} we have a C16 or C19
        // configuration.
        if (direction != 255)
        {
            // We have to check the neighboring cube, which shares the ambiguous
            // face. For this we decode the direction. This could also be done
            // with another lookup table.
            // Copy current cube coordinates into an array.
            int64_t neighborCoords[] = { x, y, z };

            // Get the dimension of the non-zero coordinate axis
            unsigned int const component = direction >> 1;

            // Get the sign of the direction
            int64_t delta = (direction & 1) == 1 ? 1 : -1;

            // Modify the correspong cube coordinate
            neighborCoords[component] += delta;

            int64_t dimension;
            if (component == 0)
                dimension = _volume->getWidth();
            else if (component == 1)
                dimension = _volume->getHeight();
            else
                dimension = _volume->getWidth();

            // Have we left the volume in this direction?
            if (neighborCoords[component] >= 0 &&
                neighborCoords[component] < (dimension - 1))
            {
                // Get the cube configuration of the relevant neighbor
                int neighborCubeCode = _getCellCode(neighborCoords[0],
                                                    neighborCoords[1],
                                                    neighborCoords[2]);

                // Look up the neighbor configuration ambiguous face direction.
                // If the direction is valid we have a C16 or C19 neighbor.
                // As C16 and C19 have exactly one ambiguous face this face is
                // guaranteed to be shared for the pair.
                if (PROBLEMATIC_CONFIGS[I2UI16(neighborCubeCode)] != 255)
                {
                    // Replace the cube configuration with its inverse.
                    cubeCode ^= 0xff;
                }
            }
        }
    }

    for (int32_t i = 0; i < 4; ++i)
    {
        if (I2UI32(DUAL_POINTS_LIST[cubeCode][i]) & edge)
            return DUAL_POINTS_LIST[cubeCode][i];
    }

    return 0;
}

void DualMarchingCubes::_calculateDualPoint(const int64_t &x, const int64_t &y, const int64_t &z,
                                            const int32_t& pointCode, Vector3f & v) const
{
    // Geth the voxel bounding box
    Vector3f pMin, pMax;
    _volume->getVoxelBoundingBox(x, y, z, pMin, pMax);

    // Initialize the point with higher voxel coordinates
    v = pMax;

    // Compute the dual point as the mean of the face vertices belonging to the
    // original marching cubes face
    Vertex p(0.f);
    int points = 0;

    // Sum edge intersection vertices using the point code
    if (I2UI32(pointCode) & EDGE0)
    {
        p.x() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x + 1, y, z) - _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE1)
    {
        p.x() += 1.0f;
        p.z() += (_isoValue - _volume->getValue(x + 1, y, z)) /
                (_volume->getValue(x + 1, y, z + 1) - _volume->getValue(x + 1, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE2)
    {
        p.x() += (_isoValue - _volume->getValue(x, y, z + 1)) /
                (_volume->getValue(x + 1, y, z + 1) - _volume->getValue(x, y, z + 1));
        p.z() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE3)
    {
        p.z() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x, y, z + 1) - _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE4)
    {
        p.x() += (_isoValue - _volume->getValue(x, y + 1, z)) /
                (_volume->getValue(x + 1, y + 1, z) - _volume->getValue(x, y + 1, z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE5)
    {
        p.x() += 1.0f;
        p.z() += (_isoValue - _volume->getValue(x + 1, y + 1, z)) /
                (_volume->getValue(x + 1, y + 1, z + 1) - _volume->getValue(x + 1, y + 1, z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE6)
    {
        p.x() += (_isoValue - _volume->getValue(x, y + 1, z + 1)) /
                (_volume->getValue(x + 1, y + 1, z + 1) - _volume->getValue(x, y + 1, z + 1));
        p.z() += 1.0f;
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE7)
    {
        p.z() += (_isoValue - _volume->getValue(x, y + 1 , z)) /
                (_volume->getValue(x, y + 1, z + 1) - _volume->getValue(x, y + 1 , z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE8)
    {
        p.y() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x, y + 1, z) - _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE9)
    {
        p.x() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x + 1, y, z)) /
                (_volume->getValue(x + 1, y + 1, z) - _volume->getValue(x + 1, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE10)
    {
        p.x() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x + 1, y, z + 1)) /
                (_volume->getValue(x + 1, y + 1, z + 1) - _volume->getValue(x + 1, y, z + 1));
        p.z() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE11)
    {
        p.z() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x, y, z + 1)) /
                (_volume->getValue(x, y + 1, z + 1) - _volume->getValue(x, y, z + 1));
        points++;
    }

    // Divide by number of accumulated points
    float invPoints = 1.0f / points;
    p.x() *= invPoints; p.y() *= invPoints; p.z() *= invPoints;
}

int64_t DualMarchingCubes::_getSharedDualPointIndex(
        const int64_t &x, const int64_t &y, const int64_t &z,
        const DMC_EDGE_CODE& edge, std::vector< Vector3f >& vertices)
{
    // Create a key for the dual point from its linearized cell ID and
    // point code
    DualPointKey key;
    key.linearizedCellID = _index(x, y, z);
    key.pointCode = _getDualPointCode(x, y, z, edge);

    // Have we already computed the dual point?
    auto iterator = pointToIndex.find(key);
    if (iterator != pointToIndex.end())
    {
        // Just return the dual point index
        return iterator->second;
    }
    else
    {
        // Create new vertex and vertex id
        int64_t newVertexId = I2I64(vertices.size());
        vertices.emplace_back();
        _calculateDualPoint(x, y, z, key.pointCode, vertices.back());

        // Insert vertex ID into map and also return it
        pointToIndex[key] = newVertexId;
        return newVertexId;
    }
}

bool DualPointKey::operator==(DualPointKey const & other) const
{
    return (linearizedCellID == other.linearizedCellID && pointCode == other.pointCode);
}

void DualMarchingCubes::_buildSharedVerticesParallel(Vertices& vertices, Triangles &triangles)
{
    // Start timer
    TIMER_SET;

    // Adding a little bit of extra voxels
    const int64_t extraVoxels = 5;
    const int64_t minValue = -1 * extraVoxels;
    const int64_t maxValue = extraVoxels;

    const int64_t maxX = _volume->getWidth() + maxValue;
    const int64_t maxY = _volume->getHeight() + maxValue;
    const int64_t maxZ = _volume->getDepth() + maxValue;

    const uint64_t sizeX = maxX + extraVoxels;

    // A list of lists of DMCVoxel's
    // This list will have reducedX entries to get filled in parallel
    DMCVoxelsList volumeDMCVoxels;
    volumeDMCVoxels.resize(static_cast< uint64_t >(sizeX));

    // Shared quad points
    int64_t i0, i1, i2, i3;

    // Clear the map before processing any items
    pointToIndex.clear();

    // Searching for non-zero voxels in parallel
    LOOP_STARTS("Searching Filled Voxels");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t x = minValue; x < maxX; ++x)
    {
        // Get a reference to the slice
        DMCVoxels& sliceDMCVoxels = volumeDMCVoxels[static_cast< uint64_t >(x + extraVoxels)];

        for (int64_t y = minValue; y < maxY; ++y)
        {
            for (int64_t z = minValue; z < maxZ; ++z)
            {
                // X-aligned edge
                if (z > minValue && y > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x + 1, y, z);

                    bool const entering = (value0 < _isoValue) && (value1 >= _isoValue);
                    bool const exiting = (value0 >= _isoValue) && (value1 < _isoValue);

                    // Create a DMCVoxel if there is one
                    if (entering || exiting)
                    {
                        sliceDMCVoxels.push_back(
                                    new DMCVoxel(x, y, z, entering, exiting, DMC_EDGE_SIDE::X));
                    }
                }

                // Y-aligned edge
                if (z > minValue && x > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x, y + 1, z);

                    bool const entering = value0 < _isoValue && value1 >= _isoValue;
                    bool const exiting = value0 >= _isoValue && value1 < _isoValue;

                    if (entering || exiting)
                    {
                        sliceDMCVoxels.push_back(
                                    new DMCVoxel(x, y, z, entering, exiting, DMC_EDGE_SIDE::Y));
                    }
                }

                // Z-aligned edge
                if (x > minValue && y > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x, y, z + 1);

                    bool const entering = value0 < _isoValue && value1 >= _isoValue;
                    bool const exiting = value0 >= _isoValue && value1 < _isoValue;

                    if (entering || exiting)
                    {
                        sliceDMCVoxels.push_back(
                                    new DMCVoxel(x, y, z, entering, exiting, DMC_EDGE_SIDE::Z));
                    }
                }
            }
        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, sizeX);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Reset the time
    TIMER_RESET;

    // Building the shared vertices
    LOOP_STARTS("Building Shared Vertices");
    for (uint64_t i = 0; i < volumeDMCVoxels.size(); i++)
    {
        LOOP_PROGRESS(i, volumeDMCVoxels.size());
        for (uint64_t j = 0; j < volumeDMCVoxels[i].size(); j++)
        {
            DMCVoxel* dmcVoxel = volumeDMCVoxels[i][j];

            // X-aligned edge
            if (dmcVoxel->side == DMC_EDGE_SIDE::X)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y, dmcVoxel->z,
                                              EDGE0, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y, dmcVoxel->z - 1,
                                              EDGE2, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y - 1, dmcVoxel->z - 1,
                                              EDGE6, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y - 1, dmcVoxel->z,
                                              EDGE4, vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel->entering)
                {
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]

                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // Quad [i0, i3, i2, i1]
                    // Triangles [i0, i3, i2] [i2, i1, i0]
                    triangles.push_back(Triangle(i0, i3, i2));
                    triangles.push_back(Triangle(i2, i1, i0));
                }
            }

            // Y-aligned edge
            else if (dmcVoxel->side == DMC_EDGE_SIDE::Y)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y, dmcVoxel->z,
                                              EDGE8, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y, dmcVoxel->z - 1,
                                              EDGE11, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel->x - 1, dmcVoxel->y, dmcVoxel->z - 1,
                                              EDGE10, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel->x - 1, dmcVoxel->y, dmcVoxel->z,
                                              EDGE9, vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel->exiting)
                {
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]
                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // Quad [i0, i3, i2, i1]
                    // Triangles [i0, i3, i2] [i2, i1, i0]
                    triangles.push_back(Triangle(i0, i3, i2));
                    triangles.push_back(Triangle(i2, i1, i0));
                }
            }

            // Z-aligned edge
            else if (dmcVoxel->side == DMC_EDGE_SIDE::Z)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y, dmcVoxel->z,
                                              EDGE3, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel->x - 1, dmcVoxel->y, dmcVoxel->z,
                                              EDGE1, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel->x - 1, dmcVoxel->y - 1, dmcVoxel->z,
                                              EDGE5, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel->x, dmcVoxel->y - 1, dmcVoxel->z,
                                              EDGE7, vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel->exiting)
                {
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]
                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // Quad [i0, i3, i2, i1]
                    // Triangles [i0, i3, i2] [i2, i1, i0]
                    triangles.push_back(Triangle(i0, i3, i2));
                    triangles.push_back(Triangle(i2, i1, i0));
                }
            }
            else
            {
                LOG_ERROR("UNKNOWN DMC EDGE!");
            }

            // Delete the DMC voxel
            delete dmcVoxel;
        }

        volumeDMCVoxels[i].clear();
        volumeDMCVoxels[i].shrink_to_fit();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clear the DMC voxel list
    volumeDMCVoxels.clear();
    volumeDMCVoxels.shrink_to_fit();
}

void DualMarchingCubes::_buildSharedVertices(Vertices &vertices, Triangles &triangles)
{
    // Start timer
    TIMER_SET;

    // Adding a little bit of extra voxels
    const int64_t extraVoxels = 5;
    const int64_t minValue = -1 * extraVoxels;
    const int64_t maxValue = extraVoxels;

    const int64_t maxX = _volume->getWidth() + maxValue;
    const int64_t maxY = _volume->getHeight() + maxValue;
    const int64_t maxZ = _volume->getDepth() + maxValue;

    const uint64_t sizeX = maxX + extraVoxels;

    // Quad points
    int64_t i0, i1, i2, i3;

    pointToIndex.clear();

    // Iterate voxels
    LOOP_STARTS("Building Shared Vertices");
    for (int64_t x = minValue; x < maxX; ++x)
    {
        for (int64_t y = minValue; y < maxY; ++y)
        {
            for (int64_t z = minValue; z < maxZ; ++z)
            {

                // Construct quads for X edge
                if (z > minValue && y > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x + 1, y, z);

                    bool const entering = (value0 < _isoValue) && (value1 >= _isoValue);
                    bool const exiting = (value0 >= _isoValue) && (value1 < _isoValue);

                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z, EDGE0, vertices);
                        i1 = _getSharedDualPointIndex(x, y, z - 1, EDGE2, vertices);
                        i2 = _getSharedDualPointIndex(x, y - 1, z - 1, EDGE6, vertices);
                        i3 = _getSharedDualPointIndex(x, y - 1, z, EDGE4, vertices);

                        if (entering)
                        {
                            triangles.push_back(Triangle(i0, i1, i2));
                            triangles.push_back(Triangle(i2, i3, i0));
                        }
                        else
                        {
                            triangles.push_back(Triangle(i0, i3, i2));
                            triangles.push_back(Triangle(i2, i1, i0));
                        }
                    }
                }

                // Construct quads for y edge
                if (z > minValue && x > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x, y + 1, z);

                    bool const entering = value0 < _isoValue && value1 >= _isoValue;
                    bool const exiting = value0 >= _isoValue && value1 < _isoValue;

                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z, EDGE8, vertices);
                        i1 = _getSharedDualPointIndex(x, y, z - 1, EDGE11, vertices);
                        i2 = _getSharedDualPointIndex(x - 1, y, z - 1, EDGE10, vertices);
                        i3 = _getSharedDualPointIndex(x - 1, y, z, EDGE9, vertices);

                        if (exiting)
                        {
                            triangles.push_back(Triangle(i0, i1, i2));
                            triangles.push_back(Triangle(i2, i3, i0));
                        }
                        else
                        {
                            triangles.push_back(Triangle(i0, i3, i2));
                            triangles.push_back(Triangle(i2, i1, i0));
                        }
                    }
                }

                // Construct quads for z edge
                if (x > minValue && y > minValue)
                {
                    const uint64_t value0 = _volume->getValue(x, y, z);
                    const uint64_t value1 = _volume->getValue(x, y, z + 1);

                    bool const entering = value0 < _isoValue && value1 >= _isoValue;
                    bool const exiting = value0 >= _isoValue && value1 < _isoValue;

                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z, EDGE3, vertices);
                        i1 = _getSharedDualPointIndex(x - 1, y, z, EDGE1, vertices);
                        i2 = _getSharedDualPointIndex(x - 1, y - 1, z, EDGE5, vertices);
                        i3 = _getSharedDualPointIndex(x, y - 1, z, EDGE7, vertices);

                        if (exiting)
                        {
                            triangles.push_back(Triangle(i0, i1, i2));
                            triangles.push_back(Triangle(i2, i3, i0));
                        }
                        else
                        {
                            triangles.push_back(Triangle(i0, i3, i2));
                            triangles.push_back(Triangle(i2, i1, i0));
                        }
                    }
                }
            }
        }

        LOOP_PROGRESS_FRACTION(x, sizeX);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

Mesh* DualMarchingCubes::generateMeshFromVolume(Volume* volume, const bool& serialExecution)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< DualMarchingCubes > workflow = std::make_unique< DualMarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateMesh(!serialExecution);
}

AdvancedMesh* DualMarchingCubes::generateAdvancedMeshFromVolume(Volume* volume,
                                                                const bool& serialExecution)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< DualMarchingCubes > workflow = std::make_unique< DualMarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateAdvancedMesh(!serialExecution);
}


}
