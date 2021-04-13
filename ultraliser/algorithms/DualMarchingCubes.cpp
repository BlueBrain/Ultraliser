/***************************************************************************************************
 * Copyright (c) 2018 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

DualMarchingCubes::DualMarchingCubes(Volume* volume,
                                     const uint8_t isoValue,
                                     const bool& generateManifold)
    : _volume(volume)
    , _isoValue(isoValue)
    , _generateManifold(generateManifold)
{
    /// EMPTY CONSTRUCTOR
}

Mesh* DualMarchingCubes::generateMesh()
{
    LOG_TITLE("Mesh Reconstruction with DMC");

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;
    _buildSharedVerticesQuadsParallel(vertices, triangles);

    Mesh* mesh = new Mesh(vertices, triangles);

    // Statistics
    _dmcGenerationTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with DMC Stats.");
    LOG_STATS(_dmcGenerationTime);

    return mesh;
}

std::unique_ptr< Mesh > DualMarchingCubes::generateMeshX()
{
    LOG_TITLE("Mesh Reconstruction with DMC");

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;
    _buildSharedVerticesQuadsParallel(vertices, triangles);

    std::unique_ptr< Mesh > mesh = std::make_unique< Mesh >(vertices, triangles);

    // Statistics
    _dmcGenerationTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with DMC Stats.");
    LOG_STATS(_dmcGenerationTime);

    return mesh;
}


AdvancedMesh* DualMarchingCubes::generateManifoldMesh()
{
    // Build the mesh
    Vertices vertices;
    Triangles triangles;
    _buildSharedVerticesQuadsParallel(vertices, triangles);

    // Return a new mesh
    return new AdvancedMesh(vertices, triangles);
}

int DualMarchingCubes::_getCellCode(const int64_t &x,
                                    const int64_t &y,
                                    const int64_t &z) const
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

int64_t DualMarchingCubes::_index(const int64_t &x,
                                  const int64_t &y,
                                  const int64_t &z) const
{
    return x + _volume->getWidth() * (y + _volume->getHeight() * z);
}

int DualMarchingCubes::_getDualPointCode(const int64_t x,
                                         const int64_t y,
                                         const int64_t z,
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

void DualMarchingCubes::_calculateDualPoint(const int64_t &x,
                                            const int64_t &y,
                                            const int64_t &z,
                                            const int32_t& pointCode,
                                            Vector3f & v) const
{
    // Initialize the point with lower voxel coordinates
    v.x() = x;
    v.y() = y;
    v.z() = z;

    // Compute the dual point as the mean of the face vertices belonging to the
    // original marching cubes face
    Vertex p;
    p.x() = 0; p.y() = 0; p.z() = 0;

    int points = 0;

    // Sum edge intersection vertices using the point code
    if (I2UI32(pointCode) & EDGE0)
    {
        p.x() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x + 1, y, z) -
                 _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE1)
    {
        p.x() += 1.0f;
        p.z() += (_isoValue - _volume->getValue(x + 1, y, z)) /
                (_volume->getValue(x + 1, y, z + 1) -
                 _volume->getValue(x + 1, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE2)
    {
        p.x() += (_isoValue - _volume->getValue(x, y, z + 1)) /
                (_volume->getValue(x + 1, y, z + 1) -
                 _volume->getValue(x, y, z + 1));
        p.z() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE3)
    {
        p.z() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x, y, z + 1) -
                 _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE4)
    {
        p.x() += (_isoValue - _volume->getValue(x, y + 1, z)) /
                (_volume->getValue(x + 1, y + 1, z) -
                 _volume->getValue(x, y + 1, z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE5)
    {
        p.x() += 1.0f;
        p.z() += (_isoValue - _volume->getValue(x + 1, y + 1, z)) /
                (_volume->getValue(x + 1, y + 1, z + 1) -
                 _volume->getValue(x + 1, y + 1, z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE6)
    {
        p.x() += (_isoValue - _volume->getValue(x, y + 1, z + 1)) /
                (_volume->getValue(x + 1, y + 1, z + 1) -
                 _volume->getValue(x, y + 1, z + 1));
        p.z() += 1.0f;
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE7)
    {
        p.z() += (_isoValue - _volume->getValue(x, y + 1 , z)) /
                (_volume->getValue(x, y + 1, z + 1) -
                 _volume->getValue(x, y + 1 , z));
        p.y() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE8)
    {
        p.y() += (_isoValue - _volume->getValue(x, y, z)) /
                (_volume->getValue(x, y + 1, z) -
                 _volume->getValue(x, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE9)
    {
        p.x() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x + 1, y, z)) /
                (_volume->getValue(x + 1, y + 1, z) -
                 _volume->getValue(x + 1, y, z));
        points++;
    }

    if (I2UI32(pointCode) & EDGE10)
    {
        p.x() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x + 1, y, z + 1)) /
                (_volume->getValue(x + 1, y + 1, z + 1) -
                 _volume->getValue(x + 1, y, z + 1));
        p.z() += 1.0f;
        points++;
    }

    if (I2UI32(pointCode) & EDGE11)
    {
        p.z() += 1.0f;
        p.y() += (_isoValue - _volume->getValue(x, y, z + 1)) /
                (_volume->getValue(x, y + 1, z + 1) -
                 _volume->getValue(x, y, z + 1));
        points++;
    }

    // Divide by number of accumulated points
    float invPoints = 1.0f / points;
    p.x() *= invPoints; p.y() *= invPoints; p.z() *= invPoints;

    // Offset point by voxel coordinates
    v.x() += p.x();
    v.y() += p.y();
    v.z() += p.z();
}

int64_t DualMarchingCubes::_getSharedDualPointIndex(
        const int64_t &x, const int64_t &y, const int64_t &z,
        const DMC_EDGE_CODE& edge, std::vector<Vector3f> & vertices)
{
    // Create a key for the dual point from its linearized cell ID and
    // point code
    DualPointKey key;
    key.linearizedCellID = _index(x, y, z);
    key.pointCode = _getDualPointCode(x, y, z, edge);

    // have we already computed the dual point?
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
    return (linearizedCellID == other.linearizedCellID &&
            pointCode == other.pointCode);
}

void DualMarchingCubes::_buildSharedVerticesQuadsParallel(Vertices& vertices,
                                                          Triangles &triangles)
{
    // Start timer
    TIMER_SET;

    // Compute the reduced dimensions of the volume that would fit
    const int64_t reducedX = I2I64(_volume->getWidth()) - 2;
    const int64_t reducedY = I2I64(_volume->getHeight()) - 2;
    const int64_t reducedZ = I2I64(_volume->getDepth()) - 2;

    // A list of lists of DMCVoxel's
    // This list will have reducedX entries to get filled in parallel
    DMCVoxelsList volumeDMCVoxels;
    volumeDMCVoxels.resize(static_cast< uint64_t >(reducedX));

    // Shared quad points
    int64_t i0, i1, i2, i3;

    // Clear the map before processing any items
    pointToIndex.clear();

    // Searching for non-zero voxels in parallel
    LOOP_STARTS("Searching Filled Voxels");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
     for (int64_t x = 0; x < reducedX; ++x)
     {
         // Get a reference to the slice
         DMCVoxels& sliceDMCVoxels =
                 volumeDMCVoxels[static_cast< uint64_t >(x)];

#ifdef ULTRALISER_USE_OPENMP
         if (omp_get_thread_num() == 0)
#endif
         {
             LOOP_PROGRESS(progress, reducedX * reducedY * reducedZ);
         }

         for (int64_t y = 0; y < reducedY; ++y)
         {
             for (int64_t z = 0; z < reducedZ; ++z)
             {
#ifdef ULTRALISER_USE_OPENMP
                #pragma omp atomic
#endif
                ++progress;

                // X-aligned edge
                if (z > 0 && y > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x + 1, y, z) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x + 1, y, z) < _isoValue;

                    // Create a DMCVoxel if there is one
                    if (entering || exiting)
                    {
                        DMCVoxel dmcVoxel;
                        dmcVoxel.x = x;
                        dmcVoxel.y = y;
                        dmcVoxel.z = z;
                        dmcVoxel.entering = entering;
                        dmcVoxel.exiting = exiting;
                        dmcVoxel.side = DMC_EDGE_SIDE::X;

                        sliceDMCVoxels.push_back(dmcVoxel);
                    }
                }

                // Y-aligned edge
                if (z > 0 && x > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x, y + 1, z) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x, y + 1, z) < _isoValue;

                    if (entering || exiting)
                    {
                        // That's a valid voxel
                        DMCVoxel dmcVoxel;
                        dmcVoxel.x = x;
                        dmcVoxel.y = y;
                        dmcVoxel.z = z;
                        dmcVoxel.entering = entering;
                        dmcVoxel.exiting = exiting;
                        dmcVoxel.side = DMC_EDGE_SIDE::Y;

                        sliceDMCVoxels.push_back(dmcVoxel);
                    }
                }

                // Z-aligned edge
                if (x > 0 && y > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x, y, z + 1) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x, y, z + 1) < _isoValue;

                    if (entering || exiting)
                    {
                        // That's a valid voxel
                        DMCVoxel dmcVoxel;
                        dmcVoxel.x = x;
                        dmcVoxel.y = y;
                        dmcVoxel.z = z;
                        dmcVoxel.entering = entering;
                        dmcVoxel.exiting = exiting;
                        dmcVoxel.side = DMC_EDGE_SIDE::Z;

                        sliceDMCVoxels.push_back(dmcVoxel);
                    }
                }
            }
        }
    }
    LOOP_DONE;

    // Statistics
    _searchingZeroVoxelsTime = GET_TIME_SECONDS;
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
            DMCVoxel& dmcVoxel = volumeDMCVoxels[i][j];

            // X-aligned edge
            if (dmcVoxel.side == DMC_EDGE_SIDE::X)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y,
                                               dmcVoxel.z,
                                               EDGE0, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y,
                                               dmcVoxel.z - 1,
                                               EDGE2, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y - 1,
                                               dmcVoxel.z - 1,
                                               EDGE6, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y - 1,
                                               dmcVoxel.z,
                                               EDGE4, vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel.entering)
                {
                    // quads.push_back(Vec4i_64(i0, i1, i2, i3));
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]

                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // quads.push_back(Vec4i_64(i0, i3, i2, i1));
                    // Quad [i0, i3, i2, i1]
                    // Triangles [i0, i3, i2] [i2, i1, i0]

                    triangles.push_back(Triangle(i0, i3, i2));
                    triangles.push_back(Triangle(i2, i1, i0));
                }
            }

            // Y-aligned edge
            else if (dmcVoxel.side == DMC_EDGE_SIDE::Y)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel.x,
                                              dmcVoxel.y,
                                              dmcVoxel.z,
                                              EDGE8, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel.x,
                                              dmcVoxel.y,
                                              dmcVoxel.z - 1,
                                              EDGE11, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel.x - 1,
                                              dmcVoxel.y,
                                              dmcVoxel.z - 1,
                                              EDGE10, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel.x - 1,
                                              dmcVoxel.y,
                                              dmcVoxel.z,
                                              EDGE9, vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel.exiting)
                {
                    // quads.push_back(Vec4i_64(i0, i1, i2, i3));
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]

                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // quads.push_back(Vec4i_64(i0, i3, i2, i1));
                    // Quad [i0, i3, i2, i1]
                    // Triangles [i0, i3, i2] [i2, i1, i0]

                    triangles.push_back(Triangle(i0, i3, i2));
                    triangles.push_back(Triangle(i2, i1, i0));
                }
            }

            // Z-aligned edge
            else if (dmcVoxel.side == DMC_EDGE_SIDE::Z)
            {
                // Generate quad
                i0 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y,
                                               dmcVoxel.z,
                                               EDGE3, vertices);
                i1 = _getSharedDualPointIndex(dmcVoxel.x - 1,
                                               dmcVoxel.y,
                                               dmcVoxel.z,
                                               EDGE1, vertices);
                i2 = _getSharedDualPointIndex(dmcVoxel.x - 1,
                                               dmcVoxel.y - 1,
                                               dmcVoxel.z,
                                               EDGE5, vertices);
                i3 = _getSharedDualPointIndex(dmcVoxel.x,
                                               dmcVoxel.y - 1,
                                               dmcVoxel.z,
                                               EDGE7,
                                               vertices);

                // Append this quad to the list based on the voxel
                if (dmcVoxel.exiting)
                {
                    // quads.push_back(Vec4i_64(i0, i1, i2, i3));
                    // Quad [i0, i1, i2, i3]
                    // Triangles [i0, i1, i2] [i2, i3, i0]

                    triangles.push_back(Triangle(i0, i1, i2));
                    triangles.push_back(Triangle(i2, i3, i0));
                }
                else
                {
                    // quads.push_back(Vec4i_64(i0, i3, i2, i1));
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
        }
    }
    LOOP_DONE;

    // Statistics
    _buildingSharedVerticesTime = GET_TIME_SECONDS;
    LOG_STATS(GET_TIME_SECONDS);

    _buildSharedVerticesQuadsParallelTime =
            _searchingZeroVoxelsTime + _buildingSharedVerticesTime;
}

void DualMarchingCubes::_buildSharedVerticesQuads(Vertices &vertices,
                                                  Triangles &triangles)
{
    // Start timer
    TIMER_SET;

    int64_t const reducedX = _volume->getWidth() - 2;
    int64_t const reducedY = _volume->getHeight() - 2;
    int64_t const reducedZ = _volume->getDepth() - 2;

    // Quad points
    int64_t i0, i1, i2, i3;

    pointToIndex.clear();

    // Iterate voxels
    LOOP_STARTS("Building Shared Vertices");
    uint64_t index = 0;
    for (int64_t x = 0; x < reducedX; ++x)
    {
        for (int64_t y = 0; y < reducedY; ++y)
        {
            for (int64_t z = 0; z < reducedZ; ++z)
            {
                index++;

                LOOP_PROGRESS(index, reducedZ * reducedY * reducedX);

                // Construct quads for X edge
                if (z > 0 && y > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x + 1, y, z) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x + 1, y, z) < _isoValue;

                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z,
                                                      EDGE0, vertices);
                        i1 = _getSharedDualPointIndex(x, y, z - 1,
                                                      EDGE2, vertices);
                        i2 = _getSharedDualPointIndex(x, y - 1, z - 1,
                                                      EDGE6, vertices);
                        i3 = _getSharedDualPointIndex(x, y - 1, z,
                                                      EDGE4, vertices);

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
                if (z > 0 && x > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x, y + 1, z) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x, y + 1, z) < _isoValue;

                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z,
                                                      EDGE8, vertices);
                        i1 = _getSharedDualPointIndex(x, y, z - 1,
                                                      EDGE11, vertices);
                        i2 = _getSharedDualPointIndex(x - 1, y, z - 1,
                                                      EDGE10, vertices);
                        i3 = _getSharedDualPointIndex(x - 1, y, z,
                                                      EDGE9, vertices);

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
                if (x > 0 && y > 0)
                {
                    bool const entering =
                            _volume->getValue(x, y, z) < _isoValue &&
                            _volume->getValue(x, y, z + 1) >= _isoValue;
                    bool const exiting  =
                            _volume->getValue(x, y, z) >= _isoValue &&
                            _volume->getValue(x, y, z + 1) < _isoValue;
                    if (entering || exiting)
                    {
                        // Generate quad
                        i0 = _getSharedDualPointIndex(x, y, z,
                                                      EDGE3, vertices);
                        i1 = _getSharedDualPointIndex(x - 1, y, z,
                                                      EDGE1, vertices);
                        i2 = _getSharedDualPointIndex(x - 1, y - 1, z,
                                                      EDGE5, vertices);
                        i3 = _getSharedDualPointIndex(x, y - 1, z,
                                                      EDGE7, vertices);

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
    }
    LOOP_DONE;

    // Statistics
    _buildSharedVerticesQuadsSerialTime = GET_TIME_SECONDS;
    LOG_STATS(GET_TIME_SECONDS);
}

}
