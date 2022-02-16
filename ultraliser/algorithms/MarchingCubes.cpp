/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include "MarchingCubes.h"
#include "MarchingCubes.hh"
#include <data/volumes/voxels/Voxels.h>

namespace Ultraliser
{

MarchingCubes::MarchingCubes(Volume* volume,
                             const uint8_t _isoValue)
    : _volume(volume)
    , _isoValue(_isoValue)
{
    /// EMPTY CONSTRUCTOR
}

double interpolateIsoValue(double _isoValue, double f1, double f2, double t1, double t2)
{
    if(f2 == f1)
        return 0.5 * (t2 + t1);

    // Return interpolated
    return (t2 - t1) * (_isoValue - f1) / (f2 - f1) + t1;
}

uint64_t addSharedVertex(double x1, double y1, double z1,
                       double c2,
                       int axis, double f1, double f2,
                       double _isoValue,
                       Vertices& vertices)
{
    size_t vertexIndex = vertices.size();

    if(axis == 0)
    {
        double x = interpolateIsoValue(_isoValue, f1, f2, x1, c2);
        vertices.push_back(Vector3f(x, y1, z1));
        return vertexIndex;
    }

    if(axis == 1)
    {
        double y = interpolateIsoValue(_isoValue, f1, f2, y1, c2);
        vertices.push_back(Vector3f(x1, y, z1));
        return vertexIndex;
    }

    if(axis == 2)
    {
        double z = interpolateIsoValue(_isoValue, f1, f2, z1, c2);
        vertices.push_back(Vector3f(x1, y1, z));
        return vertexIndex;
    }

    // Just for complication
    return 0;
}

Mesh* MarchingCubes::generateMesh(const bool &parallel)
{
    LOG_TITLE("Mesh Reconstruction with Marching Cubes");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    if (parallel)
        _buildSharedVerticesParallel(vertices, triangles);
    else
        _buildSharedVertices(vertices, triangles);

    // Reconstruct the mesh
    Mesh* mesh = new Mesh(vertices, triangles);

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

AdvancedMesh* MarchingCubes::generateAdvancedMesh(const bool &parallel)
{
    LOG_TITLE("Mesh Reconstruction with Marching Cubes");

    // Build the mesh
    Vertices vertices;
    Triangles triangles;

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Building Mesh");
    if (parallel)
        _buildSharedVerticesParallel(vertices, triangles);
    else
        _buildSharedVertices(vertices, triangles);

    // Reconstruct the mesh
    AdvancedMesh* mesh = new AdvancedMesh(vertices, triangles);

    // Statistics
    _meshExtractionTime = GET_TIME_SECONDS;
    LOG_STATUS_IMPORTANT("Mesh Reconstruction with Marching Cubes Stats.");
    LOG_STATS(_meshExtractionTime);

    return mesh;
}

void MarchingCubes::_buildSharedVerticesParallel(Vertices& vertices, Triangles &triangles)
{
    // Start the timer
    TIMER_SET;

    // Polygons
    std::vector< uint64_t > polygons;

    // Adding a little bit of extra voxels
    const int64_t extraVoxels = 2;
    const int64_t minValue = -1 * extraVoxels;
    const int64_t maxValue = extraVoxels;

    const int64_t maxX = _volume->getWidth() + maxValue;
    const int64_t maxY = _volume->getHeight() + maxValue;
    const int64_t maxZ = _volume->getDepth() + maxValue;

    const uint64_t sizeX = maxX + extraVoxels;
    const uint64_t sizeY = maxY + extraVoxels;
    const uint64_t sizeZ = maxZ + extraVoxels;

    // Total size
    const uint64_t size = (sizeX * sizeY * sizeZ * 3);

    // Mesh shared indices list
    uint64_t* sharedIndices;
    try {
        sharedIndices = new uint64_t[size];

    }  catch (...) {
        LOG_ERROR("Cannot allocate required memory for the Marching Cubes implementation!");
    }

    // A list of lists of DMCVoxel's
    // This list will have reducedX entries to get filled in parallel
    MCVoxelsList volumeMCVoxels;
    volumeMCVoxels.resize(sizeX);

    // For computations and indexing ...
    const uint64_t z3 = sizeZ * 3;
    const uint64_t yz3 = sizeY * z3;

    LOOP_STARTS("Searching Filled Voxels");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = minValue; i < maxX; ++i)
    {
        // Get a reference to the slice
        MCVoxels& sliceMCVoxels = volumeMCVoxels[static_cast< uint64_t >(i + extraVoxels)];

        for (int64_t j = minValue; j < maxY; ++j)
        {
            for (int64_t k = minValue; k < maxZ; ++k)
            {
                double v[8];
                v[0] = _volume->getValue(i, j, k);
                v[1] = _volume->getValue(i + 1, j, k);
                v[2] = _volume->getValue(i + 1, j + 1, k);
                v[3] = _volume->getValue(i, j + 1, k);
                v[4] = _volume->getValue(i, j, k + 1);
                v[5] = _volume->getValue(i + 1, j, k + 1);
                v[6] = _volume->getValue(i + 1, j + 1, k + 1);
                v[7] = _volume->getValue(i, j + 1, k + 1);

                // Get the cube index
                uint64_t cubeIndex = 0;
                for (uint8_t l = 0; l < 8; ++l)
                    if (v[l] <= _isoValue)
                        cubeIndex |= 1 << l;

                if (cubeIndex == 0 || cubeIndex == 255)
                    continue;

                sliceMCVoxels.push_back(new MCVoxel(i, j, k, cubeIndex, v));
            }
        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, sizeX);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Delta between the voxels
    const float voxelSize = _volume->getVoxelSize();

    // Building the shared vertices
    TIMER_RESET;
    LOOP_STARTS("Building Shared Vertices");
    for (uint64_t ii = 0; ii < volumeMCVoxels.size(); ii++)
    {
        LOOP_PROGRESS(ii, volumeMCVoxels.size());
        for (uint64_t jj = 0; jj < volumeMCVoxels[ii].size(); jj++)
        {
            // Reference to the voxel
            MCVoxel* mcVoxel = volumeMCVoxels[ii][jj];

            // Voxel config
            double* v = mcVoxel->vConfig;

            const int64_t i = mcVoxel->x;
            const int64_t j = mcVoxel->y;
            const int64_t k = mcVoxel->z;

            // Geth the voxel bounding box
            Vector3f pMin, pMax;
            _volume->getVoxelBoundingBox(i, j, k, pMin, pMax);

            const float x = pMax.x() - 0.5 * voxelSize;       // dx * i;
            const float x_dx = pMax.x() + 0.5 * voxelSize;    // dx * (i + 1);

            const float y = pMax.y() - 0.5 * voxelSize;       // dy * j;
            const float y_dy = pMax.y() + 0.5 * voxelSize;    // dy * (j + 1);

            const float z = pMax.z() - 0.5 * voxelSize;       // dz * k;
            const float z_dz = pMax.z() + 0.5 * voxelSize;    // dz * (k + 1);

            // Get the edges from the edge table
            int edges = MC_EDGE_TABLE[mcVoxel->cubeIndex];

            // An array of the unique indices per case
            std::array< uint64_t, 12 > uniqueIndices;

            // Generate vertices and avoid DUPLICATE vertices!
            if(edges & 0x040)
            {
                uniqueIndices[6] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 0;
                sharedIndices[idx] = uniqueIndices[6];
                addSharedVertex(x_dx, y_dy, z_dz, x, 0, v[6], v[7], _isoValue, vertices);
            }

            if(edges & 0x020)
            {
                uniqueIndices[5] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 1;
                sharedIndices[idx] = uniqueIndices[5];
                addSharedVertex(x_dx, y, z_dz, y_dy, 1, v[5], v[6], _isoValue, vertices);
            }

            if(edges & 0x400)
            {
                uniqueIndices[10] = vertices.size();
                int64_t idx = i * yz3 + j * z3 + k * 3 + 2;
                sharedIndices[idx] = uniqueIndices[10];
                addSharedVertex(x_dx, y_dy, z, z_dz, 2, v[2], v[6], _isoValue, vertices);
            }

            if(edges & 0x001)
            {
                if(j == 0 || k == 0)
                {
                    uniqueIndices[0] = vertices.size();
                    addSharedVertex(x, y, z, x_dx, 0, v[0], v[1], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = i * yz3 + (j - 1) * z3 + (k - 1) * 3 + 0;
                    uniqueIndices[0] = sharedIndices[idx];
                }
            }

            if(edges & 0x002)
            {
                if(k == 0)
                {
                    uniqueIndices[1] = vertices.size();
                    addSharedVertex(x_dx, y, z, y_dy, 1, v[1], v[2], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 1;
                    uniqueIndices[1] = sharedIndices[idx];
                }
            }

            if(edges & 0x004)
            {
                if(k == 0)
                {
                    uniqueIndices[2] = vertices.size();
                    addSharedVertex(x_dx, y_dy, z, x, 0, v[2], v[3], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 0;
                    uniqueIndices[2] = sharedIndices[idx];
                }
            }

            if(edges & 0x008)
            {
                if(i == 0 || k == 0)
                {
                    uniqueIndices[3] = vertices.size();
                    addSharedVertex(x, y_dy, z, y, 1, v[3], v[0], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + (k - 1) * 3 + 1;
                    uniqueIndices[3] = sharedIndices[idx];
                }
            }

            if(edges & 0x010)
            {
                if(j == 0)
                {
                    uniqueIndices[4] = vertices.size();
                    addSharedVertex(x, y, z_dz, x_dx, 0, v[4], v[5], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 0;
                    uniqueIndices[4] = sharedIndices[idx];
                }
            }

            if(edges & 0x080)
            {
                if(i == 0)
                {
                    uniqueIndices[7] = vertices.size();
                    addSharedVertex(x, y_dy, z_dz, y, 1, v[7], v[4], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 1;
                    uniqueIndices[7] = sharedIndices[idx];
                }
            }

            if(edges & 0x100)
            {
                if(i == 0 || j == 0)
                {
                    uniqueIndices[8] = vertices.size();
                    addSharedVertex(x, y, z, z_dz, 2, v[0], v[4], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + (j - 1) * z3 + k * 3 + 2;
                    uniqueIndices[8] = sharedIndices[idx];
                }
            }

            if(edges & 0x200)
            {
                if(j == 0)
                {
                    uniqueIndices[9] = vertices.size();
                    addSharedVertex(x_dx, y, z, z_dz, 2, v[1], v[5], _isoValue,vertices);
                }
                else
                {
                    int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 2;
                    uniqueIndices[9] = sharedIndices[idx];
                }
            }

            if(edges & 0x800)
            {
                if(i == 0)
                {
                    uniqueIndices[11] = vertices.size();
                    addSharedVertex(x, y_dy, z, z_dz, 2, v[3], v[7], _isoValue, vertices);
                }
                else
                {
                    int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 2;
                    uniqueIndices[11] = sharedIndices[idx];
                }
            }

            int64_t vertexIndex;
            int32_t* triangleTablePtr = MC_TRIANGLE_TABLE[mcVoxel->cubeIndex];
            for(uint64_t tIndex = 0;
                vertexIndex = triangleTablePtr[tIndex], vertexIndex != -1; ++tIndex)
            {
                polygons.push_back(uniqueIndices[vertexIndex]);
            }
        }
    }

    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Construct the triangles from the indices
    for (size_t i = 0; i < polygons.size() / 3; i++)
    {
        Triangle t;
        t[0] = polygons[i * 3 + 0];
        t[1] = polygons[i * 3 + 1];
        t[2] = polygons[i * 3 + 2];
        triangles.push_back(t);
    }

    // Clear
    delete [] sharedIndices;
    polygons.clear();
    polygons.shrink_to_fit();
}

void MarchingCubes::_buildSharedVertices(Vertices& vertices, Triangles &triangles)
{
    // Start the timer
    TIMER_SET;

    // Polygons
    std::vector< uint64_t > polygons;

    // Adding a little bit of extra voxels
    const int64_t extraVoxels = 5;
    const int64_t minValue = -1 * extraVoxels;
    const int64_t maxValue = extraVoxels;

    const int64_t maxX = _volume->getWidth() + maxValue;
    const int64_t maxY = _volume->getHeight() + maxValue;
    const int64_t maxZ = _volume->getDepth() + maxValue;

    const uint64_t sizeX = maxX + extraVoxels;
    const uint64_t sizeY = maxY + extraVoxels;
    const uint64_t sizeZ = maxZ + extraVoxels;

    // Total size
    uint64_t size = ((sizeX) * (sizeY) * (sizeZ) * 3);

    // Mesh shared indices list
    uint64_t* sharedIndices = new size_t[size];

    // For computations and indexing ...
    const uint64_t z3 = sizeZ * 3;
    const uint64_t yz3 = sizeY * z3;

    // Delta between the voxels
    const float voxelSize = _volume->getVoxelSize();

    LOOP_STARTS("Building Shared Vertices");
    for (int64_t i = minValue; i < maxX; ++i)
    {
        for (int64_t j = minValue; j < maxY; ++j)
        {
            for (int64_t k = minValue; k < maxZ; ++k)
            {

                // Geth the voxel bounding box
                Vector3f pMin, pMax;
                _volume->getVoxelBoundingBox(i, j, k, pMin, pMax);

                const float x = pMax.x() - 0.5 * voxelSize;       // dx * i;
                const float x_dx = pMax.x() + 0.5 * voxelSize;    // dx * (i + 1);

                const float y = pMax.y() - 0.5 * voxelSize;       // dy * j;
                const float y_dy = pMax.y() + 0.5 * voxelSize;    // dy * (j + 1);

                const float z = pMax.z() - 0.5 * voxelSize;       // dz * k;
                const float z_dz = pMax.z() + 0.5 * voxelSize;    // dz * (k + 1);

                double v[8];
                v[0] = _volume->getValue(i, j, k);
                v[1] = _volume->getValue(i + 1, j, k);
                v[2] = _volume->getValue(i + 1, j + 1, k);
                v[3] = _volume->getValue(i, j + 1, k);
                v[4] = _volume->getValue(i, j, k + 1);
                v[5] = _volume->getValue(i + 1, j, k + 1);
                v[6] = _volume->getValue(i + 1, j + 1, k + 1);
                v[7] = _volume->getValue(i, j + 1, k + 1);

                // Get the cube index
                uint64_t cubeIndex = 0;
                for (uint8_t l = 0; l < 8; ++l)
                    if (v[l] <= _isoValue)
                        cubeIndex |= 1 << l;

                // Get the edges from the edge table
                int edges = MC_EDGE_TABLE[cubeIndex];

                // An array of the unique indices per case
                std::array< uint64_t, 12 > uniqueIndices;

                // Generate vertices and avoid DUPLICATE vertices!
                if(edges & 0x040)
                {
                    uniqueIndices[6] = vertices.size();
                    int64_t idx = i * yz3 + j * z3 + k * 3 + 0;
                    sharedIndices[idx] = uniqueIndices[6];
                    addSharedVertex(x_dx, y_dy, z_dz, x, 0, v[6], v[7], _isoValue, vertices);
                }

                if(edges & 0x020)
                {
                    uniqueIndices[5] = vertices.size();
                    int64_t idx = i * yz3 + j * z3 + k * 3 + 1;
                    sharedIndices[idx] = uniqueIndices[5];
                    addSharedVertex(x_dx, y, z_dz, y_dy, 1, v[5], v[6], _isoValue, vertices);
                }

                if(edges & 0x400)
                {
                    uniqueIndices[10] = vertices.size();
                    int64_t idx = i * yz3 + j * z3 + k * 3 + 2;
                    sharedIndices[idx] = uniqueIndices[10];
                    addSharedVertex(x_dx, y_dy, z, z_dz, 2, v[2], v[6], _isoValue, vertices);
                }

                if(edges & 0x001)
                {
                    if(j == 0 || k == 0)
                    {
                        uniqueIndices[0] = vertices.size();
                        addSharedVertex(x, y, z, x_dx, 0, v[0], v[1], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = i * yz3 + (j - 1) * z3 + (k - 1) * 3 + 0;
                        uniqueIndices[0] = sharedIndices[idx];
                    }
                }

                if(edges & 0x002)
                {
                    if(k == 0)
                    {
                        uniqueIndices[1] = vertices.size();
                        addSharedVertex(x_dx, y, z, y_dy, 1, v[1], v[2], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 1;
                        uniqueIndices[1] = sharedIndices[idx];
                    }
                }

                if(edges & 0x004)
                {
                    if(k == 0)
                    {
                        uniqueIndices[2] = vertices.size();
                        addSharedVertex(x_dx, y_dy, z, x, 0, v[2], v[3], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = i * yz3 + j * z3 + (k - 1) * 3 + 0;
                        uniqueIndices[2] = sharedIndices[idx];
                    }
                }

                if(edges & 0x008)
                {
                    if(i == 0 || k == 0)
                    {
                        uniqueIndices[3] = vertices.size();
                        addSharedVertex(x, y_dy, z, y, 1, v[3], v[0], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = (i - 1) * yz3 + j * z3 + (k - 1) * 3 + 1;
                        uniqueIndices[3] = sharedIndices[idx];
                    }
                }

                if(edges & 0x010)
                {
                    if(j == 0)
                    {
                        uniqueIndices[4] = vertices.size();
                        addSharedVertex(x, y, z_dz, x_dx, 0, v[4], v[5], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 0;
                        uniqueIndices[4] = sharedIndices[idx];
                    }
                }

                if(edges & 0x080)
                {
                    if(i == 0)
                    {
                        uniqueIndices[7] = vertices.size();
                        addSharedVertex(x, y_dy, z_dz, y, 1, v[7], v[4], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 1;
                        uniqueIndices[7] = sharedIndices[idx];
                    }
                }

                if(edges & 0x100)
                {
                    if(i == 0 || j == 0)
                    {
                        uniqueIndices[8] = vertices.size();
                        addSharedVertex(x, y, z, z_dz, 2, v[0], v[4], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = (i - 1) * yz3 + (j - 1) * z3 + k * 3 + 2;
                        uniqueIndices[8] = sharedIndices[idx];
                    }
                }

                if(edges & 0x200)
                {
                    if(j == 0)
                    {
                        uniqueIndices[9] = vertices.size();
                        addSharedVertex(x_dx, y, z, z_dz, 2, v[1], v[5], _isoValue,vertices);
                    }
                    else
                    {
                        int64_t idx = i * yz3 + (j - 1) * z3 + k * 3 + 2;
                        uniqueIndices[9] = sharedIndices[idx];
                    }
                }

                if(edges & 0x800)
                {
                    if(i == 0)
                    {
                        uniqueIndices[11] = vertices.size();
                        addSharedVertex(x, y_dy, z, z_dz, 2, v[3], v[7], _isoValue, vertices);
                    }
                    else
                    {
                        int64_t idx = (i - 1) * yz3 + j * z3 + k * 3 + 2;
                        uniqueIndices[11] = sharedIndices[idx];
                    }
                }

                int64_t vertexIndex;
                int32_t* triangleTablePtr = MC_TRIANGLE_TABLE[cubeIndex];
                for(uint64_t tIndex = 0;
                    vertexIndex = triangleTablePtr[tIndex], vertexIndex != -1; ++tIndex)
                {
                    polygons.push_back(uniqueIndices[vertexIndex]);
                }
            }
        }

        LOOP_PROGRESS_FRACTION(i, sizeX);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Construct the triangles from the indices
    for (size_t i = 0; i < polygons.size() / 3; i++)
    {
        Triangle t;
        t[0] = polygons[i * 3 + 0];
        t[1] = polygons[i * 3 + 1];
        t[2] = polygons[i * 3 + 2];
        triangles.push_back(t);
    }

    // Clear
    delete [] sharedIndices;
    polygons.clear();
    polygons.shrink_to_fit();
}

Mesh* MarchingCubes::generateMeshFromVolume(Volume* volume, const bool &serialExecution)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< MarchingCubes > workflow = std::make_unique< MarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateMesh(!serialExecution);
}

AdvancedMesh* MarchingCubes::generateAdvancedMeshFromVolume(Volume* volume,
                                                            const bool &serialExecution)
{
    // Reconstruct a watertight mesh from the volume with DMC
    std::unique_ptr< MarchingCubes > workflow = std::make_unique< MarchingCubes >(volume);

    // Generate the DMC mesh
    return workflow->generateAdvancedMesh(!serialExecution);
}

}
