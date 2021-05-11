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

#include "Volume.h"
#include "../grids/VolumeGrid.h"
#include "../voxels/DMCVoxel.h"
#include <common/Common.h>
#include <geometry/Intersection.h>
#include <geometry/Utilities.h>
#include <utilities/Utilities.h>
#include <math/Functions.h>
#include <data/images/TIFFImage.h>
#include "../grids/BitVolumeGrid.h"
#include "../grids/ByteVolumeGrid.h"
#include "../grids/VoxelGrid.h"
#include <data/volumes/grids/Grids.h>
#include <data/meshes/simple/VolumeMesh.h>
#include <data/meshes/simple/MeshOperations.h>

namespace Ultraliser
{

Volume::Volume(const Vector3f& pMin,
               const Vector3f& pMax,
               const uint64_t &baseResolution,
               const float &voxelPadding,
               const VolumeGrid::TYPE& gridType)
    : _gridType(gridType)
    , _pMin(pMin)
    , _pMax(pMax)
    , _voxelPadding(voxelPadding)
    , _baseResolution(baseResolution)
{
    // Create the grid
    _createGrid();

    // Profiling
    _surfaceVoxelizationTime = 0.0;
    _solidVoxelizationTime = 0.0;
    _addingVolumePassTime = 0.0;
}

Volume::Volume(const int64_t width,
               const int64_t height,
               const int64_t depth,
               const Vector3f pMin,
               const Vector3f pMax,
               const VolumeGrid::TYPE& gridType)
    : _gridType(gridType)
{
    _gridDimensions.v[0] = width;
    _gridDimensions.v[1] = height;
    _gridDimensions.v[2] = depth;

    _pMin = pMin;
    _pMax = pMax;

    // Since we don't have any geometric bounds, use 1.0 for the voxel resolution
    // TODO: The voxel resolution should be computed from the given bounds
    _voxelSize = 1.0;

    // Allocating the grid
    _allocateGrid();
}

Volume::Volume(const std::string &prefix,
               const VolumeGrid::TYPE& gridType)
    : _gridType(gridType)
{
    // Zero zero-padding
    _voxelPadding = 0;

    // Load the volume data
    switch (gridType)
    {
    case VolumeGrid::TYPE::BIT:
        _loadBinaryVolumeData(prefix);
        break;
    case VolumeGrid::TYPE::BYTE:
        _loadByteVolumeData(prefix);
        break;

    default:
        break;
    }

    // TODO: Adjust the pMin and pMax here based on a unit cube
}

void Volume::_allocateGrid()
{
    // Create the grid
    switch (_gridType)
    {
    case VolumeGrid::TYPE::BIT:
    {
        _grid = new BitVolumeGrid(_gridDimensions.v[0], _gridDimensions.v[1], _gridDimensions.v[2]);
        break;
    }
    case VolumeGrid::TYPE::BYTE:
    {
        _grid = new ByteVolumeGrid(_gridDimensions.v[0],
                                   _gridDimensions.v[1],
                                   _gridDimensions.v[2]);
        break;
    }

    case VolumeGrid::TYPE::VOXEL:
    {
        _grid = new VoxelGrid(_gridDimensions.v[0], _gridDimensions.v[1], _gridDimensions.v[2]);
        break;
    }
    }
}

void Volume::_createGrid(void)
{
    // Compute the bounding box of the mesh
    Vector3f bbox = _pMax - _pMin;
    Vector3f delta;

//    if (isZero(_voxelPadding))
//    {
//        // If the resolution is greater than 64, use 2.5% padding
//        if (_baseResolution > 64)
//            delta = bbox * 0.025f;

//        // Otherwise, use 10%
//        else
//            delta = bbox * 0.25f;
//    }

//    // Use the specified voxel padding
//    else
//        delta = bbox * _voxelPadding;

//    // Increment the bounding box
//    _pMax += delta;
//    _pMin -= delta;

    // Compute the bounding box size
    Vector3f boundingBoxSize = (_pMax - _pMin);

    // Find the largest dimension of the mesh model to be able to create a
    // scaled grid.
    _largestDimensionIdx = getLargestDimension(boundingBoxSize);

    // Compute the voxel size
    _voxelSize = (boundingBoxSize[_largestDimensionIdx]) / I2F(_baseResolution);

    int64_t width = D2I64(std::ceil(boundingBoxSize[0] / _voxelSize));
    int64_t height = D2I64(std::ceil(boundingBoxSize[1] / _voxelSize));
    int64_t depth = D2I64(std::ceil(boundingBoxSize[2] / _voxelSize));

    // Compute the volume grid dimensions
    const int64_t extraVoxels = 0;
    _gridDimensions.v[0] = width + extraVoxels;
    _gridDimensions.v[1] = height + extraVoxels;
    _gridDimensions.v[2] = depth + extraVoxels;

    // Allocating the grid
    _allocateGrid();

    // Compute the origin of the mesh.
    // NOTE: This is not the center of the bounding box, it is just the point
    // at the base of the bounding box
    _meshOrigin = _pMin - (0.5f * extraVoxels) * Vector3f(_voxelSize);
}

void Volume::_loadHeaderData(const std::string &prefix)
{
    const std::string filePath = prefix + std::string(".hdr");

    // Open the header file
    std::ifstream hdrFileStream(filePath.c_str());

    // Volume dimensions
    hdrFileStream >> _gridDimensions.v[0];
    hdrFileStream >> _gridDimensions.v[1];
    hdrFileStream >> _gridDimensions.v[2];

    // Close the stream
    hdrFileStream.close();
}

void Volume::_loadByteVolumeData(const std::string &prefix)
{
    // Load the header data of the volume including the dimensions and the size
    _loadHeaderData(prefix);

    // We can safely allocate the grid after loading the data from teh header
    _allocateGrid();

    // Load the data from the file in the grid
    _grid->loadByteVolumeData(prefix);
}

void Volume::_loadBinaryVolumeData(const std::string &prefix)
{
    _loadHeaderData(prefix);

    // Allocate the grid with the loaded dimensions
    _allocateGrid();

    // Load the data from the file in the grid
    _grid->loadBinaryVolumeData(prefix);
}

void Volume::surfaceVoxelization(Mesh* mesh,
                                 const bool& verbose,
                                 const bool parallel)
{
    if (verbose) LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    if (verbose) LOG_STATUS("Creating Volume Shell");
    if (parallel)
        _rasterizeParallel(mesh, _grid);
    else
        _rasterize(mesh , _grid, verbose);

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    if (verbose)
    {
        LOG_STATUS_IMPORTANT("Rasterization Stats.");
        LOG_STATS(_surfaceVoxelizationTime);
    }
}

void Volume::surfaceVoxelizeVasculatureMorphology(
        VasculatureMorphology* vasculatureMorphology, const bool& verbose,
        const bool parallel)
{
    if (verbose) LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    // Get all the sections of the vascular morphology
    Sections sections = vasculatureMorphology->getSections();

    LOG_STATUS("Creating Volume Shell from Sections");
    uint64_t processedSections = 0;
    LOOP_STARTS("Rasterization");
    #pragma omp parallel for
    for (uint64_t i = 0; i < sections.size(); i++)
    {
        #pragma omp atomic
        processedSections++;

        LOOP_PROGRESS(processedSections, sections.size());

        // Construct the paths
        Paths paths = vasculatureMorphology->getConnectedPathsFromParentsToChildren(sections[i]);
        for (uint64_t j = 0; j < paths.size(); ++j)
        {
            auto mesh = new Ultraliser::Mesh(paths[j]);
            _rasterize(mesh , _grid, verbose);
        }
    }
    LOOP_DONE;

    if (verbose) LOG_STATUS("Creating Volume Shell");

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    if (verbose)
    {
        LOG_STATUS_IMPORTANT("Rasterization Stats.");
        LOG_STATS(_surfaceVoxelizationTime);
    }
}

void Volume::surfaceVoxelization(AdvancedMesh *mesh)
{
    LOG_TITLE("Surface Voxelization");

    // Start the timer
    TIMER_SET;

    LOG_STATUS("Creating Volume Shell");
    _rasterize(mesh , _grid);

    _surfaceVoxelizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATS(_surfaceVoxelizationTime);
}

void Volume::surfaceVoxelization(const std::string &inputDirectory,
                                 const std::vector< std::string>& meshFiles)
{
    LOG_TITLE("Surface Voxelization");
    TIMER_SET;

    LOG_STATUS("Creating Volume Shell");
    uint64_t processedMeshCount = 0;
    LOOP_STARTS("Rasterization");
    #pragma omp parallel for schedule( dynamic, 1 )
    for( size_t iMesh = 0; iMesh < meshFiles.size(); iMesh++ )
    {
        // Create and load the mesh from the file
        std::string meshName = meshFiles[ iMesh ];
        std::string meshFile = meshName;
        if( inputDirectory != EMPTY )
            meshFile = inputDirectory + "/" + meshFile;

        if( Ultraliser::File::exists(meshFile))
        {
            #pragma omp atomic
            processedMeshCount++;

            LOOP_PROGRESS(processedMeshCount, meshFiles.size());

            // Neuron mesh
            Mesh* mesh = new Mesh(meshFile);

            // Surface voxelization
            const bool verbose = false;
            surfaceVoxelization(mesh, verbose);

            // Free the mesh
            mesh->~Mesh();
        }
    }
    LOOP_DONE;

    // Statistics
    _surfaceVoxelizationTime = GET_TIME_SECONDS;
    LOG_STATS(_surfaceVoxelizationTime);
    LOG_DETAIL("[%zu/%zu] Meshes were Voxelized with Surface Voxelization",
               processedMeshCount, meshFiles.size());

}

void Volume::_rasterize(Mesh* mesh, VolumeGrid* grid, const bool& verbose)
{
    if (verbose) LOOP_STARTS("Rasterization");
    size_t progress = 0;
    for (size_t tIdx = 0; tIdx < mesh->getNumberTriangles(); tIdx++)
    {
        ++progress;
        if (verbose) LOOP_PROGRESS(progress, mesh->getNumberTriangles());

        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, tIdx, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ix++)
        {
            for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; iy++)
            {
                for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; iz++)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleCubeIntersection(mesh, tIdx, gi))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }
    }
    if (verbose) LOOP_DONE;
}

void Volume::_rasterizeParallel(Mesh* mesh, VolumeGrid* grid)
{
    // Start the timer
    TIMER_SET;

    LOOP_STARTS("Parallel Rasterization");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t tIdx = 0; tIdx < mesh->getNumberTriangles(); tIdx++)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, mesh->getNumberTriangles());

        // Get the pMin and pMax of the triangle within the grid
        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getBoundingBox(mesh, tIdx, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0]; ix <= pMaxTriangle[0]; ix++)
        {
            for (int64_t iy = pMinTriangle[1]; iy <= pMaxTriangle[1]; iy++)
            {
                for (int64_t iz = pMinTriangle[2]; iz <= pMaxTriangle[2]; iz++)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleCubeIntersection(mesh, tIdx, gi))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::_rasterize(AdvancedMesh* mesh, VolumeGrid* grid)
{
    // Get a an array of triangles
    AdvancedTriangle** triangles = (AdvancedTriangle **) mesh->_triangles.toArray();
    int triangleCount = mesh->_triangles.numberElements();

    LOOP_STARTS("Rasterization");
    size_t progress = 0;
    for (int i = 0; i <triangleCount; ++i)
    {
        AdvancedTriangle triangle = *triangles[i];

        ++progress;
        LOOP_PROGRESS(progress, triangleCount);

        int64_t pMinTriangle[3], pMaxTriangle[3];
        _getTriangleBoundingBox(triangle, pMinTriangle, pMaxTriangle);

        for (int64_t ix = pMinTriangle[0];
             ix <= (pMaxTriangle[0]); ix++)
        {
            for (int64_t iy = pMinTriangle[1];
                 iy <= (pMaxTriangle[1]); iy++)
            {
                for (int64_t iz = pMinTriangle[2];
                     iz <= (pMaxTriangle[2]); iz++)
                {
                    GridIndex gi(I2I64(ix), I2I64(iy), I2I64(iz));
                    if (_testTriangleGridIntersection(triangle, gi))
                        grid->fillVoxel(I2I64(ix), I2I64(iy), I2I64(iz));
                }
            }
        }
    }
    LOOP_DONE;
}

void Volume::solidVoxelization(const SOLID_VOXELIZATION_AXIS& axis)
{
    LOG_TITLE("Solid Voxelization");

    // The 2D flood filling is only supported for the solid voxelization
    LOG_STATUS("Filling Volume");
    _floodFill2D(axis);

    LOG_STATUS_IMPORTANT("Solid Voxelization Stats.");
    LOG_STATS(_solidVoxelizationTime);
}

void Volume::_floodFill2D(const SOLID_VOXELIZATION_AXIS &axis)
{
    // Start the timer
    TIMER_SET;

    switch (axis)
    {
    case X: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::X);
        break;

    case Y: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::Y);
        break;

    case Z: _floodFillAlongAxis(_grid, SOLID_VOXELIZATION_AXIS::Z);
        break;

    case XYZ:_floodFillAlongXYZ(_grid);
        break;
    }

    // Save the solid voxelization time
    _solidVoxelizationTime = GET_TIME_SECONDS;
}

void Volume::_floodFillAlongAxis(VolumeGrid* grid, const SOLID_VOXELIZATION_AXIS &axis)
{
    /// Disable buffering
    setbuf(stdout, nullptr);

    // Start the timer
    TIMER_SET;

    // The dimension with which the flood filling will happen
    int64_t dimension;

    // Flood-filling string
    std::string floodFillingString;

    // Flood-filling axis
    AXIS floodFillingAxis;

    switch (axis)
    {
    case SOLID_VOXELIZATION_AXIS::X:
    {
        dimension = getWidth();
        floodFillingAxis = AXIS::X;
        floodFillingString = "2D Slice Flood-filling (X-axis)";
    } break;

    case SOLID_VOXELIZATION_AXIS::Y:
    {
        dimension = getHeight();
        floodFillingAxis = AXIS::Y;
        floodFillingString = "2D Slice Flood-filling (Y-axis)";
    } break;

    case SOLID_VOXELIZATION_AXIS::Z:
    {
        dimension = getDepth();
        floodFillingAxis = AXIS::Z;
        floodFillingString = "2D Slice Flood-filling (Z-axis)";
    } break;

    // XYZ voxelization will be handled
    case SOLID_VOXELIZATION_AXIS::XYZ:
        break;
    }

    LOOP_STARTS(floodFillingString.c_str());
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (int64_t i = 0 ; i < dimension; ++i)
    {
        grid->floodFillSliceAlongAxis(i, floodFillingAxis);

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, dimension);
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::_floodFillAlongXYZ(VolumeGrid *grid)
{
    // Volume grids per axis
    VolumeGrid *xGrid, *yGrid, *zGrid;

    // Create the grid
    switch (_gridType)
    {
    case VolumeGrid::TYPE::BIT:
    {
        xGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));
        yGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));
        zGrid = new BitVolumeGrid(static_cast< BitVolumeGrid* >(grid));

    } break;

    case VolumeGrid::TYPE::BYTE:
    {
        xGrid = new ByteVolumeGrid(static_cast< ByteVolumeGrid* >(grid));
        yGrid = new ByteVolumeGrid(static_cast< ByteVolumeGrid* >(grid));
        zGrid = new ByteVolumeGrid(static_cast< ByteVolumeGrid* >(grid));
    } break;

    case VolumeGrid::TYPE::VOXEL:
    {
        xGrid = new VoxelGrid(static_cast< VoxelGrid* >(grid));
        yGrid = new VoxelGrid(static_cast< VoxelGrid* >(grid));
        zGrid = new VoxelGrid(static_cast< VoxelGrid* >(grid));
    } break;
    }

    // Flood fill along the three axes
    _floodFillAlongAxis(xGrid, SOLID_VOXELIZATION_AXIS::X);
    _floodFillAlongAxis(yGrid, SOLID_VOXELIZATION_AXIS::Y);
    _floodFillAlongAxis(zGrid, SOLID_VOXELIZATION_AXIS::Z);

    // Blend the three grids using AND operation and store the final result in the xGrid
    xGrid->andWithAnotherGrid(yGrid);
    xGrid->andWithAnotherGrid(zGrid);

    // Blend the xGrid with the default _grid
    _grid->orWithAnotherGrid(xGrid);

    // Release the auxiliary grids
    xGrid->~VolumeGrid();
    yGrid->~VolumeGrid();
    zGrid->~VolumeGrid();
}

int Volume::_triangleCubeSign(Mesh *mesh,
                              int tIdx, const GridIndex & gi)
{
    Vector3f boxcenter((0.5f + gi[0]) * _voxelSize + _meshOrigin[0],
                       (0.5f + gi[1]) * _voxelSize + _meshOrigin[1],
                       (0.5f + gi[2]) * _voxelSize + _meshOrigin[2]);

    Vector3f tv[3];
    for (size_t ii = 0; ii < 3; ++ii)
        tv[ii] = mesh->getVertices()[mesh->getTriangles()[tIdx][ii]];

    Vector3f e1 = tv[1] - tv[0];
    Vector3f e2 = tv[2] - tv[0];
    Vector3f n = Vector3f::cross(e1, e2).normalized();
    Vector3f d = boxcenter - tv[0];
    float dotp = Vector3f::dot(n, d);
    if (dotp > 0)
        return 2;

    // Too far away
    if (dotp < -0.9f * _voxelSize)
        return 3;

    n.normalize();
    d = boxcenter - (Vector3f::dot(n, d)) * n;
    Vector3f n0 = Vector3f::cross(tv[1] - d, tv[2] - d);

    // n0.normalize();
    float thresh = 1.0;

    if (Vector3f::dot(n0, n) < -thresh)
        return 1;

    Vector3f n1 = Vector3f::cross(tv[2] - d, tv[0] - d);
    // n1.normalize();
    if (Vector3f::dot(n1, n) < -thresh)
        return 1;

    Vector3f n2 = Vector3f::cross(tv[0] - d, tv[1] - d);
    // n2.normalize();
    if (Vector3f::dot(n2, n)< -thresh)
        return 1;

    return -1;
}

bool Volume::_testTriangleCubeIntersection(Mesh* mesh,
                                           uint64_t triangleIdx,
                                           const GridIndex& voxel)
{
    // Voxel center
    double voxelCenter[3];
    voxelCenter[0] = F2D((0.5f + voxel[0]) * _voxelSize + _meshOrigin[0]);
    voxelCenter[1] = F2D((0.5f + voxel[1]) * _voxelSize + _meshOrigin[1]);
    voxelCenter[2] = F2D((0.5f + voxel[2]) * _voxelSize + _meshOrigin[2]);

    // Voxel half size
    double voxelHalfSize[3];
    voxelHalfSize[0] = F2D(_voxelSize * 0.5f);
    voxelHalfSize[1] = F2D(_voxelSize * 0.5f);
    voxelHalfSize[2] = F2D(_voxelSize * 0.5f);

    // Triangle vertices
    double triangle[3][3];

    // For each vertex in the triangle
    for (size_t i = 0; i < 3; ++i)
    {
        // For each coordinate of the vertex
        for (size_t j = 0; j < 3; ++j)
        {
            // Load all the verticies of the selected triangle in _triangle_
            triangle[i][j] =
                    mesh->getVertices()[mesh->getTriangles()[triangleIdx][i]][j];
        }
    }

    // Test if the triangle and the voxel are intersecting or not
    return checkTriangleBoxIntersection(voxelCenter, voxelHalfSize, triangle);
}

bool Volume::_testTriangleGridIntersection(AdvancedTriangle triangle,
                                           const GridIndex& voxel)
{
    // Voxel center
    double voxelCenter[3];
    voxelCenter[0] = F2D((0.5f + voxel[0]) * _voxelSize + _meshOrigin[0]);
    voxelCenter[1] = F2D((0.5f + voxel[1]) * _voxelSize + _meshOrigin[1]);
    voxelCenter[2] = F2D((0.5f + voxel[2]) * _voxelSize + _meshOrigin[2]);

    // Voxel half size
    double voxelHalfSize[3];
    voxelHalfSize[0] = F2D(_voxelSize * 0.5f);
    voxelHalfSize[1] = F2D(_voxelSize * 0.5f);
    voxelHalfSize[2] = F2D(_voxelSize * 0.5f);

    // Construct an array to represent the data
    double triangleArray[3][3];
    triangleArray[0][0] = triangle.v1()->x;
    triangleArray[0][1] = triangle.v1()->y;
    triangleArray[0][2] = triangle.v1()->z;

    triangleArray[1][0] = triangle.v2()->x;
    triangleArray[1][1] = triangle.v2()->y;
    triangleArray[1][2] = triangle.v2()->z;

    triangleArray[2][0] = triangle.v3()->x;
    triangleArray[2][1] = triangle.v3()->y;
    triangleArray[2][2] = triangle.v3()->z;

    // Test if the triangle and the voxel are intersecting or not
    return checkTriangleBoxIntersection(voxelCenter, voxelHalfSize,
                                        triangleArray);
}

uint64_t Volume::_clampIndex(uint64_t idx, uint64_t dimension)
{
    idx = std::max(uint64_t(0) , idx);
    idx = std::min(idx, uint64_t(_grid->getDimension(dimension) - 1));
    return idx;
}

void Volume::_vec2grid(const Vector3f& v, GridIndex& grid)
{
    for (int32_t i = 0; i < DIMENSIONS; ++i)
    {
        grid[I2UI64(i)] = F2I64((v[i] - _meshOrigin[i]) / _voxelSize);
        grid[i] = _clampIndex(grid[i], i);
    }
}

void Volume::_getTriangleBoundingBox(AdvancedTriangle triangle,
                                     int64_t *tMin, int64_t *tMax)
{
    // Find the index of the voxel that intersects the triangle
    GridIndex vIdx;
    Vertex vertex(triangle.v1()->x, triangle.v1()->y, triangle.v1()->z);

    _vec2grid(vertex, vIdx);

    for (int64_t j = 0; j < DIMENSIONS; ++j)
    {
        tMin[j] = (vIdx[I2UI64(j)] - 1);
        tMax[j] = (vIdx[I2UI64(j)]);
    }

    for (uint64_t j = 1; j < DIMENSIONS; ++j)
    {
        if (j == 1)
        {
            vertex.x() = triangle.v2()->x;
            vertex.y() = triangle.v2()->y;
            vertex.z() = triangle.v2()->z;
        }
        else
        {
            vertex.x() = triangle.v3()->x;
            vertex.y() = triangle.v3()->y;
            vertex.z() = triangle.v3()->z;
        }

        _vec2grid(vertex, vIdx);

        for (uint64_t k = 0; k < DIMENSIONS; ++k)
        {
            if (vIdx[k] - 1 < (tMin[k]))
                tMin[k] = (vIdx[k] - 1);

            if (vIdx[k] > (tMax[k]))
                tMax[k] = (vIdx[k]);
        }
    }

    for (int32_t ii = 0; ii < DIMENSIONS ; ++ii)
    {
        tMin[ii] = std::max(int64_t(0), tMin[ii] - 1);
        tMax[ii] = std::min(int64_t(_grid->getDimension(ii) - 2), tMax[ii] + 2);
    }
}

void Volume::_getBoundingBox(Mesh* mesh, uint64_t i, int64_t *tMin, int64_t *tMax)
{
    // We need to get the bounding box of the voxels intersecting the given
    // triangle, so we get the triangle at index i and then we get the first
    // vertex of this triangle. Afterwards, we get the index of the voxel that
    // intersects this vertex
    GridIndex voxelIdx;
    _vec2grid(mesh->getVertices()[I2UI64(mesh->getTriangles()[i][0])], voxelIdx);

    // After getting the voxel index, we need to compute VOXEL bounding box
    // that surrounds the voxel (where the vertex is locted), simply
    for (int64_t j = 0; j < DIMENSIONS; ++j)
    {
        tMin[j] = voxelIdx[I2UI64(j)] - 1;
        tMax[j] = voxelIdx[I2UI64(j)];
    }

    for (int64_t j = 1; j < DIMENSIONS; ++j)
    {
        _vec2grid(mesh->getVertices()[(mesh->getTriangles()[i][j])], voxelIdx);

        for (uint64_t k = 0; k < DIMENSIONS; ++k)
        {
            if (voxelIdx[k] - 1 < tMin[k])
                tMin[k] = voxelIdx[k] - 1;

            if (voxelIdx[k] > tMax[k])
                tMax[k] = voxelIdx[k];
        }
    }

    for (int32_t ii = 0; ii < DIMENSIONS ; ++ii)
    {
        tMin[ii] = std::max(int64_t(0), tMin[ii] - 1);
        tMax[ii] = std::min(int64_t(_grid->getDimension(ii) - 1) , tMax[ii]);
    }
}

int64_t Volume::getWidth(void) const
{
    return _grid->getWidth();
}

int64_t Volume::getHeight(void) const
{
    return _grid->getHeight();
}

int64_t Volume::getDepth(void) const
{
    return _grid->getDepth();
}

uint64_t Volume::getNumberVoxels(void) const
{
    return I2UI64(_grid->getWidth() * _grid->getHeight() * _grid->getDepth());
}

double Volume::getSurfaceVoxelizationTime(void) const
{
    return _surfaceVoxelizationTime;
}

double Volume::getSolidVoxelizationTime(void) const
{
    return _solidVoxelizationTime;
}

double Volume::getAppendingVolumePassTime(void) const
{
    return _addingVolumePassTime;
}

bool  Volume::isFilled(const u_int64_t& index) const
{
    return _grid->isFilled(index);
}

bool  Volume::isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    return isFilled(mapToIndex(x, y, z));
}

uint64_t Volume::mapToIndex(const int64_t &x, const int64_t &y, const int64_t &z) const
{

    if(x >= getWidth()  || x < 0 || y >= getHeight() || y < 0 || z >= getDepth()  || z < 0)
    {
        LOG_ERROR("Index out of bound [%d, %d, %d]", x, y, z);
    }

    return I2UI64(x + (getWidth() * y) + (getWidth() * getHeight() * z));
}

void Volume::project(const std::string prefix,
                     const bool xy, const bool xz, const bool zy,
                     const bool &projectColorCoded) const
{
    _grid->projectVolume(prefix, xy, xz, zy, projectColorCoded);
}

void Volume::writeVolumes(const std::string &prefix,
                          const bool& binaryFormat,
                          const bool& rawFormat,
                          const bool &nrrdFormat) const
{
    // Start the timer
    TIMER_SET;

    if (binaryFormat || rawFormat || nrrdFormat)
        LOG_TITLE("Writing Volumes");

    if (binaryFormat)
    {
        LOG_STATUS("Bit Volume");
        _grid->writeBIN(prefix);
    }

    if (rawFormat)
    {
        LOG_STATUS("Byte Volume");
        _grid->writeRAW(prefix);
    }

    if (nrrdFormat)
    {
        LOG_STATUS("NRRD Volume");
        _grid->writeNRRD(prefix);
    }

    // Statictics
    LOG_STATUS_IMPORTANT("Writing Volumes Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStackXY(const std::string &outputDirectory, const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "xy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: XY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (int64_t z = 0; z < getDepth(); ++z)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getDepth());

        // Create a slice
        auto slice = std::make_unique<Image>(getWidth(), getHeight());

        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                uint64_t index = mapToIndex(I2I64(i), I2I64(j),
                                            I2I64(getDepth() - 1 - z));

                if (_grid->isFilled(index))
                    slice->setPixelColor(i , j, WHITE);
                else
                    slice->setPixelColor(i , j, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << z;
        slice->writePPM(stream.str());
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStackXZ(const std::string &outputDirectory,
                            const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "xz-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: XZ");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (int64_t y = 0; y < getHeight(); ++y)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getDepth());

        // Create a slice
        auto slice = std::make_unique<Image>(getWidth(), getHeight());

        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t k = 0; k < getDepth(); k++)
            {
                uint64_t index = mapToIndex(I2I64(i), I2I64(k),
                                            I2I64(getHeight() - 1 - y));

                if (_grid->isFilled(index))
                    slice->setPixelColor(i , k, WHITE);
                else
                    slice->setPixelColor(i , k, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << y;
        slice->writePPM(stream.str());
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStackZY(const std::string &outputDirectory,
                            const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "zy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: ZY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (int64_t i = 0; i < getWidth(); ++i)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        auto slice = std::make_unique<Image>(getDepth(), getHeight());

        for (int64_t z = 0; z < getDepth(); z++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                uint64_t index = mapToIndex(I2I64(getWidth() - 1 - i),
                                            I2I64(j), I2I64(z));

                if (_grid->isFilled(index))
                    slice->setPixelColor(z , getHeight() - j - 1, WHITE);
                else
                    slice->setPixelColor(z , getHeight() - j - 1, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << i;
        slice->writePPM(stream.str());
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::writeStacks(const std::string &outputDirectory,
                         const std::string &prefix,
                         const bool& xy,
                         const bool &xz,
                         const bool& zy) const
{
    // Start timer
    TIMER_SET;

    if (xy || xz || zy)
        LOG_TITLE("Writing Stacks");

    if (xy)
    {
        writeStackXY(outputDirectory, prefix);
    }

    if (xz)
    {
        writeStackXZ(outputDirectory, prefix);
    }

    if (zy)
    {
        writeStackZY(outputDirectory, prefix);
    }

    // Statictics
    LOG_STATUS_IMPORTANT("Writing Stacks Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Volume::exportToMesh(const std::string &prefix,
                          const bool &formatOBJ, const bool &formatPLY,
                          const bool &formatOFF, const bool &formatSTL)
{
    if (!(formatOBJ || formatPLY || formatOFF || formatSTL))
    {
        LOG_WARNING("Exporto mesh option must be enabled to export this mesh. "
                    "User one of the following: "
                    "[--export-obj, --export-ply, --export-off, --export-stl]");
        return;
    }

    TIMER_SET;
    LOG_TITLE("Constructing Volume Mesh");

    // The generated mesh from the volume
    std::unique_ptr< VolumeMesh > volumeMesh = std::make_unique< VolumeMesh >();

    // Delta value
    const Vector3f delta(1, 1, 1);

    LOOP_STARTS("Searching Filled Voxels")
    for (int64_t i = 0; i <  _grid->getWidth(); ++i)
    {
        LOOP_PROGRESS(i, _grid->getWidth());

        for (int64_t j = 0; j < _grid->getHeight(); ++j)
        {
            for (int64_t k = 0; k < _grid->getDepth(); ++k)
            {
                // Skip empty voxels
                if (_grid->isEmpty(i, j, k))
                    continue;

                Vector3f coordinate(0.5f + i, 0.5f + j, 0.5f + k);
                Vector3f pMin = _baseResolution * (coordinate - 0.5f * delta) + _meshOrigin;
                Vector3f pMax = _baseResolution * (coordinate + 0.5f * delta) + _meshOrigin;

                // A mesh representing the bounding box of the cube
                VolumeMesh* voxelCube = VolumeMesh::constructVoxelCube(pMin, pMax);

                // Append it to the volume mesh
                volumeMesh->append(voxelCube);

                // Free the voxel cube
                voxelCube->~VolumeMesh();
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_STATUS_IMPORTANT("Volume Mesh Construction Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    LOG_TITLE("Exporting Volume Mesh");
    TIMER_RESET;
    const std::string outputPrefix = prefix + VOLUME_MESH_SUFFIX;
    if (formatOBJ)
    {
        exportOBJ(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatPLY)
    {
        exportPLY(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatSTL)
    {
        exportSTL(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    if (formatOFF)
    {
        exportOFF(outputPrefix,
                  volumeMesh->vertices.data(), volumeMesh->vertices.size(),
                  volumeMesh->triangles.data(), volumeMesh->triangles.size());
    }

    LOG_STATUS_IMPORTANT("Exporting Volume Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

uint8_t Volume::getByte(const uint64_t index) const
{
    return _grid->getByte(index);
}

uint8_t Volume::getValue(const uint64_t index) const
{
    if (_grid->isFilled(index))
        return 255;
    else
        return 0;
}

uint8_t Volume::getValue(const int64_t &x,
                         const int64_t &y,
                         const int64_t &z) const
{
    if (_grid->isFilled(x, y, z))
        return 255;
    else
        return 0;
}

uint8_t Volume::getByte(const int64_t &x,
                        const int64_t &y,
                        const int64_t &z) const
{
    return _grid->getByte(mapToIndex(x, y, z));
}

void Volume::fillVoxel(const int64_t &x,
                       const int64_t &y,
                       const int64_t &z)
{

    if (x - 1 < 0 || x > _gridDimensions[0])
    {
        return;
    }
    if (y - 1 < 0 || y > _gridDimensions[1])
    {
        return;
    }
    if (z - 1 < 0 || z > _gridDimensions[2])
    {
        return;
    }

    _grid->fillVoxel(x - 1, y, z);
    _grid->fillVoxel(x - 1, y - 1, z);
    _grid->fillVoxel(x - 1, y + 1, z);
    _grid->fillVoxel(x - 1, y, z - 1);
    _grid->fillVoxel(x - 1, y - 1, z - 1);
    _grid->fillVoxel(x - 1, y + 1, z - 1);
    _grid->fillVoxel(x - 1, y, z + 1);
    _grid->fillVoxel(x - 1, y - 1, z + 1);
    _grid->fillVoxel(x - 1, y + 1, z + 1);

    // x = 0
    _grid->fillVoxel(x, y - 1, z);
    _grid->fillVoxel(x, y + 1, z);
    _grid->fillVoxel(x, y, z - 1);
    _grid->fillVoxel(x, y - 1, z - 1);
    _grid->fillVoxel(x, y + 1, z - 1);
    _grid->fillVoxel(x, y, z + 1);
    _grid->fillVoxel(x, y - 1, z + 1);
    _grid->fillVoxel(x, y + 1, z + 1);

    // x = +1
    _grid->fillVoxel(x + 1, y, z);
    _grid->fillVoxel(x + 1, y - 1, z);
    _grid->fillVoxel(x + 1, y + 1, z);
    _grid->fillVoxel(x + 1, y, z - 1);
    _grid->fillVoxel(x + 1, y - 1, z - 1);
    _grid->fillVoxel(x + 1, y + 1, z - 1);
    _grid->fillVoxel(x + 1, y, z + 1);
    _grid->fillVoxel(x + 1, y - 1, z + 1);
    _grid->fillVoxel(x + 1, y + 1, z + 1);

}

void Volume::addByte(const uint64_t &index, const uint8_t byte)
{
    _grid->addByte(index, byte);
}

void Volume::clear(void)
{
    _grid->clear();
}

void Volume::fill(const int64_t &x,
                  const int64_t &y,
                  const int64_t &z)
{
    _grid->fillVoxel(x, y, z);
}

void  Volume::fill(const u_int64_t& index)
{
    _grid->fillVoxel(index);
}

void Volume::clear(const int64_t &x,
                   const int64_t &y,
                   const int64_t &z)
{
    _grid->clearVoxel(x, y, z);
}

void  Volume::clear(const u_int64_t& index)
{
    _grid->clearVoxel(index);
}

uint64_t Volume::computeNumberNonZeroVoxels(void) const
{
    return _grid->computeNumberNonZeroVoxels();
}

std::string Volume::getFormatString() const
{
    return VolumeGrid::getTypeString(_gridType);
}

float Volume::computeVolume3()
{
    // Compute the number of non zero voxels
    const uint64_t numberNonZeroVoxels = _grid->computeNumberNonZeroVoxels();

    // Get the voxel volume in units3
    const float voxelVolume =
            _voxelSize * _voxelSize * _voxelSize;

    // Return the result
    return voxelVolume * numberNonZeroVoxels;
}

void Volume::printVolumeStats(const std::string &reference,
                              const std::string *prefix)
{
    LOG_TITLE("Volume Statistics");

    LOG_STATUS("Collecting Stats.");
    Vector3f bounds = _pMax - _pMin;
    float volumeSize = computeVolume3();

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + VOLUME_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
        {
            fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                     F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
            fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                     F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
            fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                     F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
        }

        fprintf(info, "\t* Resolution            | [%d] x [%d] x [%d] \n",
                 I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
        fprintf(info, "\t* Number of Voxels      | %" PRIu64 " \n",
                getNumberVoxels());
        fprintf(info, "\t* Volume Format         | %s \n",
                 getFormatString().c_str());
        fprintf(info, "\t* Size in Memory        | %sBytes \n",
                 FORMAT(getNumberBytes()));
        fprintf(info, "\t* Volume                | %f³ \n",
                 F2D(volumeSize));

        // Close the file
        fclose(info);
    }

    LOG_STATUS("Volume [ %s ]", reference.c_str());

    if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
    {
        LOG_INFO("\t* Bounding Box:         | [%f, %f, %f]",
                 F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
        LOG_INFO("\t* pMin:                 | [%f, %f, %f]",
                 F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
        LOG_INFO("\t* pMax:                 | [%f, %f, %f]",
                 F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
    }

    LOG_INFO("\t* Resolution            | [%d] x [%d] x [%d]",
             I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
    LOG_INFO("\t* Number of Voxels      | %" PRIu64 "",
             getNumberVoxels());
    LOG_INFO("\t* Volume Format         | %s",
             getFormatString().c_str());
    LOG_INFO("\t* Size in Memory        | %sBytes",
             FORMAT(getNumberBytes()));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volumeSize));
}

void Volume::addVolumePass(const Volume* volume)
{
    // Start the timer
    TIMER_SET;

    // If the volume dimensions are not similar to this one
    // then, print an ERROR ...
    if (!(getWidth() == volume->getWidth() &&
          getHeight() == volume->getHeight() &&
          getDepth() == volume->getDepth()))
    {
        LOG_ERROR("The dimensions of the two volumes don't match");
        return;
    }

    // Loop over the volume elements byte-by-byte and add
    // them to the corresponding one in this voume
    size_t progress = 0;
    // #pragma omp parallel for schedule(guided, 1)
    for (size_t i = 0; i < volume->getNumberBytes(); ++i)
    {
        addByte(i, volume->getByte(i));
    }

    _addingVolumePassTime = GET_TIME_SECONDS;
}

void Volume::addVolume(const std::string &volumePrefix)
{
    Volume* volume = new Volume(volumePrefix, VolumeGrid::TYPE::BYTE);

    // If the volume dimensions are not similar to this one
    // then, print an ERROR ...
    if (!(getWidth() == volume->getWidth() &&
          getHeight() == volume->getHeight() &&
          getDepth() == volume->getDepth()))
    {
        LOG_ERROR("The dimensions of the two volumes don't match");
        return;
    }

    // Loop over the volume elements byte-by-byte and add
    // them to the corresponding one in this voume
    for (size_t i = 0; i < volume->getNumberBytes(); ++i)
    {
        addByte(i, volume->getByte(i));
    }

    // Release the gird to save some memory
    volume->~Volume();
}

uint64_t Volume::getNumberBytes(void) const
{
    return _grid->getNumberBytes();
}

int32_t Volume::getLargestDimension(const Vector3f& dimensions)
{
    float value = dimensions[0];
    int index = 0;

    for (int32_t i = 1; i < DIMENSIONS; ++i)
    {
        if (value < dimensions[i])
        {
            index = i;
            value = dimensions[i];
        }
    }

    return index;
}

Volume::~Volume()
{
    _grid->~VolumeGrid();
}

Volume* Volume::constructIsoValueVolume(const Volume* volume,
                                        const uint8_t& isoValue,
                                        const int64_t &padding)
{
    Vector3f pMin(0.f), pMax(1.f);

    pMax.x() *= volume->getWidth();
    pMax.y() *= volume->getHeight();
    pMax.z() *= volume->getDepth();

    Volume* isoVolume = new Volume(volume->getWidth() + padding,
                                   volume->getHeight() + padding,
                                   volume->getDepth() + padding,
                                   pMin, pMax, VolumeGrid::TYPE::BIT);

    LOG_STATUS("Constructing Iso Volume");
    for (int64_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (int64_t y = 0; y < volume->getHeight(); ++y)
        {
            for (int64_t z = 0; z < volume->getDepth(); ++z)
            {
                int64_t xIso = x + F2I64(padding / 2.f);
                int64_t yIso = y + F2I64(padding / 2.f);
                int64_t zIso = z + F2I64(padding / 2.f);

                if (volume->getByte(x, y, z) == isoValue)
                    isoVolume->fill(xIso, yIso, zIso);
                else
                    isoVolume->clear(xIso, yIso, zIso);
            }
        }
    }
    LOOP_DONE;

    return isoVolume;
}

std::vector<uint64_t> Volume::createHistogram(const Volume* volume)
{
    // Allocation
    std::vector<uint64_t> histogram;
    histogram.resize(256);

    // Initialization
    for (uint64_t i = 0; i < 256; ++i)
        histogram[i] = 0;

    for (int64_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (int64_t y = 0; y < volume->getHeight(); ++y)
        {
            for (int64_t z = 0; z < volume->getDepth(); ++z)
            {
                histogram[volume->getByte(x, y, z)] += 1;
            }
        }
    }

    // Return the histogram array
    return histogram;
}

Volume* Volume::constructFullRangeVolume(const Volume* volume,
                                         const int64_t &padding)
{
    Vector3f pMin(0.f), pMax(1.f);

    pMax.x() *= volume->getWidth();
    pMax.y() *= volume->getHeight();
    pMax.z() *= volume->getDepth();

    Volume* isoVolume = new Volume(volume->getWidth() + padding,
                                   volume->getHeight() + padding,
                                   volume->getDepth() + padding,
                                   pMin, pMax, VolumeGrid::TYPE::BYTE);

    LOG_STATUS("Constructing Iso Volume");
    for (int64_t x = 0; x < volume->getWidth(); ++x)
    {
        LOOP_PROGRESS(x, volume->getWidth());
        for (int64_t y = 0; y < volume->getHeight(); ++y)
        {
            for (int64_t z = 0; z < volume->getDepth(); ++z)
            {
                int64_t xIso = x + F2I64(padding / 2.f);
                int64_t yIso = y + F2I64(padding / 2.f);
                int64_t zIso = z + F2I64(padding / 2.f);

                if (volume->getByte(x, y, z) > 0)
                    isoVolume->fill(xIso, yIso, zIso);
                else
                    isoVolume->clear(xIso, yIso, zIso);

            }
        }
    }
    LOOP_DONE;

    return isoVolume;
}

Volume* Volume::constructFromTiffMask(
        const std::string &maskDirectory,
        const int64_t &maskWidth, const int64_t &maskHeight,
        const Ultraliser::VolumeGrid::TYPE& gridType)
{
    // Set the timer
    TIMER_SET;

    LOG_TITLE("Loading .Tiff");
    LOG_STATUS("Mask Directory [ %s ]", maskDirectory.c_str());

    // Get a list of all the stacks of the masks
    std::vector< std::string > maskFiles;
    Ultraliser::Directory::list(maskDirectory, maskFiles, ".tif");

    // Sort the mask files
    std::sort(maskFiles.begin(), maskFiles.end());

    // Adding a little delta
    int64_t numZeroPaddingVoxels = 16;
    Volume* maskVolume = new Volume(
                maskWidth + numZeroPaddingVoxels,
                maskHeight + numZeroPaddingVoxels,
                I2I64(maskFiles.size()) + numZeroPaddingVoxels,
                Vector3f(),
           		Vector3f(),
                gridType);
    LOG_INFO("%d %d %d", maskVolume->getWidth(), maskVolume->getHeight(), maskVolume->getDepth());

    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for(int64_t i = 0; i < I2I64(maskFiles.size()); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, maskFiles.size());

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

        // Read the image
        std::string imagePath = maskDirectory + "/" + maskFiles[I2UI64(i)];
        std::unique_ptr< TiffImage > image(new TiffImage);
        image->setimageFile(imagePath);
        image->readImage();

        // Update the volume
        for(int64_t x = 0; x < maskWidth; ++x)
        {
            for(int64_t y = 0; y < maskHeight; ++y)
            {
                if(image->isPixelFilled(I2I32(x), I2I32(y)))
                {
                    maskVolume->fill(x + numZeroPaddingVoxels / 2,
                                     y + numZeroPaddingVoxels / 2,
                                     i + numZeroPaddingVoxels / 2);
                }
            }
        }
    }
    LOOP_DONE;

    // Statictics
    LOG_STATS(GET_TIME_SECONDS);

    // Return a pointer to the mask volume
    return maskVolume;
}

Volume::SOLID_VOXELIZATION_AXIS Volume::getSolidVoxelizationAxis(const std::string &argumentString)
{
    if (argumentString == "x")
    {
        return SOLID_VOXELIZATION_AXIS::X;
    }
    else if (argumentString == "y")
    {
        return SOLID_VOXELIZATION_AXIS::Y;
    }
    else if (argumentString == "z")
    {
        return SOLID_VOXELIZATION_AXIS::Z;
    }
    else if (argumentString == "xyz")
    {
        return SOLID_VOXELIZATION_AXIS::XYZ;
    }
    else
    {
        // Error, therefore terminate
        LOG_ERROR("The option [ %s ] is not valid for --solid-voxelization-axis! "
                  "Please use one of the following [x, y, z, xyz].", argumentString.c_str());

        // For the sake of compilation only.
        return SOLID_VOXELIZATION_AXIS::XYZ;
    }
}

}
