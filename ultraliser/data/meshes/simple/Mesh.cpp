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

#include "Mesh.h"
#include "MeshOperations.h"
#include "MeshStatistics.h"
#include "TriangleOperations.h"
#include <utilities/Utilities.h>
#include <math/Math.h>
#include <geometry/Intersection.h>
#include <algorithms/SectionGeometry.h>

#define ANGLE_ERROR 0.123456789f

#define MAXIMUM_FLOAT_VALUE  99999.f
#define MINIMUM_FLOAT_VALUE -99999.f

#define MAXIMUM_DOUBLE_VALUE  99999.0
#define MINIMUM_DOUBLE_VALUE -99999.0

#define VERTEX_DELETION_VALUE -998.3824223883588f

namespace Ultraliser
{

Mesh::Mesh()
{
    // Initialization of public member variables
    _numberVertices = 0;
    _numberTriangles = 0;
    _optimizationTime = 0.0;

    _vertices = nullptr;
    _triangles = nullptr;

    _neighbors = nullptr;
    _neighborList = nullptr;
}

Mesh::Mesh(const uint64_t &numVertices,
           const uint64_t &numTriangles)
{
    // Initialization of public member variables

    _numberVertices = numVertices;
    _numberTriangles = numTriangles;
    _optimizationTime = 0.0;

    _vertices = new Vertex[_numberVertices];
    _triangles = new Triangle[_numberTriangles];

    _neighbors = nullptr;
    _neighborList = nullptr;
}

Mesh::Mesh(const Samples& samples)
{
    SectionGeometry sg(samples, 16);
    _numberVertices = sg.numVertices;
    _vertices = sg.vertices;
    _numberTriangles = sg.numTriangles;
    _triangles = sg.triangles;
}

Mesh::Mesh(Vertices vertices, Triangles triangles)
{
    _numberVertices = vertices.size();
    _numberTriangles = triangles.size();
    _optimizationTime = 0.0;

    this->_vertices = new Vertex[_numberVertices];
    this->_triangles = new Triangle[_numberTriangles];

    const Vertex* vertexData = vertices.data();
    const Triangle* triangleData = triangles.data();


    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        this->_vertices[i] = vertexData[i];
    }

    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
        this->_triangles[i] = triangleData[i];
    }

    _neighbors = nullptr;
    _neighborList = nullptr;
}

Mesh::Mesh(const std::string &fileName, const bool& verbose)
{
    // Import the mesh from a given file
    import(fileName, verbose);

    _optimizationTime = 0.0;

    _neighbors = nullptr;
    _neighborList = nullptr;
}

const Vertex* Mesh::getVertices() const
{
    return _vertices;
}

const Triangle* Mesh::getTriangles() const
{
    return _triangles;
}

void Mesh::getTriangleBoundingBox(const uint64_t& triangleIndex,
                                  Vector3f& pMin, Vector3f& pMax) const
{
    pMin = Vector3f(std::numeric_limits<float>::max());
    pMax = Vector3f(std::numeric_limits<float>::min());

    Triangle triangle = getTriangles()[triangleIndex];
    for (uint64_t vertexIndex = 0; vertexIndex < 3; ++vertexIndex)
    {
        Vector3f vertex = getVertices()[triangle[vertexIndex]];
        for (int iDim = 0; iDim < DIMENSIONS; ++iDim)
        {
            if (vertex[iDim] < pMin[iDim])
                pMin[iDim] = vertex[iDim];

            if (vertex[iDim] > pMax[iDim])
                pMax[iDim] = vertex[iDim];
        }
    }
}

void Mesh::computeBoundingBox(Vector3f& pMinIn, Vector3f& pMaxIn)
{
    // Create new variables to avoid any mess if the inputs are already
    // initialized with some values
    Vector3f pMin(std::numeric_limits<float>::max());
    Vector3f pMax(std::numeric_limits<float>::min());

    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        Vertex vertex = _vertices[i];
        for (int iDim = 0; iDim < DIMENSIONS; ++iDim)
        {
            if (vertex[iDim] < pMin[iDim])
                pMin[iDim] = vertex[iDim];

            if (vertex[iDim] > pMax[iDim])
                pMax[iDim] = vertex[iDim];
        }
    }

    // Update the bounding box data
    pMinIn = pMin;
    pMaxIn = pMax;
}

void Mesh::computeRelaxedBoundingBox(Vector3f& pMinIn,
                                     Vector3f& pMaxIn,
                                     const float& relaxationPercentage)
{
    // Create new variables to avoid any mess if the inputs are already
    // initialized with some values
    Vector3f pMin(std::numeric_limits<float>::max());
    Vector3f pMax(std::numeric_limits<float>::min());

    // Compute the actual bounding box
    computeBoundingBox(pMin, pMax);

    // Extend the bounding box a little bit to avoid edge issues
    Ultraliser::Vector3f boundingBoxSize = pMax - pMin;

    // Compute the largest dimension
    float largestDimension = boundingBoxSize[0];
    if (boundingBoxSize[1] > largestDimension)
        largestDimension = boundingBoxSize[1];
    if (boundingBoxSize[2] > largestDimension)
        largestDimension = boundingBoxSize[2];

    float voxelSize = largestDimension / relaxationPercentage;

    // Stretching
    const Vector3f relaxation(5 * voxelSize);
    pMinIn = pMin - relaxation;
    pMaxIn = pMax + relaxation;
}

void Mesh::centerAtOrigin(void)
{
    // Compute the bounding box
    Vector3f pMin, pMax;
    computeBoundingBox(pMin, pMax);

    // Compute the center
    Vector3f boundingBoxSize = pMax - pMin;
    Vector3f center = pMin + (0.5 * boundingBoxSize);

    // Shift the vertices to the origin to center the mesh
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        _vertices[i] -= center;
    }
}

void Mesh::rotate(const Matrix4f& matrix)
{
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        Vector4f result = matrix * Vector4f(_vertices[i]);
        _vertices[i] = Vector3f(result.x(), result.y(), result.z());
    }
}

void Mesh::transform(const Matrix4f& matrix)
{
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        Vector4f result = matrix * Vector4f(_vertices[i], 1.0);
        _vertices[i] = Vector3f(result.x(), result.y(), result.z());
    }
}

void Mesh::translate(const Vector3f& to)
{
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        _vertices[i] += to;
    }
}

void Mesh::uniformScale(const float factor)
{
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        _vertices[i] *= factor;
    }
}

void Mesh::scale(const float x, const float y, const float z)
{
    for (uint64_t i = 0; i < _numberVertices; ++i)
    {
        _vertices[i].x() *= x;
        _vertices[i].y() *= y;
        _vertices[i].z() *= z;
    }
}

void Mesh::import(const std::string &fileName, const bool &verbose)
{
    // Start the timer
    TIMER_SET;

    std::ifstream fileStream;
    fileStream.open(fileName.c_str());
    if (!fileStream.good())
    {
        LOG_ERROR("Cannot load mesh file [ %s ]", fileName.c_str());
    }
    fileStream.close();

    if (verbose) LOG_TITLE("Importing Mesh");
    Vertices loadedVertices;
    Triangles loadedTriangles;

    // Switch to the corresponding loader based on the extension
    switch(fileName[strlen(fileName.c_str()) - 1])
    {
    case 'y':
    case 'Y':
    {
        importPLY(fileName, loadedVertices, loadedTriangles, verbose);
        break;
    }
    case 'j':
    case 'J':
    {
        importOBJ(fileName, loadedVertices, loadedTriangles, verbose);
        break;
    }
    case 'l':
    case 'L':
    {
        importSTL(fileName, loadedVertices, loadedTriangles, verbose);
        break;
    }   
    case 'f':
    case 'F':
    {
        importOFF(fileName, loadedVertices, loadedTriangles, verbose);
        break;
    }
    default:
        LOG_ERROR("The mesh [ %s ] CANNOT be read.", fileName.c_str());
    }

    _numberVertices = loadedVertices.size();
    _numberTriangles = loadedTriangles.size();

    _vertices = new Vertex[_numberVertices];
    _triangles = new Triangle[_numberTriangles];

    Vector3f* tempVertices = loadedVertices.data();
    Triangle* tempTriangles = loadedTriangles.data();

    for (uint64_t i = 0; i < _numberVertices; ++i)
        _vertices[i] = tempVertices[i];

    for (uint64_t i = 0; i < _numberTriangles; ++i)
        _triangles[i] = tempTriangles[i];

    // Statistics
    if (verbose) LOG_STATUS_IMPORTANT("Importing Mesh Stats.");
    if (verbose) LOG_STATS(GET_TIME_SECONDS);
}

Mesh* Mesh::instanciate(const uint64_t &numVertices,
                        const uint64_t &numTriangles)
{
    // Initialization of public member variables
    Mesh* instance = new Mesh(numVertices,
                              numTriangles);

    // Return a pointer to the created mesh
    return instance;
}


double Mesh::getDefaultOptimizationTime() const
{
    return _optimizationTime;
}

uint64_t Mesh::getNumberVertices() const
{
    return _numberVertices;
}

uint64_t Mesh::getNumberTriangles() const
{
    return _numberTriangles;
}

float Mesh::_cotangentAngle(const Vector3f& pivot,
                            const Vector3f& a,
                            const Vector3f& b)
{
    const auto pA = (a - pivot).normalized();
    const auto pB = (b - pivot).normalized();

    const auto sinA = Vector3f::cross(pA, pB).abs();
    const auto cosA = Vector3f::dot(pA, pB);

    return cosA / sinA;
}

void Mesh::_computeNeighborhoods(Mesh& mesh,
                                 Neighborhood& vertexNeighbors,
                                 Neighborhood& faceNeighbors)
{
    vertexNeighbors.resize(mesh.getNumberVertices());
    faceNeighbors.resize(mesh.getNumberVertices());

    for(uint64_t i = 0; i < mesh.getNumberTriangles(); ++i)
    {
        const auto& face = mesh._triangles[i];

        vertexNeighbors[I2UI64(face.x())].insert(I2UI64(face.y()));
        vertexNeighbors[I2UI64(face.x())].insert(I2UI64(face.z()));

        vertexNeighbors[I2UI64(face.y())].insert(I2UI64(face.x()));
        vertexNeighbors[I2UI64(face.y())].insert(I2UI64(face.z()));

        vertexNeighbors[I2UI64(face.z())].insert(I2UI64(face.x()));
        vertexNeighbors[I2UI64(face.z())].insert(I2UI64(face.y()));

        faceNeighbors[I2UI64(face.x())].insert(i);
        faceNeighbors[I2UI64(face.y())].insert(i);
        faceNeighbors[I2UI64(face.z())].insert(i);
    }
}

Vector3f Mesh::_computeKernel(Mesh& mesh,
                              const uint64_t vertexIndex,
                              Neighbors& verticesN,
                              Neighbors& facesN)
{
    auto sumWeights = 0.f;
    Vector3f sumPos(0.f);

    const auto& currentVertex = mesh._vertices[vertexIndex];

    for(const auto& nv : verticesN)
    {
        if(nv == vertexIndex)
            continue;

        const auto weight = _computeCotangentWeight(mesh, vertexIndex,
                                                    nv, facesN);
        sumWeights += weight;
        sumPos += mesh._vertices[nv] * weight;
    }

    const auto kernel = (sumPos / sumWeights) - currentVertex;

    return kernel;
}

float Mesh::_computeCotangentWeight(Mesh& mesh,
                                    const uint32_t vertexIndex,
                                    const uint32_t neighborIndex,
                                    Neighbors& faceN)
{
    uint64_t edge1 = 0, edge2 = 0;
    bool firstFound = false, secondFound = false;
    for(const auto & nvFace : faceN)
    {
        const Triangle& f = mesh._triangles[nvFace];
        if(_triangleContainsVertex(f, vertexIndex) &&
           _triangleContainsVertex(f, neighborIndex))
        {
            if(!firstFound)
            {
                firstFound = true;
                edge1 = f.x() != vertexIndex &&
                        f.x() != neighborIndex ? f.x() : (f.y() != vertexIndex &&
                        f.y() != neighborIndex ? f.y() : f.z());
            }
            else if(!secondFound)
            {
                secondFound = true;
                edge2 = f.x() != vertexIndex &&
                        f.x() != neighborIndex ? f.x() : (f.y() != vertexIndex &&
                        f.y() != neighborIndex ? f.y() : f.z());
            }
        }

        if(firstFound && secondFound)
            break;
    }

    if(!firstFound || !secondFound)
    {
        throw std::runtime_error("LaplacianSmoothOperator: Mesh has boundaries");
    }

    const auto& pivotA = (mesh._vertices[edge1]);
    const auto& pivotB = (mesh._vertices[edge2]);

    const auto& a = (mesh._vertices[vertexIndex]);
    const auto& b = (mesh._vertices[neighborIndex]);

    const auto ctgA = _cotangentAngle(pivotA, a, b);
    const auto ctgB = _cotangentAngle(pivotB, a, b);

    return 0.5f * (ctgA + ctgB);
}

Vertex Mesh::_smoothVertex(Mesh& mesh,
                           const uint64_t vertexIndex,
                           const Vector3f kernel,
                           const float param)
{
    auto pos = mesh._vertices[vertexIndex];
    pos += (param * kernel);
    return pos;
}

void Mesh::scaleAndTranslate(const Vector3f &center, const Vector3f &BoundingBox)
{
    // Center the reconstructed mesh at the origin
    centerAtOrigin();

    // Compute the bounding box of the mesh
    Ultraliser::Vector3f pMaxGenerated, pMinGenerated;
    computeBoundingBox(pMinGenerated, pMaxGenerated);

    // Compute the scale needed
    const Ultraliser::Vector3f currentBoundingBox = pMaxGenerated - pMinGenerated;
    const Ultraliser::Vector3f scaleValue = BoundingBox / currentBoundingBox;

    // Scale the mesh
    scale(scaleValue.x(), scaleValue.y(), scaleValue.z());

    // Translate it back to the original center of the input mesh
    translate(center);
}

void Mesh::applyLaplacianSmooth(const uint32_t& numIterations,
                                const float& smoothLambda,
                                const float& inflateMu)
{
    // If no iterations, return
    if (numIterations < 1)
        return;

    LOG_TITLE("Laplacian Smoothing");
    TIMER_SET;

    Neighborhood vertexN, faceN;
    _computeNeighborhoods(*this, vertexN, faceN);

    std::vector< Vertex > smoothedVertices(_numberVertices);
    LOOP_STARTS("Smoothing")
    for(uint32_t i = 0; i < numIterations; ++i)
    {
        #pragma omp parallel for
        for(uint64_t v = 0; v < _numberVertices; ++v)
        {
            const Vector3f kernel =
                    _computeKernel(*this, v, vertexN[v], faceN[v]);

            smoothedVertices[v] =
                    _smoothVertex(*this, v, kernel, smoothLambda);
        }

        // Update vertices
        #pragma omp parallel for
        for(uint64_t v = 0; v < _numberVertices; ++v)
        {
            _vertices[v] = smoothedVertices[v];
        }

        #pragma omp parallel for
        for(uint64_t v = 0; v < _numberVertices; ++v)
        {
            const Vector3f kernel =_computeKernel(*this, v, vertexN[v], faceN[v]);

            smoothedVertices[v] =
                    _smoothVertex(*this, v, kernel, inflateMu);
        }

        // Update vertices
        #pragma omp parallel for
        for(uint64_t v = 0; v < _numberVertices; ++v)
        {
            _vertices[v] = smoothedVertices[v];
        }

        LOOP_PROGRESS(i, numIterations);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clear the vector
    smoothedVertices.clear();

    LOG_STATUS_IMPORTANT("Laplacian Operator Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::smoothLaplacian(uint64_t numIterations)
{
    // If no iterations, return
    if (numIterations < 1)
        return;

    LOG_TITLE("Laplacian Smoothing");
    TIMER_SET;

    Neighborhood vertexN, faceN;
    _computeNeighborhoods(*this, vertexN, faceN);

    // A new list of vertices
    auto newVertices = std::make_unique< Vertex[] >(_numberVertices);

    LOOP_STARTS("Smoothing")
    for (uint64_t i = 0; i < numIterations; ++i)
    {
        #pragma omp parallel for
        for (uint64_t iVertex = 0; iVertex < _numberVertices; ++iVertex)
        {
            // New vertex
            Vector3f newVertex = _vertices[iVertex];

            // Neighbours weights
            auto neighbours = vertexN[iVertex];
            for (auto neighbour : neighbours)
            {
                newVertex += _vertices[neighbour];
            }

            // Update the new vertex
            newVertices[iVertex] = newVertex / (neighbours.size() + 1);
        }

        for (size_t iVertex = 0; iVertex < _numberVertices; ++iVertex)
        {
            Vertex tmpVertex;
            tmpVertex = newVertices[iVertex];

            newVertices[iVertex] = _vertices[iVertex];
            _vertices[iVertex] = tmpVertex;
        }

        LOOP_PROGRESS(i, numIterations);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Update vertices
    #pragma omp parallel for
    for (size_t iVertex = 0; iVertex < _numberVertices; ++iVertex)
    {
        _vertices[iVertex] = newVertices[iVertex];
    }

    // Clean
    vertexN.clear();
    vertexN.shrink_to_fit();
    faceN.clear();
    faceN.shrink_to_fit();

    LOG_STATUS_IMPORTANT("Laplacian Operator Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::exportMesh(const std::string &prefix,
                      const bool &formatOBJ,
                      const bool &formatPLY,
                      const bool &formatOFF,
                      const bool &formatSTL) const
{
    if (formatOBJ || formatPLY || formatOFF || formatSTL)
    {
        // Start timer
        TIMER_SET;

        LOG_TITLE("Exporting Mesh");

        if (formatOBJ)
        {
            exportOBJ(prefix,
                      _vertices, _numberVertices, _triangles, _numberTriangles);
        }

        if (formatPLY)
        {
            exportPLY(prefix, _vertices, _numberVertices, _triangles, _numberTriangles);
        }


        if (formatOFF)
        {
            exportOFF(prefix,
                      _vertices, _numberVertices, _triangles, _numberTriangles);
        }

        if (formatSTL)
        {
            exportSTL(prefix,
                      _vertices, _numberVertices, _triangles, _numberTriangles);
        }

        LOG_STATUS_IMPORTANT("Exporting Mesh Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void Mesh::writeDistributions(const std::string &reference, const std::string *prefix) const
{
    LOG_TITLE("Mesh Distributions");

    // Write the distributions
    MeshStatistics stats (_vertices, _triangles, _numberVertices, _numberTriangles);
    stats.writeStatsDistributions(*prefix + "-" + reference);
}

void Mesh::printStats(const std::string &reference,
                          const std::string *prefix) const
{
    LOG_TITLE("Mesh Statistics");

    LOG_STATUS("Collecting Stats.");

    float area = computeMeshSurfaceArea(_vertices, _triangles, _numberTriangles);
    float volume = computeMeshVolume(_vertices, _triangles, _numberTriangles);

    Vector3f pMin, pMax;
    computeMeshBoundingBox(_vertices, _numberVertices, pMin, pMax);
    Vector3f bounds = pMax - pMin;

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + MESH_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
        fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                F2D(pMin.x()), F2D(pMin.y()), F2D(pMin.z()));
        fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                F2D(pMax.x()), F2D(pMax.y()), F2D(pMax.z()));
        fprintf(info, "\t* Number Vertices       | %s \n",
                FORMAT(_numberVertices));
        fprintf(info, "\t* Number Triangles      | %s \n",
                FORMAT(_numberTriangles));
        fprintf(info, "\t* Surface Area          | %f² \n",
                F2D(area));
        fprintf(info, "\t* Volume                | %f³ \n",
                F2D(volume));

        if (_optimizationTime > 0.0)
            fprintf(info, "\t* Optimization Time     | %f Seconds",
                    _optimizationTime);

        // Close the file
        fclose(info);
    }

    LOG_STATUS_IMPORTANT("Mesh Stats. [ %s ]", reference.c_str());
    LOG_INFO("\t* Bounding Box:         | [%f, %f, %f]",
             F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
    LOG_INFO("\t* pMin:                 | [%f, %f, %f]",
             F2D(pMin.x()), F2D(pMin.y()), F2D(pMin.z()));
    LOG_INFO("\t* pMax:                 | [%f, %f, %f]",
             F2D(pMax.x()), F2D(pMax.y()), F2D(pMax.z()));
    LOG_INFO("\t* Number Vertices       | %s",
             FORMAT(_numberVertices));
    LOG_INFO("\t* Number Triangles      | %s",
             FORMAT(_numberTriangles));
    LOG_INFO("\t* Surface Area          | %f²",
             F2D(area));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volume));

    if (_optimizationTime > 0.0)
        LOG_INFO("\t* Optimization Time     | %f Seconds",
                 _optimizationTime);
}

void Mesh::_releaseData()
{
    // Free allocated memory
    delete [] _vertices;
    delete [] _triangles;

    // Destroy neighbor_list
    destroyNeighborlist();
}

Mesh::~Mesh()
{
    _releaseData();
}

}
