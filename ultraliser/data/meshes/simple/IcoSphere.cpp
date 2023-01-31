#include "IcoSphere.h"

namespace Ultraliser
{

/// Number of vertices and triangles at subdivision level of 1
#define VERTEX_COUNT 12
#define TRIANGLE_COUNT 20

IcoSphere::IcoSphere(const size_t& subdivisions)
{
    _subdivisions = subdivisions;

    // The vertices array for subdivision level of 1
    float vertices [VERTEX_COUNT * 3] =
    {
         0.000000,  0.000000, -1.000000,
         0.723600, -0.525720, -0.447215,
        -0.276385, -0.850640, -0.447215,
        -0.894425,  0.000000, -0.447215,
        -0.276385,  0.850640, -0.447215,
         0.723600,  0.525720, -0.447215,
         0.276385, -0.850640,  0.447215,
        -0.723600, -0.525720,  0.447215,
        -0.723600,  0.525720,  0.447215,
         0.276385,  0.850640,  0.447215,
         0.894425,  0.000000,  0.447215,
         0.000000,  0.000000,  1.000000
    };

    // The triangles array for subdivision level of 1
    size_t triangles [TRIANGLE_COUNT * 3] =
    {
        0, 1, 2,
        1, 0, 5,
        0, 2, 3,
        0, 3, 4,
        0, 4, 5,
        1, 5, 10,
        2, 1, 6,
        3, 2, 7,
        4, 3, 8,
        5, 4, 9,
        1, 10, 6,
        2, 6, 7,
        3, 7, 8,
        4, 8, 9,
        5, 9, 10,
        6, 10, 11,
        7, 6, 11,
        8, 7, 11,
        9, 8, 11,
        10, 9, 11,
    };

    _vertices = new Vertex[VERTEX_COUNT];
    _triangles = new Triangle[TRIANGLE_COUNT];

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < VERTEX_COUNT; ++i)
    {
        _vertices[i] = Vertex(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
    }

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < TRIANGLE_COUNT;++i)
    {
        _triangles[i] = Triangle(triangles[3 * i], triangles[3 * i + 1], triangles[3 * i + 2]);
    }

    _numberVertices = VERTEX_COUNT;
    _numberTriangles = TRIANGLE_COUNT;

    if (_subdivisions > 1)
    {
        for (size_t i = 2; i <= subdivisions; ++i)
        {
            _subdivideTrianglesAtMidPointsAndNormalize();
        }
    }
}

void IcoSphere::_subdivideTriangleAtMidPointsAndNormalize(const size_t &triangleIndex,
                                                          std::vector< Vector3f >& vertexList,
                                                          std::vector< Triangle >& triangleList)
{
    // Get the triangle and the indices of its vertices
    const Triangle& t = _triangles[triangleIndex];
    const size_t& v0Index = t[0];
    const size_t& v1Index = t[1];
    const size_t& v2Index = t[2];

    // Geth the vertices of the triangles
    const Vector3f& v0 = _vertices[v0Index];
    const Vector3f& v1 = _vertices[v1Index];
    const Vector3f& v2 = _vertices[v2Index];

    // Compute the midpoints along each edge
    Vector3f v3 = 0.5f * (v0 + v1); // v01
    Vector3f v4 = 0.5f * (v1 + v2); // v12
    Vector3f v5 = 0.5f * (v2 + v0); // v20

    // Ensure that all the edges are located along the unit sphere, and add them to the list
    _ensureVertexOnUnitSphere(v3); vertexList.push_back(v3);
    _ensureVertexOnUnitSphere(v4); vertexList.push_back(v4);
    _ensureVertexOnUnitSphere(v5); vertexList.push_back(v5);

    // Get the new vertex index
    const size_t v3Index = _numberVertices + vertexList.size() - 3;
    const size_t v4Index = _numberVertices + vertexList.size() - 2;
    const size_t v5Index = _numberVertices + vertexList.size() - 1;

    /// Create the new triangles
    Triangle t0, t1, t2, t3;

    // t0 -> v0 v3 v5
    t0[0] = v0Index; t0[1] = v3Index; t0[2] = v5Index;
    triangleList.push_back(t0);

    // t1-> v3 v1 v4
    t1[0] = v3Index; t1[1] = v1Index; t1[2] = v4Index;
    triangleList.push_back(t1);

    // t2-> v3 v4 v5
    t2[0] = v3Index; t2[1] = v4Index; t2[2] = v5Index;
    triangleList.push_back(t2);

    // t3-v4 v2 v5
    t3[0] = v4Index; t3[1] = v2Index; t3[2] = v5Index ;
    triangleList.push_back(t3);
}

void IcoSphere::_subdivideTrianglesAtMidPointsAndNormalize()
{
    // New lists of vertices and triangles
    std::vector< Vector3f > createdVertices;
    std::vector< Triangle > createdTriangles;

    // Subdivide all the triangles
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        _subdivideTriangleAtMidPointsAndNormalize(i, createdVertices, createdTriangles);
    }

    // Update the _vertices and _triangles lists
    const size_t totalNumberVertices = _numberVertices + createdVertices.size();
    const size_t totalNumberTriangles = createdTriangles.size();

    // Allocate the new arrays
    Vector3f* newVertices = new Vector3f[totalNumberVertices];
    Triangle* newTriangles = new Triangle[totalNumberTriangles];

    // Copy the old vertices to the new vertices list
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        newVertices[i] = _vertices[i];
    }

    // Release the old vertices list
    delete [] _vertices; _vertices = nullptr;

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < createdVertices.size(); ++i)
    {
        newVertices[_numberVertices + i] = createdVertices[i];
    }

    createdVertices.clear(); createdVertices.shrink_to_fit();

    delete [] _triangles; _triangles = nullptr;

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < createdTriangles.size(); ++i)
    {
        newTriangles[i] = createdTriangles[i];
    }

    createdTriangles.clear(); createdTriangles.shrink_to_fit();

    _numberVertices = totalNumberVertices;
    _numberTriangles = totalNumberTriangles;

    _vertices = newVertices;
    _triangles = newTriangles;
}

void IcoSphere::_ensureVertexOnUnitSphere(Vector3f& vertex) const
{
    const auto length = vertex.abs();
    vertex /= length;
}

}
