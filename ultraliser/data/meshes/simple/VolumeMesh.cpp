#include "VolumeMesh.h"

namespace Ultraliser
{

VolumeMesh::VolumeMesh(const size_t numberVertices, const size_t numberTriangles)
{
    vertices.resize(numberVertices);
    triangles.resize(numberTriangles);
}

void VolumeMesh::append(const VolumeMesh* inputMesh)
{
    // Vertex offset is the number of vertices in the current mesh
    size_t vertexCountOffset = vertices.size();

    // Number of triangles in the current mesh
    size_t triangleCountOffset = triangles.size();

    // Expand and append the vertices and the triangles of the new mesh
    vertices.insert(vertices.end(), inputMesh->vertices.begin(), inputMesh->vertices.end());
    triangles.insert(triangles.end(), inputMesh->triangles.begin(), inputMesh->triangles.end());

    // Offset the new vertices to account for the addivity
    for (size_t i = triangleCountOffset; i < triangles.size(); ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            triangles[i][j] += vertexCountOffset;
        }
    }
}

VolumeMesh* VolumeMesh::constructUnitCube(const float scale)
{
    // The cube will have 8 vertices, 6 sides each with two triangles, therefore 12 triangles
    VolumeMesh* unitCube = new VolumeMesh(8, 12);

    // Build the vertices
    unitCube->vertices[0] = Vertex(0, 0, 0) * scale;
    unitCube->vertices[1] = Vertex(1, 0, 0) * scale;
    unitCube->vertices[2] = Vertex(1, 1, 0) * scale;
    unitCube->vertices[3] = Vertex(0, 1, 0) * scale;
    unitCube->vertices[4] = Vertex(0, 0, 1) * scale;
    unitCube->vertices[5] = Vertex(1, 0, 1) * scale;
    unitCube->vertices[6] = Vertex(1, 1, 1) * scale;
    unitCube->vertices[7] = Vertex(0, 1, 1) * scale;

    // Build the triangles
    unitCube->triangles[0] = Triangle(0, 3, 1);
    unitCube->triangles[1] = Triangle(1, 3, 2);
    unitCube->triangles[2] = Triangle(5, 4, 0);
    unitCube->triangles[3] = Triangle(5, 0, 1);
    unitCube->triangles[4] = Triangle(6, 5, 1);
    unitCube->triangles[5] = Triangle(1, 2, 6);
    unitCube->triangles[6] = Triangle(3, 6, 2);
    unitCube->triangles[7] = Triangle(3, 7, 6);
    unitCube->triangles[8] = Triangle(4, 3, 0);
    unitCube->triangles[9] = Triangle(4, 7, 3);
    unitCube->triangles[10] = Triangle(7, 4, 5);
    unitCube->triangles[11] = Triangle(7, 5, 6);

    return unitCube;
}

VolumeMesh* VolumeMesh::constructVoxelCube(const Vector3f &pMin, const Vector3f &pMax)
{
    // Get the diagonal of the voxel
    Vector3f diagonal = pMax - pMin;

    // Transform for transforming the cube to the voxel location
    Matrix3f transform = Matrix3f::identity();

    for (int32_t i = 0; i < DIMENSIONS; ++i)
        transform(i, i) = diagonal[i];

    auto voxelCube = constructUnitCube();
    for (size_t i = 0; i < voxelCube->vertices.size(); ++i)
        voxelCube->vertices[i] = pMin + transform * voxelCube->vertices[i];

    return voxelCube;
}

void VolumeMesh::computeBoundingBox(Vector3f& pMinIn, Vector3f& pMaxIn)
{
    // Create new variables to avoid any mess if the inputs are already
    // initialized with some values
    Vector3f pMin(std::numeric_limits<float>::max());
    Vector3f pMax(std::numeric_limits<float>::lowest());

    for (size_t i = 0; i < vertices.size(); ++i)
    {
        Vertex vertex = vertices[i];
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

void VolumeMesh::centerAtOrigin(void)
{
    // Compute the bounding box
    Vector3f pMin, pMax;
    computeBoundingBox(pMin, pMax);

    // Compute the center
    Vector3f boundingBoxSize = pMax - pMin;
    Vector3f center = pMin + (0.5 * boundingBoxSize);

    // Shift the vertices to the origin to center the mesh
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] -= center;
    }
}

void VolumeMesh::rotate(const Matrix4f& matrix)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        Vector4f result = matrix * Vector4f(vertices[i]);
        vertices[i] = Vector3f(result.x(), result.y(), result.z());
    }
}

void VolumeMesh::transform(const Matrix4f& matrix)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        Vector4f result = matrix * Vector4f(vertices[i], 1.0);
        vertices[i] = Vector3f(result.x(), result.y(), result.z());
    }
}

void VolumeMesh::translate(const Vector3f& to)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] += to;
    }
}

void VolumeMesh::uniformScale(const float factor)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i] *= factor;
    }
}

void VolumeMesh::scale(const float x, const float y, const float z)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i].x() *= x;
        vertices[i].y() *= y;
        vertices[i].z() *= z;
    }
}

void VolumeMesh::scale(const Vector3f& factor)
{
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        vertices[i].x() *= factor.x();
        vertices[i].y() *= factor.y();
        vertices[i].z() *= factor.z();
    }
}

VolumeMesh::~VolumeMesh()
{
    // Cleaning vertices
    vertices.clear();
    vertices.shrink_to_fit();

    // Cleaning triangles
    triangles.clear();
    triangles.shrink_to_fit();
}

}
