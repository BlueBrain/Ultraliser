#include "VolumeMesh.h"

namespace Ultraliser
{

VolumeMesh::VolumeMesh(const uint64_t numberVertices, const uint64_t numberTriangles)
{
    vertices.resize(numberVertices);
    triangles.resize(numberTriangles);
}

void VolumeMesh::append(const VolumeMesh* inputMesh)
{
    // Vertex offset is the number of vertices in the current mesh
    uint64_t vertexCountOffset = vertices.size();

    // Number of triangles in the current mesh
    uint64_t triangleCountOffset = triangles.size();

    // Expand and append the vertices and the triangles of the new mesh
    vertices.insert(vertices.end(), inputMesh->vertices.begin(), inputMesh->vertices.end());
    triangles.insert(triangles.end(), inputMesh->triangles.begin(), inputMesh->triangles.end());

    // Offset the new vertices to account for the addivity
    for (uint64_t i = triangleCountOffset; i < triangles.size(); ++i)
    {
        for(uint64_t j = 0; j < 3; ++j)
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
