#pragma once

#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{

/**
 * @brief The IcoSphere class
 */
class IcoSphere : public Mesh
{
public:

    /**
     * @brief IcoSphere
     * @param subdivisions
     */
    IcoSphere(const size_t& subdivisions);
    ~IcoSphere() { };

    /**
     * @brief getNumberSubdivisions
     * @return
     */
    size_t getNumberSubdivisions() const { return _subdivisions; };

private:

    /**
     * @brief _subdivideTriangleAtMidPointsAndNormalize
     * @param triangleIndex
     * @param vertexList
     * @param triangleList
     */
    void _subdivideTriangleAtMidPointsAndNormalize(const size_t& triangleIndex,
                                                   std::vector< Vector3f >& vertexList,
                                                   std::vector< Triangle >& triangleList);

    /**
     * @brief _subdivideTrianglesAtMidPointsAndNormalize
     */
    void _subdivideTrianglesAtMidPointsAndNormalize();

    /**
     * @brief _ensureVertexOnUnitSphere
     * @param vertex
     */
    void _ensureVertexOnUnitSphere(Vector3f& vertex) const;


public:

    /**
     * @brief _subdivisions
     */
    size_t _subdivisions;
};

}
