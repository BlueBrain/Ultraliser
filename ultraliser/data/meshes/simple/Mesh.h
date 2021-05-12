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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_MESH_H
#define ULTRALISER_DATA_MESH_SIMPLE_MESH_H

#include <common/Common.h>
#include <math/Math.h>
#include <data/meshes/simple/primitives/Primitives.h>
#include <data/meshes/simple/NeighbourTriangle.h>
#include <data/structures/Neighbors.h>
#include <data/morphologies/Section.h>

namespace Ultraliser
{

/**
 * @brief The Mesh class
 */
class Mesh
{

public:
    /**
     * @brief OptimizationMesh
     */
    Mesh();
    ~Mesh();

    /**
     * @brief OptimizationMesh
     * @param vertices
     * @param triangles
     */
    Mesh(Vertices _vertices, Triangles _triangles);

    /**
     * @brief OptimizationMesh
     * @param numVertices
     * @param numTriangles
     */
    Mesh(const uint64_t &numVertices, const uint64_t &numTriangles);

    /**
     * @brief OptimizationMesh
     * @param fileName
     */
    Mesh(const std::string &fileName);

    /**
     * @brief Mesh
     * @param section
     */
    Mesh(const Samples& samples);

    /**
     * @brief import
     * @param filename
     * @param verbose
     */
    void import(const std::string &fileName);

    /**
     * @brief instanciate
     * Create a new instance of the mesh, but with a different number of
     * vertices and triangles.
     * @param numVertices
     * @param numTriangles
     * @return
     */
    Mesh* instanciate(const uint64_t &numVertices,
                      const uint64_t &numTriangles);

    /**
     * @brief getVertices
     * @return
     */
    const Vertex* getVertices() const;

    /**
     * @brief getTriangles
     * @return
     */
    const Triangle* getTriangles() const;

    /**
     * @brief translate
     * @param to
     */
    void translate(const Vector3f& to);

    /**
     * @brief rotate
     * @param matrix
     */
    void rotate(const Matrix4f &matrix);

    /**
     * @brief transform
     * @param matrix
     */
    void transform(const Matrix4f &matrix);

    /**
     * @brief uniformScale
     * @param factor
     */
    void uniformScale(const float factor);

    /**
     * @brief scale
     * @param x
     * @param y
     * @param z
     */
    void scale(const float x = 1.0, const float y = 1.0, const float z = 1.0);

    /**
     * @brief computeBoundingBox
     * @param pMin
     * @param pMax
     */
    void computeBoundingBox(Vector3f& pMin, Vector3f& pMax);

    /**
     * @brief computeRelaxedBoundingBox
     * @param pMinIn
     * @param pMaxIn
     * @param relaxationPercentage
     */
    void computeRelaxedBoundingBox(Vector3f& pMinIn,
                                   Vector3f& pMaxIn,
                                   const float& relaxationPercentage = 1.0);

    Mesh & operator = (const Mesh mesh);

    /**
     * @brief centerAtOrigin
     */
    void centerAtOrigin(void);

    /**
     * @brief destroyNeighborlist
     */
    void destroyNeighborlist();

    /**
     * @brief createNeighbourList
     */
    void createNeighbourList();

    /**
     * @brief getPositionSurfaceOnly
     * @param x
     * @param y
     * @param z
     * @param a
     * @param b
     * @param c
     * @return
     */
    Vector3f getPositionSurfaceOnly(const float& x,
                                    const float& y,
                                    const float& z,
                                    const int64_t &a,
                                    const int64_t &b,
                                    const int64_t &c);
    /**
     * @brief computeDotProduct
     * Computes the dot product of a triangle vertices
     * @param a
     * First vertex index.
     * @param b
     * Second vertex index.
     * @param c
     * Third vertex index.
     * @return
     */
    float computeDotProduct(const int64_t &a,
                            const int64_t &b,
                            const int64_t &c);

    /**
     * @brief computeCrossProduct
     * @param a
     * @param b
     * @param c
     * @return
     */
    Vector3f computeCrossProduct(const int64_t &a,
                                 const int64_t &b,
                                 const int64_t &c);

    /**
     * @brief rotate
     * @param sx
     * @param sy
     * @param sz
     * @param theta
     * @param phi
     * @param angle
     * @return
     */
    Vector3f rotate(const float& sx, const float& sy, const float& sz,
                    const float& theta, const float& phi, const float& angle);

    /**
     * @brief edgeFlipping
     * @param index
     */
    void edgeFlipping(int64_t index);

    /**
     * @brief checkFlipAction
     * @param a
     * @param b
     * @param c
     * @param d
     * @return
     */
    bool checkFlipAction(const int64_t &a, const int64_t &b,
                         const int64_t &c, const int64_t &d);

    /**
     * @brief computeEigenVector
     * @param index0
     * @param eigenValue
     * @param maxAngle
     * @return
     */
    EigenVector computeEigenVector(const int64_t &index0,
                                   Vector3f *eigenValue,
                                   float* maxAngle);

    /**
     * @brief moveVertexAlongSurface
     * Moves a vertex along the surface.
     * @param index The index of the vertex in the mesh.
     */
    void moveVertexAlongSurface(int64_t index);

    /**
     * @brief isValidVertex
     * @param v
     * @return
     */
    bool isValidVertex(const Vector3f& v);

    /**
     * @brief getAngleSurfaceOnly
     * @param a
     * @param b
     * @param c
     * @return
     */
    float getAngleSurfaceOnly(const int64_t &a,
                              const int64_t &b,
                              const int64_t &c,
                              bool &angleError);

    /**
     * @brief getVertexNormal
     * Computes the normal of a vertex.
     * @param index
     * Vertex index.
     * @return
     * The computed normal.
     */
    Vector3f getVertexNormal(const int64_t &index);

    /**
     * @brief computeAngles
     * @param computedMinAngle
     * @param computedMaxAngle
     * @param numSmall
     * @param numLarge
     * @param maxMinAngle
     * @param minMaxAngle
     */
    void computeAngles(float* computedMinAngle,
                       float* computedMaxAngle,
                       int64_t* numSmall, int64_t* numLarge,
                       int64_t maxMinAngle, int64_t minMaxAngle);

    /**
     * @brief subdividePolygin
     * @param starNGR
     * @param triangleAvailableList
     * @param triangleAvailableIndex
     */
    void subdividePolygin(NeighborTriangle *starNGR,
                          int64_t *triangleAvailableList,
                          int64_t *triangleAvailableIndex);

    /**
     * @brief smoothNormal
     * @param index
     */
    void smoothNormal(const int64_t index);

    /**
     * @brief smooth
     * Smoothes the mesh.
     * @param maxMinAngle
     * @param minMaxAngle
     * @param maxIterations
     * @param flipEdges
     * @return
     */
    bool smooth(const int64_t &maxMinAngle = 15,
                const int64_t &minMaxAngle = 150,
                const int64_t &maxIterations = 6,
                const bool& flipEdges = true);

    /**
     * @brief smoothNormals
     * @return
     */
    void smoothNormals();

    /**
     * @brief refine
     * Refine the durface mesh.
     */
    void refine();

    /**
     * @brief coarse
     * @param coarseRate
     * @param flatnessRate
     * @param densenessWeight
     * @param maxNormalAngle
     * @param iteration
     * @return
     */
    bool coarse(const float& coarseRate,
                const float& flatnessRate,
                const float& densenessWeight,
                const float& maxNormalAngle,
                const int64_t &iteration = 1);

    /**
     * @brief coarseDense
     * @param denseRate
     * @param iterations
     */
    void coarseDense(const float &denseRate, const int64_t &iterations = 1);

    /**
     * @brief coarseFlat
     * @param flatnessRate
     * @param iterations
     */
    void coarseFlat(const float &flatnessRate, const int64_t &iterations = 1);

    /**
     * @brief optimizeUsingDefaultParameters
     */
    void optimizeUsingDefaultParameters();

    /**
     * @brief optimize
     * @param smoothingIterations
     * @param denseFactor
     */
    void optimize(const uint64_t &optimizationIterations,
                  const int64_t &smoothingIterations,
                  const float& denseFactor);

    /**
     * @brief optimizeAdaptively
     * @param optimizationIterations
     * @param smoothingIterations
     * @param flatFactor
     * @param denseFactor
     */
    void optimizeAdaptively(const uint64_t &optimizationIterations,
                            const uint64_t &smoothingIterations,
                            const float &flatFactor,
                            const float &denseFactor);

    /**
     * @brief exportMesh
     * @param prefix
     * @param formatOBJ
     * @param formatPLY
     * @param formatOFF
     * @param formatSTL
     */
    void exportMesh(const std::string &prefix,
                    const bool &formatOBJ = false,
                    const bool &formatPLY = false,
                    const bool & formatOFF = false,
                    const bool & formatSTL = false) const;

    /**
     * @brief writeDistributions
     * @param reference
     * @param prefix
     */
    void writeDistributions(const std::string &reference,
                            const std::string *prefix) const;

    /**
     * @brief printStats
     * @param reference
     * @param prefix
     */
    void printStats(const std::string &reference,
                        const std::string* prefix = nullptr) const;

    /**
     * @brief getDefaultOptimizationTime
     * @return
     */
    double getDefaultOptimizationTime() const;

    /**
     * @brief removeFloatingFaces
     */
    void removeFloatingFaces();

    /**
     * @brief getNumberVertices
     * @return
     */
    uint64_t getNumberVertices() const;

    /**
     * @brief getNumberTriangles
     * @return
     */
    uint64_t getNumberTriangles() const;

    /**
     * @brief applyLaplacianSmooth
     * Smooth the mesh surface using Taubian smoothing (Uniform Laplacian
     * Kernel with cotangent weight + inflation).
     * @param numIterations
     * The number of smoothing iterations to apply.
     * @param smoothLambda
     * Smoothing kernel multiplier (must be equal or greater than 0).
     * @param inflateMu
     * Inflate kernel multiplier (must be equal or less than 0).
     */
    void applyLaplacianSmooth(const uint32_t &numIterations = 1,
                              const float& smoothLambda = 0.2,
                              const float& inflateMu = 0.1);

    /**
     * @brief scaleAndTranslateGeneratedMesh
     * Scale the translate the generated mesh to fit the dimensions of the input mesh.
     * @param center
     * The center of the input mesh.
     * @param BoundingBox
     * The bounding box of the input mesh.
     */
    void scaleAndTranslate(const Vector3f &center, const Vector3f &BoundingBox);

private:

    /**
     * @brief _releaseData
     * Releases the data of the mesh.
     */
    void _releaseData();

    /**
     * @brief _cotangentAngle
     * @param pivot
     * @param a
     * @param b
     * @return
     */
    float _cotangentAngle(const Vector3f& pivot,
                           const Vector3f& a,
                           const Vector3f& b);

    void _computeNeighborhoods(Mesh& mesh,
                               Neighborhood& vertexNeighbors,
                               Neighborhood& faceNeighbors);

    Vector3f _computeKernel(Mesh& mesh,
                            const uint64_t vertexIndex,
                            Neighbors& verticesN,
                            Neighbors& facesN);

    float _computeCotangentWeight(Mesh& mesh,
                                  const uint32_t vertexIndex,
                                  const uint32_t neighborIndex,
                                  Neighbors& faceN);

    Vector3f _smoothVertex(Mesh& mesh,
                           const uint64_t vertexIndex,
                           const Vector3f kernel,
                           const float param);

    inline bool _triangleContainsVertex(const Triangle & t, const uint64_t vIndex)
    {
        return (I2UI64(t.x()) == vIndex ||
                I2UI64(t.y()) == vIndex ||
                I2UI64(t.z()) == vIndex);
    }

public:

    /**
     * @brief vertices
     * A list of all the vertices in the mesh.
     */
    Vector3f* _vertices;

    /**
     * @brief triangles
     * A list of all the triangles in the mesh.
     */
    Triangle* _triangles;

    /**
     * @brief neighbors
     * A point32_ter to the neighbors (triangles).
     */
    Triangle* _neighbors;

private:

    /**
     * @brief neighborList
     * Point32_ter to neighbour list.
     */
    NeighborTriangle** _neighborList;

    /**
     * @brief numberVertices
     * Number of vertices in the mesh.
     */
    uint64_t _numberVertices;

    /**
     * @brief numberTriangles
     * Number of triangles (or faces) in the meshes.
     */
    uint64_t _numberTriangles;

    /**
     * @brief defaultOptimizationTime
     */
    double _optimizationTime;

    // Mesh Statistics
    friend class MeshStatistics;
};

}

#endif // ULTRALISER_DATA_MESH_SIMPLE_MESH_H
