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

#pragma once

#include <common/Common.h>
#include <math/Math.h>
#include <data/meshes/simple/NeighbourTriangle.h>
#include <data/meshes/simple/primitives/Primitives.h>
#include <data/morphologies/NeuronMorphology.h>
#include <data/morphologies/AstrocyteMorphology.h>

#include <data/morphologies/Section.h>
#include <data/common/ROI.h>
#include <data/structures/Neighbors.h>
#include <math/Math.h>

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
    Mesh(const size_t &numVertices, const size_t &numTriangles);

    /**
     * @brief OptimizationMesh
     * @param fileName
     */
    Mesh(const std::string &fileName, const bool &verbose = false);

    /**
     * @brief Mesh
     * @param samples
     * @param bevelSides
     */
    Mesh(const Samples& samples, const size_t &bevelSides = 16);

    /**
     * @brief Mesh
     * @param morphology
     */
    Mesh(const NeuronMorphology* morphology);

    /**
     * @brief Mesh
     * @param morphology
     */
    Mesh(const AstrocyteMorphology* morphology);

    /**
     * @brief append
     * Appends a new mesh to this mesh (or joint them)
     * @param inputMesh
     * The input mesh to be added to this mesh.
     */
    void append(const Mesh* inputMesh);

    /**
     * @brief import
     * @param filename
     * @param verbose
     */
    void import(const std::string &fileName, const bool& verbose = true);

    /**
     * @brief instanciate
     * Create a new instance of the mesh, but with a different number of
     * vertices and triangles.
     * @param numVertices
     * @param numTriangles
     * @return
     */
    Mesh* instanciate(const size_t &numVertices, const size_t &numTriangles);

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
     * @brief getTriangleBoundingBox
     * @param pMin
     * @param pMax
     */
    void getTriangleBoundingBox(const size_t &triangleIndex, Vector3f& pMin, Vector3f& pMax) const;

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
     * @brief getPositionSurfaceOnly
     * @param x
     * @param y
     * @param z
     * @param a
     * @param b
     * @param c
     * @return
     */
    Vector3f getPositionSurfaceOnly(const float& x, const float& y, const float& z,
                                    const int64_t &a, const int64_t &b, const int64_t &c);
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
    float computeDotProduct(const int64_t &a, const int64_t &b, const int64_t &c);

    /**
     * @brief computeCrossProduct
     * @param a
     * @param b
     * @param c
     * @return
     */
    Vector3f computeCrossProduct(const int64_t &a, const int64_t &b, const int64_t &c);

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
    bool checkFlipAction(const int64_t &a, const int64_t &b, const int64_t &c, const int64_t &d);

    /**
     * @brief computeEigenVector
     * @param index0
     * @param eigenValue
     * @param maxAngle
     * @return
     */
    EigenVector computeEigenVector(const int64_t &index0, Vector3f *eigenValue, float* maxAngle);

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
    float getAngleSurfaceOnly(const int64_t &a, const int64_t &b, const int64_t &c,
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


    void createEdgeList();
    /**
     * @brief subdivideTriangleAtCentroid
     * @param triangleIndex
     * @param vertexList
     * @param triangleList
     */
    void subdivideTriangleAtCentroid(const size_t& triangleIndex,
                                     std::vector< Vector3f >& vertexList,
                                     std::vector< Triangle >& triangleList);

    /**
     * @brief subdivideTriangleAtEdges
     * @param triangleIndex
     * @param vertexList
     * @param triangleList
     */
    void subdivideTriangleAtEdges(const size_t& triangleIndex,
                                  std::vector< Vector3f >& vertexList,
                                  std::vector< Triangle >& triangleList);

    /**
     * @brief refineSelectedTriangles
     * @param trianglesIndices
     */
    void refineSelectedTriangles(const std::vector< size_t > &trianglesIndices);

    /**
     * @brief refineROIs
     * @param regions
     */
    void refineROIs(const ROIs& regions);

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
    void optimize(const size_t &optimizationIterations,
                  const int64_t &smoothingIterations,
                  const float& denseFactor);

    /**
     * @brief optimizeAdaptively
     * @param optimizationIterations
     * @param smoothingIterations
     * @param flatFactor
     * @param denseFactor
     */
    void optimizeAdaptively(const size_t &optimizationIterations,
                            const size_t &smoothingIterations,
                            const float &flatFactor,
                            const float &denseFactor);

    /**
     * @brief optimizeAdapttivelyWithROI
     * @param optimizationIterations
     * @param smoothingIterations
     * @param flatCoarseFactor
     * @param denseCoarseFactor
     * @param regions
     */
    void optimizeAdapttivelyWithROI(const size_t &optimizationIterations,
                                    const size_t &smoothingIterations,
                                    const float &flatFactor,
                                    const float &denseFactor,
                                    const ROIs &regions);

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
                    const bool &formatOFF = false,
                    const bool &formatSTL = false) const;

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
                        const std::string *prefix = nullptr) const;

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
    size_t getNumberVertices() const;

    /**
     * @brief getNumberTriangles
     * @return
     */
    size_t getNumberTriangles() const;

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
     * @brief shrinkSurface
     * @param numIterations
     * @param lambda
     * @param mu
     */
    void shrinkOrInflateSurface(const Neighborhood &vertexNeighbours,
                                const Neighborhood faceNeighbors,
                                const bool &shrink,
                                const float& lambda = 0.5f,
                                const float& mu = -0.4);

    /**
     * @brief smoothSurface
     * Another simple implementation to the Laplacian filter.
     * @param numIterations
     * Number of iterations.
     */
    void smoothSurface(size_t numIterations);

    /**
     * @brief Mesh simplification using edge collapse. One vertex per iteration is collapsed,
     * until either the requested percentage is satisfied or no further edges can be collapsed.
     * @param vertexPercentage Maximum percentage of vertices to remove.
     */
    void collapseEdges(float vertexPercentage);

    /**
     * @brief scaleAndTranslateGeneratedMesh
     * Scale the translate the generated mesh to fit the dimensions of the input mesh.
     * @param center
     * The center of the input mesh.
     * @param BoundingBox
     * The bounding box of the input mesh.
     */
    void scaleAndTranslate(const Vector3f &center, const Vector3f &BoundingBox);

    /**
     * @brief updateData
     * @param vertices
     * @param numberVertices
     * @param triangles
     * @param numberTriangles
     */
    void updateData(const Vertices vertices, const Triangles triangles);

    /**
     * @brief relaseData
     * Releases data temporarily.
     * NOTE: Use this function carefully. If you release the data, make sure not to use any of
     * the public functions of the Mesh untill the updateData() function is called again.
     */
    void relaseData() { _releaseData(); }

private:

    /**
     * @brief INitializes the mesh from a list of vertices and a list of triangles
     * @param vertices Mesh vertices
     * @param triangles Mesh triangles
     */
    void _initFromVertexAndTriangleList(Vertices vertices, Triangles triangles);

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

    /**
     * @brief _computeNeighborhoods
     * @param mesh
     * @param vertexNeighbors
     * @param faceNeighbors
     */
    void _computeNeighborhoods(Mesh& mesh,
                               Neighborhood& vertexNeighbors,
                               Neighborhood& faceNeighbors);

    /**
     * @brief _computeKernel
     * @param mesh
     * @param vertexIndex
     * @param verticesN
     * @param facesN
     * @return
     */
    Vector3f _computeKernel(Mesh& mesh,
                            const size_t vertexIndex,
                            Neighbors& verticesN,
                            Neighbors& facesN);

    /**
     * @brief _computeCotangentWeight
     * @param mesh
     * @param vertexIndex
     * @param neighborIndex
     * @param faceN
     * @return
     */
    float _computeCotangentWeight(Mesh& mesh,
                                  const uint32_t vertexIndex,
                                  const uint32_t neighborIndex,
                                  Neighbors& faceN);

    /**
     * @brief _smoothVertex
     * @param mesh
     * @param vertexIndex
     * @param kernel
     * @param param
     * @return
     */
    Vector3f _smoothVertex(Mesh& mesh,
                           const size_t vertexIndex,
                           const Vector3f kernel,
                           const float param);

    /**
     * @brief _triangleContainsVertex
     * @param t
     * @param vIndex
     * @return
     */
    inline bool _triangleContainsVertex(const Triangle & t, const size_t vIndex)
    {
        return (I2UI64(t.x()) == vIndex || I2UI64(t.y()) == vIndex || I2UI64(t.z()) == vIndex);
    }

    /**
     * @brief destroyNeighborlist
     */
    void _destroyNeighborlist();

    /**
     * @brief _createNeighbourList
     */
    void _createNeighbourList();

    /**
     * @brief _updateVertexMarkers
     */
    void _resetVertexMarkers();

    /**
     * @brief _destroyVertexMarkers
     */
    void _destroyVertexMarkers();

    /**
     * @brief _updateTriangleMarkers
     */
    void _resetTriangleMarkers();

    /**
     * @brief _destroyTriangleMarkers
     */
    void _destroyTriangleMarkers();

    /**
     * @brief _selectVerticesInROI
     * @param regions
     */
    void _selectVerticesInROI(const ROIs& regions);

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
     * A Pointer to the neighbors (triangles).
     */
    Triangle* _neighbors;

private:

    /**
     * @brief neighborList
     * Pointer to neighbour list.
     */
    NeighborTriangle** _neighborList;

    /**
     * @brief _vertexMarkers
     * The evrtex markers are used to preserve specific vertices in the mesh.
     */
    std::vector< uint8_t > _vertexMarkers;

    /**
     * @brief numberVertices
     * Number of vertices in the mesh.
     */
    size_t _numberVertices;

    /**
     * @brief numberTriangles
     * Number of triangles (or faces) in the meshes.
     */
    size_t _numberTriangles;

    /**
     * @brief defaultOptimizationTime
     */
    double _optimizationTime;

    // Mesh Statistics
    friend class MeshStatistics;
};

}

