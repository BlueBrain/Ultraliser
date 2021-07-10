/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marco Attene < IMATI-GE / CNR >
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
 * The content of this file is based on MeshFix. The code has been modified under the terms of
 * the GNU General Public License as published by the Free Software Foundation either version 3 of
 * the License, or (at your option) any later version.
 * MeshFix has a dual license for free and commercial use. For further information, please refer
 * to the original repository at < https://github.com/MarcoAttene/MeshFix-V2.1>.
 **************************************************************************************************/

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_MESH_H
#define ULTRALISER_DATA_MESHES_ADVANCED_MESH_H

#include <data/meshes/advanced/math/Matrix.h>
#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/primitives/AdvancedEdge.h>
#include <data/meshes/advanced/primitives/AdvancedTriangle.h>
#include <data/meshes/advanced/primitives/AdvancedIndexedMesh.hh>
#include <data/meshes/simple/primitives/Primitives.h>
#include <data/structures/Neighbors.h>
#include <data/meshes/simple/Mesh.h>

namespace Ultraliser
{

inline void swapPointers(void **a, void **b)
{
    void *t = *a;
    *a = *b;
    *b = t;
}

/**
 * @brief The AdvancedMesh class
 * This class represents a advanced and oriented triangle mesh.
 * Vertices, Edges and Triangles are stored in the Lists V, E and T
 * respectively. Methods boundaries(), handles() and shells() may
 * be used to retrieve the respective topological entities.
 * Navigation of the mesh is based on the topological relationships
 * stored in each Vertex, Edge and Triangle.
 * Some methods would require a global update only to maintain
 * consistent values of the protected fields n_boundaries, n_handles
 * and n_shells. On the other hand, the same methods would work in
 * constant time if these values would not need to be updated.
 * To keep complexity as low as possible,
 * we make use of the 'dirty bits' d_boundaries, d_handles and
 * d_shells to mark that the respective entities must be updated.
 * The access functions boundaries(), handles() and shells()
 * check the status of the respective dirty bit and do a global
 * update (i.e., eulerUpdate()) only if necessary.
 * The complexity of the methods is provided as a function of a
 * generic 'N', which is O(V.numels()) = O(E.numels()) = O(T.numels()).
 */
class AdvancedMesh
{

public:

    /**
     * @brief _vertices
     * Vertex set
     */
    List _vertices;

    /**
     * @brief _edges
     * Edge set
     */
    List _edges;

    /**
     * @brief _triangles
     * Triangle set
     */
    List _triangles;

protected:

    /**
     * @brief _loadedNumberVertices
     */
    uint64_t _loadedNumberVertices;

    /**
     * @brief _loadedNumberTriagnles
     */
    uint64_t _loadedNumberTriagnles;

    /**
     * @brief _nBoundaries
     * Number of boundary loops
     */
    int _nBoundaries;

    /**
     * @brief _nHandles
     * Number of handles
     */
    int _nHandles;

    /**
     * @brief _nShells
     * Number of connected components
     */
    int _nShells;

    /**
     * @brief _dBoundaries
     * Dirty bit for n_boundaries
     */
    bool _dBoundaries;

    /**
     * @brief _dHandles
     * Dirty bit for n_handles
     */
    bool _dHandles;

    /**
     * @brief _dShells
     * Dirty bit for n_shells
     */
    bool _dShells;

public:

    /**
     * @brief info
     * Generic information attached to this mesh
     */
    void *info;

public:

    /**
     * @brief Mesh
     * Empty triangulation. Should be used only prior to a call to load().
     */
    AdvancedMesh();

    /**
     * @brief Mesh
     * Loads a file directly from a file.
     * @param filePath
     */
    AdvancedMesh(const std::string filePath);

    /**
     * @brief Mesh
     * @param vertices
     * @param triangles
     */
    AdvancedMesh(Vertices vertices, Triangles triangles);

    /**
     * @brief Mesh
     * @param vertices
     * @param numberVertices
     * @param triangles
     * @param numberTriangles
     */
    AdvancedMesh(const Vertex *vertices, const uint64_t &numberVertices,
                 const Triangle *triangles, const uint64_t &numberTriangles);

    /**
     * @brief Mesh
     * Pre-defined triangulation. Currently, only "triangle" and "tetrahedron"
     * are recognized.
     */
    AdvancedMesh(const char *inputMeshDefinition);

    /**
     * @brief init
     * @param inputMeshDefinition
     */
    void init(const char *inputMeshDefinition);

    /**
     * @brief Mesh
     *  Clones an existing Trianglation.
     * @param clone_info
     */
    AdvancedMesh(const AdvancedMesh *input, const bool cloneInfo = false);

    /**
     * @brief init
     * @param cloneInfo
     */
    void init(const AdvancedMesh *input, const bool cloneInfo = false);

    /**
     * @brief Mesh
     *  Clones an existing connected component.

     * Creates a new Mesh out of a connected component of an existing
     * Mesh. 't' is a triangle of the connected component that must
     * be copied. If 'keepReferences' is TRUE, each element of the existing mesh
     * keeps a pointer to the corresponding new element in the 'info' field.
     * @param t
     * @param keepReferences
     */
    AdvancedMesh(const AdvancedTriangle *t, const bool keepReferences = false);

    /**
     * @brief init
     * @param t
     * @param keepReferences
     */
    void init(const AdvancedTriangle *t, const bool keepReferences = false);

    /**
     * Destructor. Frees the memory allocated for all the mesh elements.
     * Warning! This method uses the freeNodes() method of the class List,
     * which is not guaranteed to work correctly on systems other than
     * Linux (see the documentation of List for details).
     * Assuming that the method works correctly, however, the calling
     * function is responsible of freeing the memory that was possibly
     * allocated for objects pointed to by the 'info' field of mesh
     * elements. Clearly, this must be done before calling the destructor.
     */
    ~AdvancedMesh();

    /**
     * @brief isBaseType
     * Returns true only if object is a basic Mesh.
     * All the reimplementations must return false.
     * @return
     */
    //  bool isBaseType() const { return true; }

    /**
     * @brief boundaries
     * Get the number of boundary loops of the triangle mesh. O(1) or O(N).
     * @return
     */
    int boundaries()
    {
        if (_dBoundaries)
        {
            eulerUpdate();
        }

        return _nBoundaries;
    }

    /**
     * @brief handles
     * Get the number of handles of the triangle mesh. O(1) or O(N).
     * @return
     */
    int handles()
    {
        if (_dHandles)
        {
            eulerUpdate();
        }

        return _nHandles;
    }

    /**
     * @brief shells
     * Get the number of connected components of the triangle mesh. O(1) or O(N)
     * @return
     */
    int shells()
    {
        if (_dShells)
        {
            eulerUpdate();
        }

        return _nShells;
    }

    /**
     * @brief load
     * Initialize the triangle mesh from the file 'filename'.

     * The file format is automatically deduced from the magic number
     * or the filename extension. If 'update' is FALSE, the
     * global update for the topological entities is prevented.
     * Currently, the following file formats are supported:
     * Open Inventor (IV), VRML 1.0 and 2.0 (WRL), Object File Format (OFF),
     * IMATI Ver-Tri (VER, TRI), PLY, OBJ, STL.
     * A non-zero value is returned in case of error. Specifically,
     * IO_CANTOPEN means that the file couldn't be opened for reading.
     * IO_FORMAT means that the file format was not recognized by the loader.
     * IO_UNKNOWN represents all the other errors.
     * The calling function is responsible of verifying that the mesh is
     * empty before calling this method.
     * @param fileName
     * @param update
     * @return
     */
    void importMesh(const std::string &filePath, const bool update = 0);

    /**
     * @brief loadOFF
     * Loads OFF
     * @return
     */
    void importOFF(const std::string &filePath);

    /**
     * @brief loadPLY
     * Loads PLY
     * @return
     */
    void importPLY(const std::string &fileName);

    /**
     * @brief loadOBJ
     * @param fileName
     * @return
     */
    void importOBJ(const std::string &fileName);

    /**
     * @brief AdvancedMesh::importSTL
     *
     * ASCII STL
     * An ASCII STL file begins with the line
     *      solid name
     * where name is optional string (though if name is omitted there must still be space after solid).
     * The file continues with any number of triangles, each represented as follows:
     *
     * facet normal ni nj nk
     *      outer loop
     *          vertex v1x v1y v1z
     *          vertex v2x v2y v2z
     *          vertex v3x v3y v3z
     *      endloop
     *  endfacet
     * where each n or v is a floating-point number in sign-mantissa-"e"-sign-exponent format.
     * The file concludes with
     *  endsolid name
     *
     * Binary STL
     * UINT8[80]    – Header                 -     80 bytes
     * UINT32       – Number of triangles    -      4 bytes
     * foreach triangle                      -  50 bytes:
     *      REAL32[3] – Normal vector             - 12 bytes
     *      REAL32[3] – Vertex 1                  - 12 bytes
     *      REAL32[3] – Vertex 2                  - 12 bytes
     *      REAL32[3] – Vertex 3                  - 12 bytes
     *      UINT16    – Attribute byte count      -  2 bytes
     * end
     *
     * @param filePath
     * @return
     */
    void importSTL(const std::string &fileName);

    /**
     * @brief cutAndStitch
     * Convert to advanced.
     * @return
     */
    int cutAndStitch();

    /**
     * @brief createIndexedTriangle
     * @return
     */
    AdvancedTriangle * createIndexedTriangle(ExtendedVertex **vertexList, int64_t, int64_t, int64_t);

    /**
     * @brief CreateTriangleFromVertices
     * @return
     */
    AdvancedTriangle * createTriangleFromVertices(ExtendedVertex *vertex1,
                                          ExtendedVertex *vertex2,
                                          ExtendedVertex *vertex3);

    /**
     * @brief coordBackApproximation
     * This function approximates the vertex coordinates with the values
     * that can be represented in an ASCII file.
     */
    void coordBackApproximation();

    /**
     * @brief getNumberVertices
     * @return
     */
    uint64_t getNumberVertices() const
    {
        return _vertices.numberElements();
    }

    /**
     * @brief getNumberEdges
     * @return
     */
    uint64_t getNumberEdges() const
    {
        return _edges.numberElements();
    }

    /**
     * @brief getNumberTriangles
     * @return
     */
    uint64_t getNumberTriangles() const
    {
        return _triangles.numberElements();
    }

    /**
     * @brief getNumberSelfIntersectingFaces
     * @return
     */
    uint64_t getNumberSelfIntersectingFaces();

protected:

    /**
     * @brief _closeLoadingSession
     */
    void _closeLoadingSession(FILE *file, int loadedFaces, ExtendedVertex **, bool);

    /**
     * @brief _pinch
     * Implements the cutting and stitching procedure to convert to advanced mesh.
     * Assumes that singular edges to be cut and stitched are marked as BIT5.
     * @param inputEdge
     * @param withCommonVertex
     * @return
     */
    bool _pinch(AdvancedEdge *inputEdge, bool withCommonVertex);

    /**
     * @brief duplicateEdge
     * If the 'e' is internal, creates a copy of the edge and associates e->t2
     * to this copy. After this operation, 'e' and its copy are boundary edges
     * with the same vertices. If 'e' is already on boundary, nothing is done
     * and nullptr is returned. Otherwise the newly created copy is returned.
     * @param e
     * @return
     */
     AdvancedEdge *duplicateEdge(AdvancedEdge *e);

public:

    /**
     * @brief newVertex
     * @return
     */
     AdvancedVertex* newVertex();

    /**
     * @brief newVertex
     * @return
     */
     AdvancedVertex* newVertex(const double &, const double &, const double &);

    /**
     * @brief newVertex
     * @return
     */
     AdvancedVertex* newVertex(AdvancedPoint *);

    /**
     * @brief newVertex
     * @return
     */
     AdvancedVertex* newVertex(AdvancedPoint &);

    /**
     * @brief newVertex
     * @return
     */
     AdvancedVertex* newVertex(AdvancedVertex *);

    /**
     * @brief newEdge
     * @return
     */
     AdvancedEdge* newEdge(AdvancedVertex *, AdvancedVertex *);

    /**
     * @brief newEdge
     * @return
     */
     AdvancedEdge* newEdge(AdvancedEdge *);

    /**
     * @brief newTriangle
     * @return
     */
     AdvancedTriangle* newTriangle();

    /**
     * @brief newTriangle
     * @return
     */
     AdvancedTriangle* newTriangle(AdvancedEdge *, AdvancedEdge *, AdvancedEdge *);

    /**
     * @brief newObject
     * @return
     */
     AdvancedMesh* newObject() const { return new AdvancedMesh(); }

    /**
     * @brief newObject
     * @param tm
     * @param ci
     * @return
     */
     AdvancedMesh* newObject(const AdvancedMesh *tm, const bool ci = false) const
     {
         return new AdvancedMesh(tm, ci);
     }

    /**
     * @brief newObject
     * @param s
     * @return
     */
     AdvancedMesh* newObject(const char *s) const { return new AdvancedMesh(s); }

     /**
      * @brief toSimpleMesh
      * Convertes the advanced mesh to a simple mesh.
      * @return
      */
     Mesh* toSimpleMesh() const;

    /**
     * @brief exportMesh
     * @param prefix
     * @param exportOBJ
     * @param exportPLY
     * @param exportOFF
     * @param exportSTL
     */
    void exportMesh(const std::string &prefix,
                    const bool& formatOBJ = false,
                    const bool& formatPLY = false,
                    const bool& formatOFF = false,
                    const bool& formatSTL = false);
    /**
     * @brief saveOFF
     *  Saves OFF 1.0
     * @return
     */
    void exportOFF(const std::string &filePrefix);

    /**
     * @brief saveOBJ
     * Saves OBJ
     * @return
     */
    void exportOBJ(const std::string &fileName);

    /**
     * @brief saveSTL
     * Saves STL
     * @return
     */
    void exportSTL(const std::string &filePath);

    /**
     * @brief savePLY
     * Saves PLY 1.0 (ascii or binary)
     * @param ascii
     * @return
     */
    void exportPLY(const std::string &filePath, bool writeASCII = true);


    /**
     * @brief append
     * Append another triangle mesh to the existing one.

     * This method works exactly as the 'load()' method, except for the fact
     * that it does not assume that the mesh is empty.
     * @param filename
     * @param update
     * @return
     */
    void append(const std::string &filePath, const bool update = 1);

    /**
     * @brief moveMeshElements
     * Move all the elements of 't' to this mesh and delete 't' itself.
     * @param t
     * @param delInput
     */
     void moveMeshElements(AdvancedMesh *inputMesh, bool deleteInput = true);

    /**
     * @brief createEdge
     * Creates an Edge connecting two existing mesh vertices.

     * Returns the newly created edge. If an edge connecting the two vertices
     * already exists in the mesh, then no new edge is created and the old one
     * is returned.
     * @param v1
     * @param v2
     * @return
     */
    AdvancedEdge* createEdge(AdvancedVertex *v1, AdvancedVertex *v2);

    /**
     * @brief createEdge
     * Creates an Edge connecting two existing mesh Extended vertices.

     * Returns the newly created edge. If an edge connecting the two vertices
     * already exists in the mesh, then no new edge is created and the old one
     * is returned.
     * If 'check' is FALSE, the check for previously existing edges is skipped.
     * @param v1
     * @param v2
     * @param check
     * @return
     */
    AdvancedEdge* createEdge(ExtendedVertex *v1, ExtendedVertex *v2,
                             const bool check = 1);

    /**
     * @brief createTriangle
     *  Creates a properly oriented Triangle bounded by three existing mesh edges.

     * Returns the newly created Triangle. If e1, e2 and e3
     * are not suitable for creating a properly oriented and
     * advanced triangle, the creation fails and nullptr is returned.
     * @param e1
     * @param e2
     * @param e3
     * @return
     */
     AdvancedTriangle * createTriangle(AdvancedEdge *e1, AdvancedEdge *e2, AdvancedEdge *e3);

    /**
     * @brief createUnorientedTriangle
     * Creates an arbitrarily oriented Triangle bounded by three existing mesh edges.

     * Returns the newly created Triangle. If either e1, e2 or e3
     * has already two incident triangles, the creation fails and nullptr is returned.
     * This method assumes that e1, e2 and e3 are incident to exactly three vertices.
     * @return
     */
     AdvancedTriangle * createUnorientedTriangle(AdvancedEdge *, AdvancedEdge *, AdvancedEdge *);


    /**
     * @brief eulerEdgeTriangle
     *  Creates a newEdge 'e' and an oriented Triangle bounded by
     * 'e', 'e1' and 'e2'.

     * The newly created triangle is returned, unless 'e1' and 'e2' do
     * not share a vertex or they are not boundary edges.
     * In this cases, nullptr is returned.
     * @param e1
     * @param e2
     * @return
     */
     AdvancedTriangle *eulerEdgeTriangle(AdvancedEdge *e1, AdvancedEdge *e2);


    /**
     * @brief splitEdge
     * Splits and edge at a given point and returns the newly created Vertex.
     * If the boolean parameter is set to true, the 'mask' fields of edges and
     * triangles are propagated to the new elements.
     * @return
     */
     AdvancedVertex *splitEdge(AdvancedEdge *, AdvancedPoint *, bool copyMask = 0);

     /**
      * @brief getVerticesAndTrianglesArray
      * @param vertexArray
      * @param triangleArray
      * @param numberVertices
      * @param numberTriangles
      */
     void getVerticesAndTrianglesArray(Vertex *&vertexArray,
                                       Triangle *&triangleArray,
                                       uint64_t &numberVertices,
                                       uint64_t &numberTriangles) const;

     /**
      * @brief printStats
      * Prints the statistics of the mesh
      * @param reference
      * @param prefix
      */
     void printStats(const std::string &reference, const std::string *prefix= nullptr);

     /**
      * @brief writeDistributions
      * @param reference
      * @param prefix
      */
     void writeDistributions(const std::string &reference, const std::string *prefix) const;

    /**
     * @brief splitTriangle
     * Splits a triangle at a given point and returns the newly created Vertex.
     * If the boolean parameter is set to true, the 'mask' field of the
     * triangle is propagated to the new triangle.
     * @return
     */
     AdvancedVertex *splitTriangle(AdvancedTriangle *triangle, AdvancedPoint *point, bool copyMask = 0);

    /**
     * @brief bridgeBoundaries
     * Creates two new triangles connecting the boundary edges e1 and e2
     * and returns their common edge.
     * If e1 and e2 share a vertex, then only one triangle is created and
     * e1 is returned.
     * Returns nullptr if either e1 or e2 are not boundary edges.
     * @param e1
     * @param e2
     * @return
     */
     AdvancedEdge *bridgeBoundaries(AdvancedEdge *edge1, AdvancedEdge *edge2);

    /**
     * @brief unlinkTriangle
     * Unlinks a triangle from the mesh. O(1).

     * Resulting isolated vertices and edges are unlinked too.
     * If necessary, this method duplicates non-advanced vertices that can
     * occur due to the removal of the triangle.
     * The unlinked triangle, along with the other possible unlinked elements,
     * must be removed from the List T through removeUnlinkedElements().
     */
    void unlinkTriangle(AdvancedTriangle *triangle);

    /**
     * @brief unlinkTriangleNoManifold
     * Unlinks a triangle from the mesh. O(1).

     * No check is performed on the resulting topology, which may be inconsistent.
     * The unlinked triangle, along with the other possible unlinked elements,
     * must be removed from the List T through removeUnlinkedElements().
     */
    void unlinkTriangleNoManifold(AdvancedTriangle *triangle);


    /**
     * @brief removeTriangle
     * Removes a triangle from the mesh. O(N).

     * This is equivalent to an unlinkTriangle(t) followed by a
     * removeUnlinkedElements().
     * @param t
     */
    void removeTriangle(AdvancedTriangle *t)
    {
        unlinkTriangle(t);
        removeUnlinkedElements();
    }

    /**
     * @brief removeTriangles
     *  Removes all the unlinked triangles from List T.
     * Returns the number of removed triangles. O(N).
     * @return
     */
    int removeTriangles();

    /**
     * @brief removeEdges
     * Removes all the unlinked edges from List E.
     * Returns the number of removed edges. O(N).
     * @return
     */
    int removeEdges();

    /**
     * @brief removeVertices
     * Removes all the unlinked vertices from the mesh.
     * Remove extra vertices that does NOT belong to an edge, and just flying.
     * O(N).
     * @return
     * Returns the number of removed vertices from the mesh.
     */
    int removeVertices();

    /**
     * @brief removeUnlinkedElements
     * Removes all the unlinked elements from the lists.
     * Returns the number of removed elements. O(N).
     * @return
     */
    int removeUnlinkedElements()
    {
        return removeTriangles() + removeEdges() + removeVertices();
    }

    /**
     * @brief removeRedundantVertices
     * Removes all the vertices that can be deleted without changing the
     * geometric realization. O(N).
     * @return
     */
    int removeRedundantVertices();

    /**
     * @brief deselectTriangles
     * Deselects all the triangles. O(N).
     */
    void deselectTriangles();

    /**
     * @brief removeSelectedTriangles
     * Removes all the selected triangles. O(N).
     */
    void removeSelectedTriangles();

    /**
     * @brief selectBoundaryTriangles
     * Selects all the triangles having at least one boundary vertex. O(N).
     * Returns the number of selected triangles.
     * @return
     */
    int selectBoundaryTriangles();

    uint64_t getNumberBoundaryEdges();

    /**
     * @brief growSelection
     * Enlarges the current selection of one triangle in width. O(N).

     * Each triangle sharing at least one vertex with a currently selected
     * triangle becomes selected.
     * Returns the number of newly selected triangles.
     * @return
     */
    int growSelection();

    /**
     * @brief shrinkSelection
     * Shrinks the current selection of one triangle in width. O(N).

     * Each triangle sharing at least one vertex with a currently unselected
     * triangle becomes unselected.
     */
    void shrinkSelection();

    /**
     * @brief invertSelection
     * Inverts the selection status of all the triangles. O(N).

     * If 't0' is not nullptr, then only the connected component containing 't0'
     * is inverted.
     * @param t0
     */
    void invertSelection(AdvancedTriangle *inputTriangles = nullptr);

    /**
     * @brief reselectSelection
     * If 't0' is selected, deselect everything but the selected triangles
     * connected to 't0'
     * @param t0
     */
    void reselectSelection(AdvancedTriangle *t0);

    /**
     * @brief createSubMeshFromSelection
     * Creates a new Mesh out of an existing selection containing 't0'. O(output).

     * If necessary, non-advanced vertices are properly duplicated.
     * If 'keep_ref' is set to TRUE, then elements of the original mesh point
     * (through their info field) to corresponding elements of the newly created copy.
     * @param t0
     * @param keep_ref
     * @return
     */
     AdvancedMesh *createSubMeshFromSelection(AdvancedTriangle *selection = nullptr,
                                              bool keepReference = 0,
                                              const bool &verbose = true);

    /**
     * @brief createSubMeshFromTriangle
     * Creates a new Mesh out of an existing triangle 't0'. O(output).
     * @param t0
     * @return
     */
     AdvancedMesh *createSubMeshFromTriangle(AdvancedTriangle *t0);

    /**
     * @brief selectSphericalRegion
     * Marks all the triangles within distance L from 'p' as selected. O(output).

     * A triangle is considered to be within distance L from 'p' only if all
     * its three vertices are so.
     * Point 'p' is assumed to belong to triangle 't0', which is required to
     * limit the complexity to the size of the selected region.
     * Returns the number of selected triangles.
     * @param t0
     * @param L
     * @param p
     * @return
     */
    int  selectSphericalRegion(AdvancedTriangle *inputMesh, const double &distance, const AdvancedPoint *p);

    /**
     * @brief deselectSphericalRegion
     * Marks all the triangles within distance L from 'p' as deselected. O(output).

     * A triangle is considered to be within distance L from 'p' only if all
     * its three vertices are so.
     * Point 'p' is assumed to belong to triangle 't0', which is required to
     * limit the complexity to the size of the selected region.
     * Returns the number of deselected triangles.
     * @param t0
     * @param L
     * @param p
     * @return
     */
    int  deselectSphericalRegion(AdvancedTriangle *inputMesh, const double &distance, const AdvancedPoint *p);

    /**
     * @brief reselectSphericalRegion
     * Deselects all the triangles farther than L from 'p'. O(N).

     * A triangle is considered to be farther than L from 'p' if at least
     * one of its three vertices is so.
     * Point 'p' is assumed to belong to triangle 't0'. Passing 't0' avoids
     * the non robust and expensive computation of point-in-triangle.
     * @param t0
     * @param L
     * @param p
     */
    void reselectSphericalRegion(AdvancedTriangle *t0, const double &distance, const AdvancedPoint *p);

    /**
     * @brief retriangulateSelectedRegion
     * Re-triangulates the currently selected region using a Delaunay-like approach. O(SlogS).

     * A common plane is computed as the average of the planes of the triangles
     * selected; then, the vertices of the region are projected on the plane
     * and edges are iteratively swapped up to convergence (which is
     * guaranteed on planar and simple domains).
     * Finally, the vertices are moved back to their original positions.
     * This operation is particularly useful to improve the quality of nearly
     * flat regions.
     * The selection must be simple and its Gauss map must be entirely
     * contained in a semi-sphere.
     * Returns TRUE on success, FALSE otherwise.
     * @return
     */
    bool retriangulateSelectedRegion();

    /**
     * @brief isSelectionSimple
     * TRUE iff the set of selected triangles in 'l' is simply connected.
     * O(l->numels()).
     * @param l
     * @return
     */
    bool isSelectionSimple(List *selectionList);

    /**
     * @brief unmarkEverythingButSelections
     * Unmarks all the elements but leaves the selection status of triangles as is. O(N).
     */
    void unmarkEverythingButSelections();

    /**
     * @brief selectConnectedComponent
     * Selects all the triangles of the connected component containing t0. O(N).

     * If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
     * Returns the number of selected triangles.
     * @param t0
     * @param stop_on_sharp
     * @return
     */
    int selectConnectedComponent(AdvancedTriangle *t0, bool stop_on_sharp = 0);

    /**
     * @brief deselectConnectedComponent
     * Deselects all the triangles of the connected component containing t0. O(N).

     * If 'stop_on_sharp', expansion from 't0' brakes at tagged sharp edges.
     * Returns the number of deselected triangles.
     * @param t0
     * @param stop_on_sharp
     * @return
     */
    int deselectConnectedComponent(AdvancedTriangle *t0, bool stopOnSharp = 0);

    /**
     * @brief append
     * Append to the current mesh a copy of all the elements of 't'.
     * The newly created elements form a new selection.
     * @param t
     */
     void append(AdvancedMesh *input);

    /**
     * @brief split
     *
     * This method removes one connected component from the mesh and creates
     * a separate new mesh out of it. The components to be removed is the one
     * containing the first triangle in the list T.
     * Possible selection flags are deleted by this method.
     * @return
     */
     AdvancedMesh *split(const bool &verbose = true);

     /**
      * @brief splitPartitions
      * Splits a mesh with more than one parition into a group of meshes for further processing.
      * @return
      * A list of all the mesh paritions, each of them is an independent mesh.
      */
     std::vector< AdvancedMesh* > splitPartitions(const bool &verbose = true);

     /**
      * @brief appendMeshes
      * Appends a list of meshes to this mesh.
      * @param listMeshes
      * A list of meshes to be appended to this mesh.
      */
     void appendMeshes(std::vector< AdvancedMesh* > listMeshes);

     void toSimpleMesh(Mesh* mesh);

    /**
     * @brief getRegion
     * Make a list of triangles within distance L from 'p'. O(output).

     * Starting from 't0', which is assumed to contain 'p', add a triangle at a
     * time to the list as long as all the vertices stay within distance L from 'p'.
     * @param t0
     * @param L
     * @param p
     * @return
     */
    List* getRegion(AdvancedTriangle *inputMesh, const double &distance, const AdvancedPoint *p);

    /**
     * @brief removeRegion
     * Removes triangles within distance L from 'p'. O(N).

     * Starting from 't0', which is assumed to contain 'p', remove a triangle at a
     * time as long as all its vertices stay within distance L from 'p'.
     * @param t0
     * @param L
     * @param p
     */
    void removeRegion(AdvancedTriangle *inputTriangulaton, const double &distance, const AdvancedPoint *p);

    /**
     * @brief nextVertexOnRegionBoundary
     * Get the vertex next to 'v' on the boundary of the region. O(1).
     * @param v
     * @return
     */
    AdvancedVertex *nextVertexOnRegionBoundary(AdvancedVertex *vertex) const;

    /**
     * @brief getRegionInternalVertices
     * Retrieve internal vertices of a region. O(l->numels()).

     * This method returns a list containing an edge of the region's boundary
     * as its first element, and all the internal vertices as the remaining elements.
     * @param l
     * @return
     */
    List* getRegionInternalVertices(List *region);

    /**
     * @brief transformShell
     * Transform the vertices of the shell containing 't0' using the matrix m. O(S).
     * @param t0
     * @param m
     */
    void transformShell(AdvancedTriangle *inputTriangulation, const Matrix4x4& m);

    /**
     * @brief translate
     * Translate the mesh by a vector
     * @param t_vec
     */
    void translate(const AdvancedPoint& translationVector);

    /**
     * @brief getCenter
     * Return the center of mass of the mesh.
     * @return
     */
    AdvancedPoint getCenter() const;

    /**
     * @brief removeShell
     * Remove all the triangles belonging to the shell containing 'inputTriangles'. O(N).
     * @param inputTriangles
     */
    void removeShell(AdvancedTriangle *inputTriangles);

    /**
     * @brief sharpEdgeTagging
     * Tag sharp edges based on threshold angle. O(N).

     * Tag as sharp all the edges in which the normals of the two incident
     * triangles form an angle greater than 't'.
     * @param t
     */
    void  sharpEdgeTagging(const double taggingAngle);

    /**
     * @brief unmarkEverything
     * Unmark all the elements. O(N).
     */
    void   unmarkEverything();

    /**
     * @brief getBoundingBox
     * Bounding box longest edge. 'b' and 't' are set as the longest diagonal end-points. O(N).
     * @param b
     * @param t
     * @return
     */
    double getBoundingBox(AdvancedPoint& pMin, AdvancedPoint& pMax) const;

    /**
     * @brief bboxLongestDiagonal
     Bounding box longest diagonal. O(N).
     * @return
     */
    double bboxLongestDiagonal()
    {
        AdvancedPoint a, b;
        getBoundingBox(a, b);
        return a.distance(b);
    }

    /**
     * @brief getBoundingBallRadius
     * Approximate bounding ball radius. O(N).
     * @return
     */
    double getBoundingBallRadius() const;

    /**
     * @brief area
     * Total area of the mesh. O(N).
     * @return
     */
    double area() const;

    /**
     * @brief volume
     * Total volume of the mesh assuming that boundaries() = 0. O(N).
     * @return
     */
    double volume() const;

    /**
     * @brief normalize
     * Scale the mesh to make it fit within a cube [0,0,0]-[s,s,s]. O(N).
     * @param s
     */
    void normalize(const double s = 1.0);

    /**
     * @brief quantize
     * Scale the mesh to make it fit within a cube [0,0,0]-[s,s,s] and snap
     * coordinates on grid points O(N).
     * @param s
     */
    void quantize(const int s = 65536);

    /**
     * @brief transform
     * Transform the mesh geometry using the transformation matrix m. O(N).
     * @param m
     */
    void transform(const Matrix4x4& m);

    /**
     * @brief addNormalNoise
     * Randomly move vertices along their normals. O(N).
     * Displacement is bounded by 'p'% of the bounding ball radius.
     * @param p
     */
    void addNormalNoise(const double p);

    /**
     * @brief iterativeEdgeSwaps
     * Iteratively swaps edges to minimize the Delaunay minimum angle. O(N).

     * Edges tagged as sharp are constrained not to swap.
     * On generically curved advanceds this process is not guaranteed to converge.
     * This method returns TRUE if convergence is reached, FALSE otherwise.
     * @return
     */
    bool iterativeEdgeSwaps();

    /**
     * @brief isInnerPoint
     * True if the mesh properly contains 'p' (exact when using rationals).
     * The mesh is assumed to be a well-defined polyhedron (i.e. no boundary,
     * no degenerate triangles, no self-intersections) and to have a correct
     * orientation.
     * result is undetermined otherwise.
     * @param p
     * @return
     */
    bool isInnerPoint(AdvancedPoint& p) const;

    /**
     * @brief loopSubdivision
     * This function performs one step of loop subdivision over the entire mesh.
     * Noe that if 'midPoint' is set, only the connectivity is subdivided, whereas the surface
     * shape is kept unchanged.
     *
     * @param midPoint
     * If this flag is set, the connectivity is subdivided, whereas the surface shape itself is kept
     * unchanged.
     */
    void loopSubdivision(const bool &midPoint = false);


    /**
     * @brief flipNormals
     * Invert all the triangle and edge orientations. O(N).
     */
    void flipNormals();

    /**
     * @brief flipNormals
     * Invert the orientation of triangles and edges belonging to the shell
     * containing 't0'. O(S).
     * @param t0
     */
    void flipNormals(AdvancedTriangle *t0);

    /**
     * @brief topTriangle
     * Return the triangle with the maximum 'z' coordinate in the shell
     * containing 't0'. O(N).

     * Useful for orienting meshes bounding solids.
     * @param t0
     * @return
     */
    AdvancedTriangle *topTriangle(AdvancedTriangle *t0);

    /**
     * @brief eulerUpdate
     * Updates the values of nBoundaries, nHandles and nShells.
     * Complexity O(N).

     * The relative dirty bits are set to zero.
     */
    void eulerUpdate(const bool &verbose = true);

    /**
     * @brief openToDisk
     * Duplicates edges and vertices to make the mesh homeomorphic to a disk.
     * O(N).
     */
    void openToDisk();

    /**
     * @brief removeSmallestComponents
     * Removes all connected components but the one having most triangles.
     * Returns the number of components removed.
     * @return
     */
    int removeSmallestComponents();

    /**
     * @brief forceNormalConsistence
     * Checks that triangles are consistently oriented and, if they are not,
     * invert some of them to achieve an overall consistency. If the mesh is
     * not orientable cut it. Returns: 0 if mesh was already oriented; 1 if
     * the mesh could be oriented without cuts; >1 if cuts were necessary.
     * @return
     */
    int forceNormalConsistence();

    /**
     * @brief forceNormalConsistence
     * Checks that triangles are consistently, but acts on a single connected component and uses
     * one specific triangle from which the orientation is propagated.
     * @return
     */
    int forceNormalConsistence(AdvancedTriangle *inputTriangle);

    /**
     * @brief duplicateNonManifoldVertices
     * This function detects singular vertices in the mesh and duplicte them.
     * If a vertex is topologically non-manifold, this data structure does not guarantee its
     * functionality. Therefore, in order to use the same triangle mesh, this method allows to
     * duplicate such vertices.
     * Note that the data-structure cannot code non-manifold edges.
     *
     * @return
     * Returns the number of singular vertices after being duplicated.
     */
    int duplicateNonManifoldVertices();

    /**
     * @brief removeDuplicatedTriangles
     * Remove redundant triangles (i.e. having the same vertices as others)
     * and return their number.
     * @return
     */
    int removeDuplicatedTriangles();

    /**
     * @brief checkConnectivity
     * This function checks and validates the connectivity of the mesh. If everything is fine in
     * the mesh , a nullptr is returned, otherwise the issue will be returned.
     * This method should be used when implementing new algorithms to check the consistency of
     * the connectivity graph. Because an inconsistent graph is not assumed by all the other
     * methods, such a flaw is considered critical and the program should terminate.
     *
     * @return
     * If connectivity is ok, nullptr is returned, otherwise a string describing the error is
     * returned.
     */
    const char *checkConnectivity();

    /**
     * @brief fixConnectivity
     * If called in rebuildConnectivity(bool) fix the connectivity between
     * geometric elements
     * @return
     */
    bool fixConnectivity();

    /**
     * @brief rebuildConnectivity
     * Considers triangles as purely geometric entities and recomputes their
     * connectivity based on vertex correspondence.
     * Returns false if mesh was not an oriented advanced.
     * @return
     */
     bool rebuildConnectivity(bool fixMeshConnectivity = true);

    /**
     * @brief mergeCoincidentEdges
     * Looks for topologically different edges having the same geometry (i.e. coincident vertices)
     * and unify them.
     *
     * @return
     * Return the number number.
     */
    int mergeCoincidentEdges();

    /**
     * @brief removeDegenerateTriangles
     * Looks for zero-area triangles and resolve them by splits and collapses.
     * If, for any reason, some degeneracies cannot be removed
     * a warning is issued and the degenerate triangles that could not
     * be resolved are selected.
     * The absolute value of the integer returned is the number of
     * collapses performed; the return value is negative if some
     * degenerate triangles could not be resolved.
     * @return
     */
    int removeDegenerateTriangles();

    /**
     * @brief strongDegeneracyRemoval
     * Calls 'removeDegenerateTriangles()' and, if some degeneracies remain,
     * removes them and fills the resulting holes. Then tries again and, if
     * some degeneracies still remain, removes them and their neighbors and
     * fills the resulting holes, and so on, until the neighborhood growth
     * reaches maxIterations. If even in this case some degeneracies remain,
     * returns false, otherwise returns true.
     * @param maxIterations
     * @return
     */
    bool strongDegeneracyRemoval(int maxIterations);

    /**
     * @brief strongIntersectionRemoval
     *  Removes all the intersecting triangles and patches the resulting holes.
     * If the patches still produce intersections, iterates again on a larger
     * neighborhood. Tries up to max_iters times before giving up. Returns
     * true only if all the intersections could be removed.
     * @param max_iters
     * @return
     */
    bool strongIntersectionRemoval(int max_iters);

    /**
     * @brief cleanMesh
     * Iteratively call strongDegeneracyRemoval and strongIntersectionRemoval
     * to produce an eventually clean mesh without degeneracies and intersections.
     * The two aforementioned methods are called up to maxIterations times and
     * each of them is called using 'innerLoops' as a parameter.
     * Returns true only if the mesh could be completely cleaned.
     * @param maxIterations
     * @param innerLoops
     * @return
     */
    bool cleanMesh(uint64_t maxIterations = 25, int innerLoops = 3);

    /**
     * @brief ensureWatertightness
     * Ensures that this mesh is watertightight with no holes, no self
     * intersections and two-advanced. This function must be called before
     * optimizing the mesh to avoid failure during the optimization process.
     * @return True if the optimization works and false if it fails.
     */
    bool ensureWatertightness();

    /**
     * @brief removeBoundaryTriangles
     * Removes the boundary triangles to ensure that the watertightness is
     * truely achieved.
     * @return
     */
    int removeBoundaryTriangles();

    /**
     * @brief removeOverlappingTriangles
     * Removes overlapping triangles and return their number.
     * If possible, swap edges to remove overlaps. When it is not enough, remove the
     * overlapping triangles from the mesh.
     * @return
     * Returns the number of triangles that was necessary to remove.
     */
    int removeOverlappingTriangles();

    /**
     * @brief checkGeometry
     * This function checks the triangulation of the geometry. It validates the mesh for
     * degeneracies, concident vertices and overlaps.
     * If such a flaw is found returns its the closest vertex.
     * If no issue is detected it return a nullptr.
     *
     * @return
     * If something is wrong in the mesh, it returns the closest vertex where the issue happened.
     */
    AdvancedVertex *checkGeometry();

    /**
     * @brief selectIntersectingTriangles
     * Selects all the triangles that unproperly intersect other parts of
     * the mesh and return their number. The parameter 'tris_per_cell'
     * determines the depth of the recursive space subdivision used to keep
     * the complexity under a resonable threchold. The default value is safe
     * in most cases.
     * if 'justproper' is true, coincident edges and vertices are not regarded
     * as intersections even if they are not common subsimplexes.
     * @param tris_per_cell
     * @param justproper
     * @return
     */
    int selectIntersectingTriangles(uint16_t trianglesPerCell = 50,
                                    bool justProper = false);

    /**
     * @brief safeCoordBackApproximation
     * This is as coordBackApproximation() but it also checks for
     * intersections and, if any, it tries different approximations.
     * Returns true if no intersections remain.
     * @return
     */
    bool safeCoordBackApproximation();

    /**
     * @brief removeSmallestComponents
     * Removes all the connected components whose area is less than 'epsilon'.
     * Returns the number of components removed.
     * @param epsilon
     * @return
     */
    int removeSmallestComponents(double epsilonArea);

    /**
     * @brief StarTriangulateHole
     * Computes the barycenter of the boundary loop and connects it with
     * all the edges of the loop. Returns the number of triangles
     * created.
     * @return
     */
     int starTriangulateHole(AdvancedEdge *);

    /**
     * @brief TriangulateHole
     * Creates a triangulation whose projection on the plane with normal 'nor'
     * is Delaunay. Returns the number of triangles created.
     * @param nor
     * @return
     */
    int TriangulateHole(AdvancedEdge *edge, AdvancedPoint *normal);

    /**
     * @brief TriangulateHole
     * Creates a triangulation of the hole based on heuristics.
     * Returns the number of triangles created.
     * @return
     */
    int TriangulateHole(AdvancedEdge *edge);

    /**
     * @brief TriangulateFlatHole
     * Creates a triangulation of the hole based on the assumption that it is flat.
     * Returns the number of triangles created.
     * @return
     */
     // int TriangulateFlatHole(Edge *);

    /**
     * @brief TriangulateHole
     *  Creates a triangulation and inserts the additional points in the List.
     * Returns the number of triangles created.
     * @return
     */
     int TriangulateHole(AdvancedEdge*edge, List*vertexList);

    /**
     * @brief fillHole
     *  Hole filling algorithm. Performs a triangulation based on heuristics and,
     * if 'refine' is set to true, adds inner vertices to reproduce the sampling
     * density of the surroundings.
     * @param refine
     */
    void fillHole(AdvancedEdge *edge, bool refine = true);

    /**
     * @brief fillHoles
     * Fills all the holes in the mesh having at least 'minNumberBoundaryEdges'
     * boundary edges. If minNumberBoundaryEdges is 0, it will fill all the
     * holes in the mesh. This function is extremely important to use before
     * processing any given mesh to avoid optimization issues.
     * If 'refinePatches' is true, it adds inner vertices to reproduce the
     * sampling density of the surroundings.
     * If 'numberEdges' is 0 (default), all the holes are patched.
     * @param minNumberBoundaryEdges
     * The minimum number of boundary edges in a hole that must be filled.
     * @param refinePatches
     * A flag to add inner vertices to reporoduce the sampling density of the
     * surrounding triangles.
     * @return
     * Returns number of holes patched.
     */
    uint64_t fillHoles(const uint64_t minNumberBoundaryEdges = 0,
                       const bool refinePatches = true);

    /**
     * @brief refineSelectedHolePatches
     * Takes a selected region and inserts inner vertices to reproduce the
     * sampling density of the surroundings. If 't0' is not nullptr, only the
     * selected region containing 't0' is refined. Returns the number of
     * vertices inserted.
     * @param t0
     * @return
     */
     int refineSelectedHolePatches(AdvancedTriangle *t0 =nullptr);

    /**
     * @brief retriangulateVT
     * Retriangulates the vertex neghborhood based on heuristics.
     * @return
     */
    int retriangulateVT(AdvancedVertex *vertex);

    /**
     * @brief joinBoundaryLoops
     * Joins the two boundary vertices gv and gw through an edge.
     * A pair of triangles is added to properly change the topology of the mesh.
     * @return
     */
    /**
     * @brief joinBoundaryLoops
     * Joins the two boundary vertices gv and gw through an edge. A pair of triangles is added to
     * properly change the topology of the mesh. On success, the return value is the new edge
     * connecting the two vertices. nullptr is returned on failure.
     * If 'justconnect' is false, the remaining hole is filled with new triangles, unless gv and gw
     * are contiguous on the boundary loop (failure).
     * If 'justConnect' is true, gv and gw must not belong to the same boundary loop (failure).
     * If 'refine', the patching triangles are refined to reproduce neighboring density.
     *
     * @param justConnect
     * @return
     */
    AdvancedEdge* joinBoundaryLoops(AdvancedVertex *gv, AdvancedVertex *gw,
                                    const bool &justConnect, const bool &refine);

    /**
     * @brief printReport
     */
    void printReport();

    /**
     * @brief getNumberSingularEdges
     * @return
     */
    int getInitialNumberSingularEdges() const { return _loadedNumberSingularEdges; }

    /**
     * @brief getInitialNumberNonManifoldvertices
     * @return
     */
    int getInitialNumberNonManifoldvertices() const { return _loadedNonManifoldVertices; }

    /**
     * @brief getInitialNumberDuplicateTriangles
     * @return
     */
    int getInitialNumberDuplicateTriangles() const { return _loadedDuplicateTriangles; }

    /**
     * @brief getInitialNumberFloatingVertices
     * @return
     */
    int getInitialNumberFloatingVertices() const { return _loadedNumberFloatingVertices; }

private:

    /**
     * @brief cotangentAngle
     * @param pivot
     * @param a
     * @param b
     * @return
     */
    double _cotangentAngle(const Vector3f& pivot, const Vector3f& a, const Vector3f& b);

    /**
     * @brief _buildIndexedMesh
     * @param mesh
     * @return
     */
    AdvancedIndexedMesh _buildIndexedMesh(AdvancedMesh& mesh);

    /**
     * @brief _computeNeighborhoods
     * @param mesh
     * @param vertexNeighbors
     * @param faceNeighbors
     */
    void _computeNeighborhoods(AdvancedIndexedMesh& mesh,
                               Neighborhood& vertexNeighbors, Neighborhood& faceNeighbors);

    /**
     * @brief _computeKernel
     * @param mesh
     * @param vertexIndex
     * @param verticesN
     * @param facesN
     * @return
     */
    Vector3f _computeKernel(AdvancedIndexedMesh& mesh, const uint32_t &vertexIndex,
                            Neighbors& verticesN, Neighbors& facesN);

    /**
     * @brief _computeCotangentWeight
     * @param mesh
     * @param vertexIndex
     * @param neighborIndex
     * @param faceN
     * @return
     */
    float _computeCotangentWeight(AdvancedIndexedMesh& mesh,
                                  const uint32_t &vertexIndex, const uint32_t &neighborIndex,
                                  Neighbors& faceN);

    /**
     * @brief _smoothVertex
     * @param mesh
     * @param vertexIndex
     * @param kernel
     * @param smoothParam
     * @param inflateParam
     * @return
     */
    AdvancedVertex _smoothVertex(AdvancedIndexedMesh& mesh,
                                 const uint32_t vertexIndex, const Vector3f kernel,
                                 const float &smoothParam, const float &inflateParam);

    /**
     * @brief _convertVertexToVector3f
     * @param vertex
     * @return
     */
    Vector3f _convertVertexToVector3f(const AdvancedVertex* vertex);

public:

    /**
     * @brief applyLaplacianSmooth
     * @param numIterations
     * @param smoothLambda
     * @param inflateMu
     */
    void applyLaplacianSmooth(const uint64_t &numIterations,
                              const float &smoothLambda, const float &inflateMu);

protected:

    /**
     * @brief _initialValuesUpdated
     */
    bool _initialValuesUpdated;

    /**
     * @brief _loadedNumberSingularEdges
     */
    int _loadedNumberSingularEdges;

    /**
     * @brief _loadedNumberFloatingVertices
     */
    int _loadedNumberFloatingVertices;

    /**
     * @brief _loadedOrientationFlag
     */
    int _loadedOrientationFlag;

    /**
     * @brief _loadedNonManifoldVertices
     */
    int _loadedNonManifoldVertices;

    /**
     * @brief _loadedDuplicateTriangles
     */
    int _loadedDuplicateTriangles;

protected:

    /**
     * @brief _watsonInsert
     * @return
     */
    AdvancedVertex* _watsonInsert(AdvancedPoint *point, List *delaunayTriangulation, int numberTriangles);
};

// Defines
#include <data/meshes/advanced/Macros.hh>
}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_MESH_H
