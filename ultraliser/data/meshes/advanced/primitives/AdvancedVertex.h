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

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_VERTEX_H
#define ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_VERTEX_H

#include <common/Common.h>
#include <data/structures/Heap.h>
#include <data/structures/List.h>
#include <data/meshes/advanced/primitives/AdvancedPoint.h>

namespace Ultraliser
{

// Forward declarations
class AdvancedEdge;
class AdvancedTriangle;

/**
 * @brief The AdvancedVertex class
 * This AdvancedVertex is part of the AdvancedMesh.
 * AdvancedMesh is based on TMesh of MeshFix, which preserves the connectivity information.
 *
 * This class represents a vertex of an oriented triangulation. The base-class AdvacedPoint
 * describes geometrical and additional attributes of the vertex. The field 'e0' is sufficient to
 * retrieve all the neighboring elements in optimal time, while the field 'mask' is useful for
 * assigning up to 256 different states to the vertex.
 */
class AdvancedVertex : public AdvancedPoint
{
public :

    /**
     * @brief e0
     * One of the incident edges
     */
    AdvancedEdge *e0;

    /**
     * @brief mask
     * bit-mask for marking purposes
     */
    unsigned char mask;

    /**
     * @brief Vertex
     * Constructor
     * Creates a new vertex with coordinates (0, 0, 0).
     */
    AdvancedVertex();

    /**
     * @brief Vertex
     * Constructor
     * Creates a new vertex with the same coordinates (x, y, z).
     *
     * @param x
     * X-coordinate
     * @param y
     * Y-coordinate
     * @param z
     * Z-coordinate
     */
    AdvancedVertex(const double &x, const double &y, const double &z);


    /**
     * @brief Vertex
     * Constructor
     * Creates a new vertex with the same coordinates as 'p'. The info field is not copied.
     *
     * @param p
     * Input point
     */
    AdvancedVertex(const AdvancedPoint *p);

    /**
     * @brief Vertex
     * Constructor
     * Creates a new vertex with the same coordinates as 'p'. The info field is not copied.
     *
     * @param p
     * Input point
     */
    AdvancedVertex(const AdvancedPoint &p);
    ~AdvancedVertex();

    /**
     * @brief isBaseType
     * Returns true only if object is a basic Vertex. All the reimplementations must return false.
     *
     * @return
     * Returns true only if object is a basic Vertex. All the reimplementations must return false.
     */
    bool isBaseType() const { return true; }

    /**
     * @brief isLinked
     * TRUE iff vertex is not isolated.
     * @return
     */
    bool isLinked() const { return (e0 != 0); }

    /**
     * @brief getAdjacentVertices
     * List of adjacent vertices.
     *
     * Returns the list of vertices which are linked to this through an edge.
     * The  list is counter-clockwise ordered. In the case of an internal
     * vertex the list starts from the opposite vertex of e0. If the vertex is on
     * the  boundary,  the  list starts from the opposite vertex of the clock-
     * wise-most boundary edge.
     */
    List *getAdjacentVertices() const;

    /**
     * @brief getAdjacentEdges
     * Returns the list of all the incident edges over this vertex.
     *
     * @return
     * Returns the list of edges incident at this vertex. The list is
     * counter-clockwise ordered.
     * In the case of an internal vertex the list starts from 'e0'.
     * If the vertex is on the boundary, the list starts from its clockwise-most incident boundary
     * edge.
     */
    List *getIncidentEdges() const;

    /**
     * @brief VT
     * List of incident triangles.
     * Returns  the  list  of triangles incident at this. The list is counter-clockwise ordered.
     * In the case of an internal vertex the list starts from the triangle
     * on the left of e0, when looking from this.
     * If the vertex is on the boundary, the list starts from the clockwise-most
     * boundary triangle.
     *
     * @return
     * Returns  the  list  of triangles incident at this. The list is counter-clockwise ordered.
     */
    List* VT() const;

    /**
     * @brief getEdge
     * Returns the edge connecting this vertex to 'v'. NULL if such an edge does not exist.
     *
     * @param v
     * A given vertex.
     * @return
     * Returns the edge connecting this vertex to 'v'. NULL if such an edge does not exist.
     */
    class AdvancedEdge *getEdge(const AdvancedVertex *v) const;

    /**
     * @brief valence
     * Returns the number of incident edges.
     *
     * @return
     * Returns the number of incident edges.
     */
    int valence() const;

    /**
     * @brief isOnBoundary
     * TRUE iff vertex is on boundary.
     *
     * @return
     * TRUE iff vertex is on boundary.
     */
    int isOnBoundary() const;

    /**
     * @brief nextBoundaryEdge
     * Returns the edge following this vertex on the boundary.
     * This edge is the counterclockwise-most incident edge.
     * Returns nullptr if this vertex is not on the boundary.
     *
     * @return
     * Returns the edge following this vertex on the boundary, nullptr otherwise.
     */
    AdvancedEdge *nextBoundaryEdge() const;

    /**
     * @brief prevBoundaryEdge
     * Returns the edge preceeding this vertex on the boundary.
     * This edge is the clockwise-most incident edge.
     * Returns nullptr if this vertex is not on the boundary.
     *
     * @return
     * Returns the edge preceeding this vertex on the boundary, nullptr otherwise.
     */
    AdvancedEdge *prevBoundaryEdge() const;

    /**
     * @brief nextOnBoundary
     * Returns the vertex following this one on the boundary.
     * If the vertex is on the boundary, this is equivalent to nextBoundaryEdge()->oppositeVertex(v)
     * otherwise returns otherwise.
     * @return
     * Returns the vertex following this one on the boundary, nullptr otherwise.
     */
    AdvancedVertex *nextOnBoundary() const;

    /**
     * @brief prevOnBoundary
     * Returns the vertex preceeding this one on the boundary.
     * If the vertex is on the boundary, this is equivalent to prevBoundaryEdge()->oppositeVertex(v)
     * otherwise returns nullptr otherwise.
     *
     * @return
     * Returns the vertex preceeding this one on the boundary, nullptr otherwise.
     */
    AdvancedVertex *prevOnBoundary() const;

    /**
     * @brief isFlat
     * TRUE iff vertex neighborhood is a flat disk. Always FALSE for boundary vertices.
     *
     * @return
     * TRUE iff vertex neighborhood is a flat disk. Always FALSE for boundary vertices.
     */
    bool isFlat() const;

    /**
     * @brief isDoubleFlat
     * TRUE iff vertex neighborhood is made of either two flat halfdisks or one
     * flat halfdisk on a rectilinear boundary.
     * When TRUE, the two edges may be initialized with the two non-flat
     * incident edges.
     * If the two halfdisks are also coplanar, returns TRUE and e1 and e2 are set to nullptr.
     *
     * @param e1
     * @param e2
     * @return
     */
    bool isDoubleFlat(AdvancedEdge **e1, AdvancedEdge **e2) const;

    /**
     * @brief removeIfRedundant
     * Unlinks the vertex if it is either Flat() or DoubleFlat(). On success,
     * the function returns TRUE, the vertex neighborhood is retriangulated,
     * and the geometric realization does not change.
     * If 'check_neighborhood' is 'false', incident triangles are assumed to
     * be neither degenerate nor overlapping.
     * Attention! The method unlinks the vertex but does not remove it from
     * the Basic_TMesh list.
     * The calling function must clear the lists through
     * TMesh::removeUnlinkedElements().
     * @param check_neighborhood
     * @return
     */
    bool removeIfRedundant(bool checkNeighborhood = true);

    /**
     * @brief getNormal
     * Normal at the vertex computed as the sum of incident triangle normals
     * weighted on their incidence angle.
     * @return
     */
    AdvancedPoint getNormal() const;

    /**
     * @brief getBoundaryAngle
     * Returns  the angle between the two incident boundary edges.
     * If the vertex is not on the boundary, returns -1.
     * @return
     */
    double getBoundaryAngle() const;

    /**
     * @brief getAngleForTriangulation
     * Discriminant Angle for triangulating 3D polygons.
     * This method is useful when patching holes, and represents  a  heuristic
     * for choosing which vertex of the hole's boundary must be patched first.
     * Several cases are considered, including degenerate ones. Let 'v1' and 'v2'
     * be  the two boundary vertices which are linked to this one through a boundary
     * edge. If 'v1' and 'v2' coincide, the method returns a negative  number.
     * If  'v1'  ,  'this'  ,  'v2'  form  a flat angle, the method returns 3PI (270
     * degrees). If the angle formed by 'v1' ,  'this'  ,  'v2'  is  0,  the  method
     * returns 0. If the vertex is not on the boundary, the method returns the
     * limit number DBL_MAX. In all the other cases the method returns the sum
     * of three angles A + D1 + D2, where A is the angle formed by v1 , this ,
     * v2 , while D1 is the angle between the  normal  of  the  clockwise-most
     * incident  boundary  triangle and the normal of the triangle v1 , this ,
     * v2; D2 is the analogous for the counterclockwise-most incident boundary
     triangle.
     * @return
     */
    double getAngleForTriangulation() const;

    /**
     * @brief getAngleOnAveragePlane
     * Discriminant Angle for triangulating flat (or nearly flat) polygons.
     * This  method  returns the angle between the two incident boundary edges
     * when projected onto the plane whose normal is 'n'. This  angle
     * may be more than PI, because it represents the aperture of the non-tri-
     * angulated region around the vertex when projected on the plane. If  the
     * vertex  is  not  on  the  boundary,  the method returns the limit value
     * DBL_MAX.
     * @param n
     * @return
     */
    double getAngleOnAveragePlane(AdvancedPoint *normal) const;

    /**
     * @brief totalAngle
     *  Sum of incident angles. Returns -1 if on boundary.
     * @return
     */
    double totalAngle() const;

    /**
     * @brief gaussianCurvature
     * Excess angle. Returns DBL_MAX if on boundary.
     * @return
     */
    double gaussianCurvature() const { double t = totalAngle(); return (t >= 0) ? (t) : (DBL_MAX); }

    /**
     * @brief totalDihedralAngle
     * Sum of signed dihedral angles. Returns DBL_MAX if on boundary.
     * @return
     */
    double totalDihedralAngle() const;

    /**
     * @brief voronoiArea
     * A third of the total area of incident triangles.
     * @return
     */
    double voronoiArea() const;

    /**
     * @brief zip
     * Zips the gap starting from here.
     * @return
     */
    int zip(const bool =1);

    /**
     * @brief inverseCollapse
     * @return
     */
    AdvancedEdge *inverseCollapse(AdvancedVertex *,
                                  AdvancedEdge *, AdvancedEdge *, AdvancedEdge *,
                                  AdvancedEdge *, AdvancedEdge *,
                                  AdvancedTriangle *, AdvancedTriangle *);
};

#define TO_VERTEX(POINTER)      static_cast< AdvancedVertex* > (POINTER)
#define TO_EDGE(POINTER)        static_cast< AdvancedEdge* > (POINTER)

// Scans the nodes 'n' of a list 'l' of vertices 'v'.
#define FOR_EACH_VV_VERTEX(l, v, n)                                                                 \
    for (n = l->head(), v = (n) ? (TO_VERTEX(n->data)) : nullptr;                                   \
         n != nullptr; n = n->next(), v = (n) ? (TO_VERTEX(n->data)) : nullptr)

// Scans the nodes 'n' of a list 'l' of edges 'e'.
#define FOR_EACH_VE_EDGE(l, e, n)                                                                   \
    for (n = l->head(), e = (n) ? (TO_EDGE(n->data)) : nullptr;                                     \
         n != nullptr; n = n->next(), e = (n) ? (TO_EDGE(n->data)) : nullptr)

// Scans the nodes 'n' of a list 'l' of triangles 't'.
#define FOR_EACH_VT_TRIANGLE(l, t, n)                                                               \
    for (n = l->head(), t = (n) ? ((AdvancedTriangle*) n->data) : nullptr;                          \
         n != nullptr; n = n->next(), t = (n) ? ((AdvancedTriangle*) n->data) : nullptr)

 /**
 * @brief The ExtendedVertex class
 * Extended vertex for temporary use during connectivity creation.
 * This class is used to allow the reconstruction of the connectivity
 * in linear time (average case) and to handle badly oriented input files.
 * It provides a complete VE relation.
 */
class ExtendedVertex
{
public:

    /**
     * @brief v
     */
    AdvancedVertex *v;

    /**
     * @brief VE
     */
    List VE;

    /**
     * @brief ExtVertex
     * @param a
     */
    ExtendedVertex(AdvancedVertex* a)
    {
        v = a;
    }
};

}

#endif //ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_VERTEX_H

