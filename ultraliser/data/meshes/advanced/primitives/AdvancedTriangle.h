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

#pragma once

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_TRIANGLE_H
#define ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_TRIANGLE_H

#include "AdvancedEdge.h"

namespace Ultraliser
{

/**
 * @brief The Triangle class
 * Triangle of a TMesh.
 *
 * This  class represents a triangle of a triangulation. Each Triangle has
 * an orientation (clockwise or counter-clockwise) due  to  the  order in
 * which  its  edges are stored in the class. When looking at the triangle
 * so that (e1, e2, e3) are sorted counter-clockwise, the  normal  at  the
 * triangle  points  towards  the  observer.  The field mask is useful for
 * assigning up to 256 different states to the edge.
 */
class AdvancedTriangle
{
public :

    /**
     * @brief edge1
     * First edge of the triangle.
     */
    AdvancedEdge *edge1;

    /**
     * @brief edge2
     * Second edge of the triangle.
     */
    AdvancedEdge *edge2;

    /**
     * @brief edge3
     * Third edge of the triangle.
     */
    AdvancedEdge *edge3;

    // Further information
    void *info;

    /**
     * @brief mask
     * Bit-mask for marking purposes
     */
    unsigned char mask;

    /**
     * @brief Triangle
     */
    AdvancedTriangle();

    /**
     * @brief Triangle
     */
    AdvancedTriangle(AdvancedEdge *, AdvancedEdge *, AdvancedEdge *);

    /**
     * @brief isBaseType
     * Returns true only if object is a basic Triangle.
     * All the reimplementations must return false.
     * @return
     */
    bool isBaseType() const { return true; }

    /**
     * @brief isLinked
     * TRUE if properly linked.
     * @return
     */
    bool isLinked() const { return (edge1 != nullptr); }

    /**
     * @brief invert
     * Inverts the orientation of the triangle.
     */
    void invert() { p_swap((void **)(&edge2), (void **)(&edge3)); }

    /**
     * @brief v1
     * First vertex
     * @return
     */
    AdvancedVertex *v1() const { return edge1->commonVertex(edge2); }

    /**
     * @brief v2
     * Second vertex
     * @return
     */
    AdvancedVertex *v2() const { return edge2->commonVertex(edge3); }

    /**
     * @brief v3
     * Third vertex
     * @return
     */
    AdvancedVertex *v3() const { return edge3->commonVertex(edge1); }

    /**
     * @brief t1
     * First adjacent triangle. nullptr if boundary.
     * @return
     */
    AdvancedTriangle *t1() const { return edge1->oppositeTriangle(this); }

    /**
    * @brief t2
    * Second adjacent triangle. nullptr if boundary.
    * @return
    */
    AdvancedTriangle *t2() const {return edge2->oppositeTriangle(this); }

    /**
     * @brief t3
     * Third adjacent triangle. nullptr if boundary.
     * @return
     */
    AdvancedTriangle *t3() const {return edge3->oppositeTriangle(this); }

    /**
     * @brief hasEdge
     * TRUE iff 'e' is an edge of the triangle.
     * @param e
     * @return
     */
    bool hasEdge(const AdvancedEdge *e) const {return (e == edge1 || e == edge2 || e == edge3); }

    /**
     * @brief hasVertex
     * TRUE iff 'v' is a vertex of the triangle.
     * @param v
     * @return
     */
    bool hasVertex(const AdvancedVertex *v) const
    {
        return (edge1->hasVertex(v) || edge2->hasVertex(v) || edge3->hasVertex(v));
    }

    /**
     * @brief oppositeEdge
     * Triangle's edge opposite to 'v'.
     * nullptr if 'v' is not a vertex of the triangle.
     * @param v
     * @return
     */
    AdvancedEdge *oppositeEdge(const AdvancedVertex *v) const
    {
        return ((!edge1->hasVertex(v)) ? (edge1) :
               ((!edge2->hasVertex(v)) ? (edge2) :
               ((!edge3->hasVertex(v)) ? (edge3) : (nullptr))));
    }

    /**
     * @brief oppositeTriangle
     * Adjacent triangle opposite to 'v'.
     * nullptr if 'v' is not a vertex of the triangle.
     * @param v
     * @return
     */
    AdvancedTriangle *oppositeTriangle(const AdvancedVertex *v) const
    {
        return ((!edge1->hasVertex(v)) ? (t1()) :
               ((!edge2->hasVertex(v)) ? (t2()) :
               ((!edge3->hasVertex(v)) ? (t3()):nullptr)));
    }

    /**
     * @brief oppositeVertex
     * Triangle's vertex opposite to 'e'.
     * nullptr if 'e' is not an edge of the triangle.
     * @param e
     * @return
     */
    AdvancedVertex *oppositeVertex(const AdvancedEdge *e) const
    {
        return (e == edge1) ? (v2()) :
              ((e == edge2) ? (v3()) :
              ((e == edge3) ? (v1()) : nullptr));
    }

    /**
     * @brief rightTriangle
     * Triangle adjacent to the next edge of 'e'.
     * nullptr if 'e' is not an edge of the triangle.
     * @param e
     * @return
     */
    AdvancedTriangle *rightTriangle(const AdvancedEdge *e) const
    {
        return (e == edge1) ? (t2()) :
              ((e == edge2) ? (t3()) :
              ((e == edge3) ? (t1()) : nullptr));
    }

    /**
     * @brief leftTriangle
     * Triangle adjacent to the previous edge of 'e'.
     * nullptr if 'e' is not an edge of the triangle.
     * @param e
     * @return
     */
    AdvancedTriangle *leftTriangle(const AdvancedEdge *e) const
    {
        return (e == edge1) ? (t3()) :
              ((e == edge2) ? (t1()) :
              ((e == edge3) ? (t2()) : nullptr));
    }

    /**
     * @brief nextEdge
     * Edge next to 'e' in the ordering or the triangle.
     * nullptr if 'e' is not an edge of the triangle.
     * @param e
     * @return
     */
    AdvancedEdge *nextEdge(const AdvancedEdge *e) const
    {
        return ((e == edge1) ? (edge2) :
               ((e == edge2) ? (edge3) :
               ((e == edge3) ? (edge1) : nullptr)));
    }

    /**
     * @brief prevEdge
     * Edge preceeding 'e' in the ordering or the triangle.
     * nullptr if 'e' is not an edge of the triangle.
     * @param e
     * @return
     */
    AdvancedEdge *prevEdge(const AdvancedEdge *e) const
    {
        return ((e == edge1) ? (edge3) :
               ((e == edge2) ? (edge1) : ((e == edge3) ? (edge2):nullptr)));
    }

    /**
     * @brief nextVertex
     * Vertex next to 'v' in the ordering or the triangle.
     * nullptr if 'v' is not a vertex of the triangle.
     * @param v
     * @return
     */
    AdvancedVertex *nextVertex(const AdvancedVertex *v) const
    {
        return (!edge1->hasVertex(v)) ? (v3()) :
              ((!edge2->hasVertex(v)) ? (v1()) :
              ((!edge3->hasVertex(v)) ? (v2()) : nullptr));
    }

    /**
     * @brief prevVertex
     * Vertex preceeding 'v' in the ordering or the triangle.
     * nullptr if 'v' is not a vertex of the triangle.
     * @param v
     * @return
     */
    AdvancedVertex *prevVertex(const AdvancedVertex *v) const
    {
        return (!edge1->hasVertex(v)) ? (v1()) :
              ((!edge2->hasVertex(v)) ? (v2()) :
              ((!edge3->hasVertex(v)) ? (v3()) : nullptr));
    }

    /**
     * @brief nextEdge
     * Edge next to 'e' in the ordering or the triangle.
     * nullptr if 'e' is not an edge of the triangle.
     * @param v
     * @return
     */
    AdvancedEdge *nextEdge(const AdvancedVertex *v) const
    {
        return (v == v1()) ? (edge2) :
              ((v == v2()) ? (edge3) :
              ((v == v3()) ? (edge1) : nullptr));
    }

    /**
     * @brief commonEdge
     * If this triangle shares an edge with 'b', then such an edge is returned. nullptr otherwise.
     * @param b
     * @return
     */
    AdvancedEdge *commonEdge(const AdvancedTriangle *b) const
    {
        return ((edge1 == b->edge1 || edge1 == b->edge2 || edge1 == b->edge3) ? (edge1) :
              (((edge2 == b->edge1 || edge2 == b->edge2 || edge2 == b->edge3) ? (edge2) :
              (((edge3 == b->edge1 || edge3 == b->edge2 || edge3 == b->edge3) ? (edge3) : nullptr)))));
    }

    /**
     * @brief commonVertex
     * If this triangle shares a vertex with 'b', then such a vertex is
     * returned, nullptr otherwise.
     * @param b
     * @return
     */
    AdvancedVertex *commonVertex(const AdvancedTriangle *b) const;

    /**
     * @brief replaceEdge
     * Replace edge 'a' with edge 'b' in the triangle and return TRUE.
     * If 'a' is not an edge of the triangle return FALSE.
     * @param a
     * @param b
     * @return
     */
    bool replaceEdge(const AdvancedEdge *a, AdvancedEdge *b)
    {if (edge1==a) edge1=b; else if (edge2==a) edge2=b; else if (edge3==a) edge3=b; else return 0; return 1; }

    /**
     * @brief checkAdjNor
     *  TRUE if the oriantation is consistent with the one of 't' OR if
     * this and 't' do not share any edge.
     * @param t
     * @return
     */
    bool checkAdjNor(const AdvancedTriangle *t) const;

    /**
     * @brief getVector
     Return a vector orthogonal to the plane of the triangle.
If triangle is degenerate return a null vector.
     * @return
     */
    AdvancedPoint  getVector() const;

    /**
     * @brief getCenter
     * Return the triangle's barycenter.
     * @return
     */
    AdvancedPoint  getCenter() const;

    /**
     * @brief getCircleCenter
     * Return the center of the triangle's bounding sphere.
     * @return
     */
    AdvancedPoint  getCircleCenter() const;

    /**
     * @brief inSphere
     * TRUE iff 'p' is inside the triangle's bounding sphere.
     * @param p
     * @return
     */
    bool inSphere(const AdvancedPoint *p) const;

    /**
     * @brief squaredDistanceFromPoint
     * Squared distance of 'p' from the plane of the triangle.
     * Return -1 if triangle is degenerate.
     * @param p
     * @return
     */
    double squaredDistanceFromPoint(const AdvancedPoint *p) const;

    /**
     * @brief pointTriangleSquaredDistance
     * Squared distance of 'p' from the closest point of the triangle.
     * Return -1 if triangle is degenerate.
     * If closest point is in the interior of the triangle, *closest_edge
     * and *closest_vertex are set to nullptr.
     * If the closest point is in the interior of an edge, *closest_edge
     * is initialized with that edge.
     * If the closest point is a vertex, *closest_vertex is initialized with it.
     * @param p
     * @param closest_edge
     * @param closest_vertex
     * @return
     */
    double pointTriangleSquaredDistance(const AdvancedPoint *p,
                                       AdvancedEdge **closestEdge = nullptr,
                                       AdvancedVertex **closestVertex = nullptr) const;

    /**
     * @brief project
     * Projection of 'p' on the plane of the triangle.
     * Return INFINITE_POINT if triangle is degenerate.
     * @param p
     * @return
     */
    AdvancedPoint project(const AdvancedPoint *p) const;

    /**
     * @brief getLongestEdge
     * Returns the longest edge of the triangle.
     * @return
     */
    AdvancedEdge *getLongestEdge() const;

    /**
     * @brief isExactlyDegenerate
     * Degeneracy check using exact predicates.
     * Return TRUE iff triangle has zero area.
     * @return
     */
    bool isExactlyDegenerate() const;

    /**
     * @brief getCapEdge
     * Returns the edge opposite to a 'cap' vertex
     * @return
     */
    AdvancedEdge *getCapEdge() const;

    /**
     * @brief printTriangle
     * Print the coordinates of the three vertices to the file handler pointed
     * to by 'f' (stdout by default).
     * @param f
     */
    void printTriangleToFile(FILE *file =stdout) const;

    /**
     * @brief intersects
     * True if this triangle itersects 't' other than on common subsimplexes
     * if 'justproper' is true, coincident edges and vertices are not regarded
     * as intersections even if they are not common subsimplexes.
     * @param t
     * @param justproper
     * @return
     */
    bool intersects(const AdvancedTriangle *t, bool justproper =false) const;

    /**
     * @brief getNormal
     * Return a normal vector with direction (v1-v2) cross (v2-v3).
     * If triangle is degenerate return a null vector.
     * @return
     */
    AdvancedPoint  getNormal() const;

    /**
     * @brief area
     * Area of the triangle (Heron's formula).
     * @return
     */
    double area() const;

    /**
     * @brief perimeter
     * Perimeter of the triangle.
     * @return
     */
    double perimeter() const;

    /**
     * @brief getAngle
     * Angle at vertex 'v'. Return -1 if 'v' is not a vertex of the triangle.
     * @param v
     * @return
     */
    double getAngle(const AdvancedVertex *v) const;

    /**
     * @brief getDAngle
     * Angle between the normal vector of this and the one of 't'.
     * Return -1 if one or both the triangles are degenerate.
     * @param t
     * @return
     */
    double getDAngle(const AdvancedTriangle *t) const;

    /**
     * @brief distanceFromPoint
     *  Distance of 'p' from the plane of the triangle.
     * Return -1 if triangle is degenerate.
     * @param p
     * @return
     */
    double distanceFromPoint(const AdvancedPoint *p) const;

    /**
    * @brief pointTriangleDistance
    * Distance of 'p' from the closest point of the triangle.
    * Return -1 if triangle is degenerate.
    * If 'c' is not nullptr, its coordinates are set to the ones of the
    * closest point.
    * @param p
    * @param c
    * @return
    */
    double pointTriangleDistance(const AdvancedPoint *p, AdvancedPoint *c = nullptr) const;

    /**
      * @brief overlaps
      * Return TRUE iff one of the adjacent triangles overlaps with this one
      * @return
      */
    bool overlaps() const;
};

#define FOR_EACH_TRIANGLE_EDGE(t, e)                                                                \
    for ((e) = (t)->e1; (e) != nullptr;                                                             \
         (e) = ((e) == (t)->e3) ? (nullptr) : ((t)->nextEdge(e)))

}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_TRIANGLE_H

