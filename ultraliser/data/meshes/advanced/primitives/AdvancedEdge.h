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

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_EDGE_H
#define ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_EDGE_H

#include "AdvancedPoint.h"
#include "AdvancedVertex.h"

namespace Ultraliser
{

/**
 * @brief The Edge class
 * This  class  represents an edge of a triangulation. An Edge is the main
 * part of the Basic_TMesh data structure. Each edge has an  orientation
 * (i.e. from v1 to v2) which forces the order in which incident triangles
 * (t1 and t2) are stored in the class. When looking the edge so  that  it
 * points  "upwards", if the normal of t1 points towards the observer then
 * t1 must be on the left of the  edge.  The  field  mask  is  useful  for
 * assigning up to 256 different states to the edge.
 */
class AdvancedEdge
{
public:

    /**
     * @brief v1
     * Edge vertex 1.
     */
    AdvancedVertex* v1;

    /**
     * @brief v2
     * Edge vertex 2.
     */
    AdvancedVertex* v2;

    /**
     * @brief t1
     * Incident triangle 1.
     */
    class AdvancedTriangle* t1;

    /**
     * @brief t2
     * Incident triangle 2.
     */
    class AdvancedTriangle* t2;

    /**
     * @brief mask
     * A bit-mask for marking purposes
     */
    unsigned char mask;

    //!< Further information
    void *info;

    /**
     * @brief Edge
     */
    AdvancedEdge();

    /**
     * @brief Edge
     * @param s
     * @param d
     */
    AdvancedEdge(AdvancedVertex *s, AdvancedVertex *d);
    ~AdvancedEdge();

    /**
     * @brief isBaseType
     * Returns true only if object is a basic Edge. All the reimplementations
     * must return false.
     * @return
     */
    bool isBaseType() const { return true; }

    /**
     * @brief isLinked
     * TRUE iff edge is properly linked to a TMesh.
     * @return
     */
    bool isLinked()	const 	{ return (v1 != nullptr);}

    /**
     * @brief hasVertex
     * TRUE iff 'v' is an end-point of the edge.
     * @param v
     * @return
     */
    bool hasVertex(const AdvancedVertex *v) const { return (v1==v || v2==v);}

    /**
     * @brief hasTriangle
     * TRUE iff 't' is incident to the edge.
     * @param t
     * @return
     */
    bool hasTriangle(const AdvancedTriangle *t) const { return (t1==t || t2==t);}

    /**
     * @brief hasVertices
     * TRUE if both 'va' and 'vb' are vertices of the edge.
     * @param va
     * @param vb
     * @return
     */
    bool hasVertices(const AdvancedVertex *va, const AdvancedVertex *vb) const
    {
        return ((v1 == va && v2 == vb) || (v2 == va && v1 == vb));
    }

    /**
     * @brief squaredLength
     * Squared length of the edge.
     * @return
     */
    double squaredLength() const 	{ return v1->squaredDistance(v2);}

    /**
     * @brief toVector
     * Convert to vector v2 - v1.
     * @return
     */
    AdvancedPoint toVector() const
    {
        return (*v2) - (*v1);
    }

    /**
     * @brief getMidPoint
     * Return the edge's mid-point.
     * @return
     */
    AdvancedPoint getMidPoint()	const
    {
        return ((*v1) + (*v2)) / 2.0;
    }

    /**
     * @brief invert
     * Invert the edge's orientation.
     */
    void invert()
    {
        p_swap((void **)(&v1), (void **)(&v2));
        p_swap((void **)(&t1), (void **)(&t2));
    }

    /**
     * @brief leftTriangle
     * Triangle on the left of the edge when looking from 'v'. nullptr if 'v'
     * is not a vertex of the edge.
     * @param v
     * @return
     */
    AdvancedTriangle *leftTriangle(const AdvancedVertex *v) const
    {
        return ((v1 == v)?(t1):((v2 == v)?(t2):(nullptr)));
    }

    /**
     * @brief rightTriangle
     * Triangle on the right of the edge when looking from 'v'. nullptr if 'v'
     * is not a vertex of the edge.
     * @param v
     * @return
     */
    AdvancedTriangle *rightTriangle(const AdvancedVertex *v) const
    {
        return ((v1 == v)?(t2):((v2 == v)?(t1):(nullptr)));
    }

    /**
     * @brief oppositeVertex
     * Vertex opposite to 'v'. nullptr if 'v' is not a vertex of the edge.
     * @param v
     * @return
     */
    AdvancedVertex *oppositeVertex(const AdvancedVertex *v) const
    {
        return ((v1 == v) ? (v2) : ((v2 == v) ? (v1) : (nullptr)));
    }


    /**
     * @brief oppositeTriangle
     * Incident triangle opposite to 't'. nullptr if 't' is not incident
     * to the edge.
     * @param t
     * @return
     */
    AdvancedTriangle *oppositeTriangle(const AdvancedTriangle *t) const
    {
        return ((t1 == t) ? (t2) : ((t2 == t) ? (t1) : (nullptr)));
    }

    /**
     * @brief replaceVertex
     * Replace vertex 'a' with vertex 'b' in the edge and return TRUE.
     * If 'a' is not a vertex of the edge return FALSE.
     * @param a
     * @param b
     * @return
     */
    bool replaceVertex(const AdvancedVertex *a, AdvancedVertex *b)
    {
        if (v1 == a)
            v1 = b;
        else if (v2 == a)
            v2 = b;
        else
            return 0;
        return 1;
    }

    /**
     * @brief replaceTriangle
     * Replace incident triangle 'a' with 'b' and return TRUE. If 'a' is not
     * incident to the edge return FALSE.
     * @param a
     * @param b
     * @return
     */
    bool replaceTriangle(const AdvancedTriangle *a, AdvancedTriangle *b)
    {
        if (t1 == a)
            t1 = b;
        else if (t2 == a)
            t2=b;
        else
            return 0;
        return 1;
    }

    /**
     * @brief commonVertex
     * Vertex shared with edge 'b'. nullptr if this and 'b' do not share any
     * vertex.
     * @param b
     * @return
     */
    AdvancedVertex *commonVertex(const AdvancedEdge *b) const
    {
        return ((v1 == b->v1 || v1 == b->v2) ?
                    (v1) : ((v2 == b->v1 || v2 == b->v2) ? (v2) : (nullptr)));
    }

    /**
     * @brief isOnBoundary
     * TRUE if and only if the edge is on the boundary (i.e., one of the two
     * incident triangles is simply a nullptr).
     * @return
     */
    bool isOnBoundary() const { return (t1 == nullptr || t2 == nullptr);}

    /**
     * @brief isIsolated
     * TRUE iff edge is isolated (i.e., both the two incident triangles
     * are nullptr).
     * @return
     */
    bool isIsolated() const { return (t1 == nullptr && t2 == nullptr);}

    /**
     * @brief isDegenerate
     * TRUE iff the two endpoints coincide exactly
     * @return
     */
    bool isDegenerate() const { return ((*v1)==(*v2));}

    /**
     * @brief getBoundaryTriangle
     * If the edge is on boundary return its only incident triangle,
     * nullptr otherwise.
     * @return
     */
    AdvancedTriangle *getBoundaryTriangle() const
    {
        return (t2 == nullptr)?(t1):((t1 == nullptr)?(t2):(nullptr));
    }

    /**
     * @brief printEdge
     * Print the coordinates of the end-ponts to the file handler pointed to
     * by 'f' (stdout by default).
     * @param f
     */
    void printEdge(FILE *f =stdout) const {v1->printPointToFile(f); v2->printPointToFile(f);}

    /**
     * @brief swap
     * Combinatorial edge-swap.
     * Vertices of the edge are replaced with vertices of the two incident
     * triangles which are opposite to this edge.
     * Connectivity information is updated properly.
     * If the edge is on boundary or if the edge after the swap already exists
     * return FALSE and do not change anything.
     * Return TRUE on success.
     * If 'fast' is set, no topological check is performed.
     * @param fast
     * @return
     */
    bool swap(const bool fast=0);

    /**
     * @brief collapseOnV1
     * Edge collapse.
     * This method collapses the edge and  updates  the  connectivity  of  the
     * neighboring  elements consistently. The edge will be contracted into its
     * first vertex v1 (or v2, for collapseOnV2()).
     * This method returns v1 (or v2) on success, nullptr otherwise.
     * Failure occurs when the collapse would produce an invalid connectivity graph.
     * Caution! If the collapse succeeds the edge, its incident triangles and
     * the second vertex are unlinked, but they are still present in the lists
     * of the TMesh.
     * The calling function is responsible of removing them from the lists using
     * the method removeUnlinkedElements().
     * @return
     */
    AdvancedVertex *collapseOnV1();

    /**
     * @brief collapseOnV2
     * * Edge collapse.
     * This method collapses the edge and  updates  the  connectivity  of  the
     * neighboring  elements consistently. The edge will be contracted into its
     * first vertex v1 (or v2, for collapseOnV2()).
     * This method returns v1 (or v2) on success, nullptr otherwise.
     * Failure occurs when the collapse would produce an invalid connectivity graph.
     * Caution! If the collapse succeeds the edge, its incident triangles and
     * the second vertex are unlinked, but they are still present in the lists
     * of the TMesh.
     * The calling function is responsible of removing them from the lists using
     * the method removeUnlinkedElements().
     * @return
     */
    AdvancedVertex *collapseOnV2();

    /**
     * @brief collapse
     * This method collapses the edge and  updates  the  connectivity  of  the
     * neighboring  elements consistently. The edge will be transformed into a
     * vertex with the coordinates of 'p'.
     * This method returns TRUE on success, FALSE otherwise.
     * Failure occurs when the collapse would produce an invalid connectivity graph.
     * Caution! If the collapse succeeds the edge, its incident triangles
     * and the second vertex are unlinked, but they are still present in the
     * lists of the Basic.
     * The calling function is responsible of removing them from the lists using
     * the method removeUnlinkedElements().
     * @param p
     * @return
     */
    bool collapse(const AdvancedPoint& p);

    /**
     * @brief collapse
     * Edge collapse.
     * @return
     * This method collapses the edge and  updates  the  connectivity  of  the
     * neighboring  elements consistently. The edge will be transformed into a
     * vertex placed at the edge's mid-point.
     * This method returns TRUE on success, FALSE otherwise.
     * Failure occurs when the collapse would produce an invalid connectivity
     * graph.
     * Caution! If the collapse succeeds the edge, its incident triangles
     * and the second vertex are unlinked, but they are still present in the
     * lists of the TMesh.
     * The calling function is responsible of removing them from the lists using
     * the method removeUnlinkedElements().
     */
    bool collapse();

    /**
     * @brief merge
     * Merge with another boundary edge.
     * @param e
     * @return
     * If both this and 'e' are boundary edges, the edge 'e' is identified with
     * this one, and the connectivity of the neighboring elements is updated
     * consistently.
     * This method returns TRUE on success, FALSE otherwise.
     * Failure occurs when the merge would produce an invalid
     * connectivity graph (i.e., non orientable).
     * Caution! If the merge succeeds the edge 'e' and its two  end-points
     * are  unlinked,  but  they are still present in the lists of the
     * TMesh. It's responsibility of the calling  function  to  remove
     * them from the lists using the method removeUnlinkedElements().
     */
    bool merge(AdvancedEdge *e);

    /**
     * @brief stitch
     * Stitching primitive.
     * If there is a copy of this edge incident to one of the end-points,
     * identify it with this edge, and update the connectivity properly.
     * @return
     * This method returns TRUE on success, FALSE otherwise.
     * Caution! If the stitch succeeds, the duplicated edge
     * is unlinked, but it is still present in the lists of the
     * TMesh. It's responsibility of the calling function to remove
     * it from the lists using the method removeEdges().
     */
    bool stitch();

    /**
     * @brief overlaps
     * Returns TRUE if edge is not on boundary and its two incident
     * triangles overlap.
     * @return
     */
    bool overlaps() const;

    /**
     * @brief intersects
     * True if this edge itersects 't' other than on common subsimplexes
     * @param t
     * @return
     */
    bool intersects(const AdvancedTriangle *t) const;

    /**
     * @brief getConvexity
     * Returns a positive value if dihedral angle is less than flat (convex),
     * 0 if flat, negative if more than flat (concave). Returns DBL_MAX if
     * edge is on boundary.
     * @return
     */
    double getConvexity() const;

    /// NOTE: FUNCTIONS BELOW THIS LINE MAY RETURN APPROXIMATED/NOT ROBUST
    /// RESULTS EVEN WHEN USING RATIONALS

    /**
     * @brief length
     * Euclidean length of the edge.
     * @return
     */
    double length() const
    {
        return v1->distance(v2);
    }

    /**
     * @brief toUnitVector
     * Convert to normalized vector (v2-v1)/|v2-v1|.
     * @return
     */
    AdvancedPoint toUnitVector() const;

    /**
     * @brief getNormal
     * Return the normal at the edge as the average of the normals of the two
     * incident triangles.
     * @return
     * A null (0,0,0) vector is returned if the edge is on boundary.
     */
    AdvancedPoint getNormal() const;

    /**
     * @brief delaunayMinAngle
     * Return the minimum among the six angles of the two incident triangles.
     * @return
     * 2PI if on boundary.
     */
    double delaunayMinAngle() const;

    /**
     * @brief dihedralAngle
     * Dihedral angle at the edge.
     * @return
     * Dihedral angle at the edge.
     */
    double dihedralAngle() const;

    /**
     * @brief curvature
     * Angle between the normals of the two incident triangles.
     * Angle between the normals of the two incident triangles. If
     * the edge is on boundary or one or both the incident triangles are
     * degenerate, return -1.
     * @return
     */
    double curvature() const;
};

/**
 * @brief edgeCompare
 * Edge comparison based on length to be used with jqsort() or abstractHeap.
 * @param a
 * @param b
 * @return
 */
int edgeCompare(const void *a, const void *b);

/**
 * @brief lexEdgeCompare
 * Lexycographic edge comparison to be used with jqsort() or abstractHeap.
 * @param a
 * @param b
 * @return
 */
int lexEdgeCompare(const void *a, const void *b);

/**
 * @brief vtxEdgeCompare
 * Vertex-based edge comparison for qsort.
 * Duplicated edges are contiguous in this sorting.
 * @param a
 * @param b
 * @return
 */
int vtxEdgeCompare(const void *a, const void *b);

}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_EDGE_H

