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

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_POINT_H
#define ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_POINT_H

#include <common/Common.h>

namespace Ultraliser
{

/// Orientation predicates using filtering on doubles
extern "C" double orient2d(double *, double *, double *);
extern "C" double orient3d(double *, double *, double *, double *);

/**
 * @brief orient2D
 * Orientation predicates on doubles
 * orient2D: >0 =0 <0 if (p, q, r) are CCW, aligned, CW respectively
 * @param px
 * @param py
 * @param qx
 * @param qy
 * @param rx
 * @param ry
 * @return
 */
double orient2D(const double& px, const double& py,
                const double& qx, const double& qy,
                const double& rx, const double& ry);

/**
 * @brief The Point class
 * Geometric point definition.

 * This class represents a point in the Euclidean 3D space. It can be used
 * to  represent  3D vectors originating at (0,0,0) and terminating at the
 * corresponding point. Several methods of  this  class  are  intended  to
 * manipulate  vectors  rather  than  points;  for  example, a call of the
 * method normalize is an actual normalization if the object is a  vector,
 * but  it  has  to  be intended as a projection on the unit sphere if the
 * object is intended to be a point. An object of type Point is a  triplet
 * (x,y,z)  of  coordinates  endowed with a pointer 'info' to possible additional
 * information. Each coordinate is a number of type 'coord' which, by
 * default,  is  a standard double. Operations on points include addition,
 * subtraction, cross and dot product, and many others. This class  implements
 * several useful operations using vector arithmethic. For example,
 * the simple piece of code "A = B*C;" assignes to A the value of the  dot
 * product of B and C.
 * Nearly zero or nearly flat angles are automatically snapped to
 * exactly zero and exactly flat angles if the difference is smaller
 * than the global variable _acos_tolerance. This is the very basic application
 * of our version of the epsilon geometry for robust computation.
 */
class AdvancedPoint
{
public :

    /**
     * @brief x
     * Coordinates
     */
    double x;

    /**
     * @brief y
     * Coordinates
     */
    double y;

    /**
     * @brief z
     * Coordinates
     */
    double z;

    /**
     * @brief info
     */
    void *info;

    /**
     * @brief Point
     * Creates a new point with coordinates (0,0,0).
     */
    AdvancedPoint()
    {
        x = y = z = 0;
        info = nullptr;
    }

    /**
     * @brief Point
     * Creates a new point with the same coordinates as 's'.
     * The info field is not copied.
     * @param s
     */
    AdvancedPoint(const AdvancedPoint *s)
    {
        x = s->x;
        y = s->y;
        z = s->z;
        info = nullptr;
    }

    /**
     * @brief Point
     * Creates a new point with the same coordinates as 's'.
     * The info field is not copied.
     * @param s
     */
    AdvancedPoint(const AdvancedPoint& s)
    {
        x = s.x;
        y = s.y;
        z = s.z;
        info = nullptr;
    }

    /**
     * @brief Point
     * Creates a new point with coordinates (a, b, c).
     * @param a
     * @param b
     * @param c
     */
    AdvancedPoint(const double& a, const double& b, const double& c)
    {
        x = a;
        y = b;
        z = c;
        info = nullptr;
    }

    /**
     * @brief isPoint
     * Do not remove this.
     * It makes the compiler produce a vtable for this object.
     * @return
     */
    bool isPoint() const { return true; }

    /**
     * @brief setValue
     * Set the coordinates to (a, b, c).
     * @param a
     * @param b
     * @param c
     */
    void setValue(const double& a, const double& b, const double& c)
    {
        x = a;
        y = b;
        z = c;
    }

    /**
     * @brief setValue
     * Set the coordinates as those of 'p'
     * @param p
     */
    void setValue(const AdvancedPoint& p)
    {
        x = p.x;
        y = p.y;
        z = p.z;
    }

    /**
     * @brief setValue
     * Set the coordinates as those of '*p'.
     * @param p
     */
    void setValue(const AdvancedPoint *p)
    {
        x = p->x;
        y = p->y;
        z = p->z;
    }

    /**
     * @brief operator -
     * Returns the vector difference
     * @param p
     * @return
     */
    AdvancedPoint operator-(const AdvancedPoint& p) const
    {
        return AdvancedPoint(x - p.x, y - p.y, z - p.z);
    }

    /**
     * @brief operator +
     * Returns the vector sum
     * @param p
     * @return
     */
    AdvancedPoint operator+(const AdvancedPoint& p) const
    {
        return AdvancedPoint(x + p.x, y + p.y, z + p.z);
    }

    /**
     * @brief operator +=
     * Sums another point.
     * @param p
     */
    void operator+=(const AdvancedPoint& p)
    {
        x += p.x;
        y += p.y;
        z += p.z;
    }

    /**
     * @brief operator -=
     * Subtracts another point
     * @param p
     */
    void operator-=(const AdvancedPoint& p)
    {
        x -= p.x;
        y -= p.y;
        z -= p.z;
    }

    /**
     * @brief operator &
     * Returns the Cross Product
     * @param p
     * @return
     */
    AdvancedPoint operator&(const AdvancedPoint& p) const
    {
        return AdvancedPoint(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
    }

    /**
     * @brief operator *
     * Returns the Dot Product
     * @param p
     * @return
     */
    double operator*(const AdvancedPoint& p) const
    {
        return (x * p.x + y * p.y + z * p.z);
    }

    /**
     * @brief operator *
     * Returns the product with a scalar
     * @param d
     * @return
     */
    AdvancedPoint  operator*(const double& d) const
    {
        return AdvancedPoint(x * d, y * d, z * d);
    }

    /**
     * @brief operator *=
     * Multiplies by a scalar
     * @param m
     */
    void operator*=(const double& m)
    {
        x *= m;
        y *= m;
        z *= m;
    }

    /**
     * @brief operator /=
     * Divides by a scalar
     * @param m
     */
    void operator/=(const double& m)
    {
        x /= m;
        y /= m;
        z /= m;
    }

    /**
     * @brief operator /
     * Returns the vector divided by the scalar
     * @param d
     * @return
     */
    AdvancedPoint operator/(const double& d) const
    {
        return AdvancedPoint(x / d, y / d, z / d);
    }

    /**
     * @brief operator ==
     * TRUE iff coordinates are equal
     * @param p
     * @return
     */
    bool operator==(const AdvancedPoint& p) const
    {
        return (x == p.x && y == p.y && z == p.z);
    }

    /**
     * @brief operator !=
     * FALSE iff coordinates are equal
     * @param p
     * @return
     */
    bool operator!=(const AdvancedPoint& p) const
    {
        return (x != p.x || y != p.y || z != p.z);
    }

    /**
     * @brief operator <
     * TRUE iff this is lexycographically smaller than s
     * @param s
     * @return
     */
    bool operator<(const AdvancedPoint& s) const;

    /**
     * @brief at
     * Returns the i'th coordinate
     * @param i
     * @return
     */
    inline double& at(unsigned char i)
    {
        return (i == 0) ? (x) : ((i == 1) ? (y) : (z));
    }

    /**
     * @brief operator []
     * Returns the i'th coordinate
     * @param i
     * @return
     */
    inline double& operator[](unsigned char i)
    {
        return (i == 0) ? (x) : ((i == 1) ? (y) : (z));
    }

    /**
     * @brief inverse
     * Returns the inverse vector
     * @return
     */
    AdvancedPoint inverse() const { return AdvancedPoint(-x, -y, -z); }

    /**
     * @brief invert
     * Inverts the vector
     */
    void invert()
    {
        x = -x;
        y = -y;
        z = -z;
    }

    /**
     * @brief isNull
     * TRUE if vector is (0,0,0)
     * @return
     */
    bool isNull() const { return (x == 0 && y == 0 && z == 0); }

    /**
     * @brief squaredLength
     * Squared distance from origin
     * @return
     */
    double squaredLength() const { return (x * x + y * y + z * z); }

    /**
     * @brief squaredDistance
     * Squared distance from '*b'
     * @param b
     * @return
     */
    double squaredDistance(const AdvancedPoint *b) const
    {
        return (((*(this)) - (*b)).squaredLength());
    }

    /**
     * @brief linearSystem
     * Returns the solution of the linear system Ax = d, where A is a 3x3 matrix
     * whose rows are row1, row2 and row3, d = this
     * @param row1
     * @param row2
     * @param row3
     * @return
     */
    AdvancedPoint linearSystem(const AdvancedPoint& row1, const AdvancedPoint& row2, const AdvancedPoint& row3);

    /**
     * @brief project
     * Projects the vector on the plane with normal 'n' passing through the origin.
     * @param n
     */
    void project(const AdvancedPoint *n);

    /**
     * @brief projection
     * Returns the projection of the point on the straight line though 'a' and 'b'.
     * @param a
     * @param b
     * @return
     */
    AdvancedPoint projection(const AdvancedPoint *a, const AdvancedPoint *b) const;

    /**
     * @brief printPoint
     * Prints the coordinates of the point to a file handler. stdout is the default.
     * @param fp
     */
    void printPointToFile(FILE *filePointer = stdout) const
    {
        fprintf(filePointer, "%f %f %f,\n", float(x), float(y), float(z));
    }

    /**
     * @brief exactOrientation
     * Exact orientation test.
     * Return value is positive iff the tetrahedron (this,a, b, c) has a positive volume;
     * It is negative iff the tetrahedron (this,a, b, c) has a negative volume;
     * It is zero iff the tetrahedron (this,a, b, c) has a zero volume.
     * @param a
     * @param b
     * @param c
     * @return
     */
    double exactOrientation(const AdvancedPoint *a, const AdvancedPoint *b, const AdvancedPoint *c) const;

    /**
     * @brief side3D
     * @param p1
     * @param p2
     * @param p3
     * @return
     */
    double side3D(const AdvancedPoint *p1, const AdvancedPoint *p2, const AdvancedPoint *p3) const
    {
        return exactOrientation(p1, p2, p3);
    }

    /**
     * @brief exactMisalignment
     * Exact misalignment test. Returns TRUE iff points are not aligned.
     * @param a
     * @param b
     * @return
     */
    bool exactMisalignment(const AdvancedPoint *a, const AdvancedPoint *b) const;

    /**
     * @brief notAligned
     * @param a
     * @param b
     * @return
     */
    bool notAligned(const AdvancedPoint *a, const AdvancedPoint *b) const
    {
        return exactMisalignment(a,b);
    }

    /**
     * @brief exactSameSideOnPlane
     * Exact planar side test. Returns TRUE iff 'this', Q, A and B are coplanar
     * and 'this' and Q are (properly) on the same side of A-B.
     * Warning! Coplanarity is not checked, result is undetermined if
     * 'this', Q, A and B are not coplanar.
     * @param Q
     * @param A
     * @param B
     * @return
     */
    bool exactSameSideOnPlane(const AdvancedPoint *Q, const AdvancedPoint *A, const AdvancedPoint *B) const;

    /**
     * @brief lineLineIntersection
     * Itersection point between lines p-q and r-s.
     * Return INFINITE_POINT if lineas are either parallel or degenerate.
     * @param p
     * @param q
     * @param r
     * @param s
     * @return
     */
    static AdvancedPoint lineLineIntersection(const AdvancedPoint& p, const AdvancedPoint& q,
                                              const AdvancedPoint& r, const AdvancedPoint& s);

    /**
     * @brief linePlaneIntersection
     * Itersection point between line p-q and plane r-s-t.
     * Return INFINITE_POINT for parallel/degenerate args.
     * @param p
     * @param q
     * @param r
     * @param s
     * @param t
     * @return
     */
    static AdvancedPoint linePlaneIntersection(const AdvancedPoint& p, const AdvancedPoint& q,
                                               const AdvancedPoint& r, const AdvancedPoint& s,
                                               const AdvancedPoint& t);

    /**
     * @brief squaredTriangleArea3D
     *  Squared area of the triangle p-q-r.
     * @param p
     * @param q
     * @param r
     * @return
     */
    static double squaredTriangleArea3D(const AdvancedPoint& p,
                                        const AdvancedPoint& q,
                                        const AdvancedPoint& r);

    /**
     * @brief pointInInnerSegment
     * True if 'p' is a point of the segment v1-v2 (endpoints excluded)
     * @param p
     * @param v1
     * @param v2
     * @return
     */
    static bool pointInInnerSegment(const AdvancedPoint *p,
                                    const AdvancedPoint *v1,
                                    const AdvancedPoint *v2);

    /**
     * @brief pointInSegment
     * True if 'p' is a point of the segment v1-v2 (endpoints included).
     *
     * @param p
     * @param v1
     * @param v2
     * @return
     * Returns true if 'p' is a point of the segment v1-v2 (endpoints included).
     */
    static bool pointInSegment(const AdvancedPoint *p,
                               const AdvancedPoint *v1,
                               const AdvancedPoint *v2);

    /**
     * @brief pointInInnerTriangle
     * True if the coplanar point 'p' is in the inner area of v1-v2-v3.
     * Undetermined if points are not coplanar.
     * @param p
     * @param v1
     * @param v2
     * @param v3
     * @return
     */
    static bool pointInInnerTriangle(const AdvancedPoint *p,
                                     const AdvancedPoint *v1,
                                     const AdvancedPoint *v2,
                                     const AdvancedPoint *v3);

    /**
     * @brief pointInTriangle
     * True if the coplanar point 'p' is either in the inner area of v1-v2-v3 or on its border.
     * Undetermined if points are not coplanar.
     *
     * @param p
     * @param v1
     * @param v2
     * @param v3
     * @return
     * Returns true if the coplanar point 'p' is either in the inner area of 't' or on its border.
     * Undetermined if p and t are not coplanar.
     */
    static bool pointInTriangle(const AdvancedPoint *p,
                                const AdvancedPoint *v1,
                                const AdvancedPoint *v2,
                                const AdvancedPoint *v3);

    /**
     * @brief segmentsIntersect
     * True if (p1 - p2) properly intersects (sp1 - sp2) at any point (endpoints included).
     * Collinear overlapping segments are not considered to be properly intersecting.
     * @param p1
     * @param p2
     * @param sp1
     * @param sp2
     * @return
     */
    static bool segmentsIntersect(const AdvancedPoint *p1,
                                  const AdvancedPoint *p2,
                                  const AdvancedPoint *sp1,
                                  const AdvancedPoint *sp2);

    /**
     * @brief innerSegmentsCross
     * True if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
     * Collinear overlapping segments are not considered to be properly intersecting.
     *
     * @param p1
     * @param p2
     * @param sp1
     * @param sp2
     * @return
     * Returns true if the interior of (p1-p2) properly intersects the interior of (sp1-sp2).
     */
    static bool innerSegmentsCross(const AdvancedPoint& p1,
                                   const AdvancedPoint& p2,
                                   const AdvancedPoint& sp1,
                                   const AdvancedPoint& sp2);

    /**
     * @brief segmentIntersectsTriangle
     * True if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
     * @param s1
     * @param s2
     * @param v1
     * @param v2
     * @param v3
     * @return
     */
    static bool segmentIntersectsTriangle(const AdvancedPoint *s1,
                                          const AdvancedPoint *s2,
                                          const AdvancedPoint *v1,
                                          const AdvancedPoint *v2,
                                          const AdvancedPoint *v3);

    /**
     * @brief segmentIntersectsTriangle
     * True if segment (s1-s2) intersects the triangle v1-v2-v3 (border included).
     * Accelerated version - relative orientations are passed as parameters.
     * @param s1
     * @param s2
     * @param v1
     * @param v2
     * @param v3
     * @param o1
     * @param o2
     * @return
     */
    static bool segmentIntersectsTriangle(const AdvancedPoint *s1,
                                          const AdvancedPoint *s2,
                                          const AdvancedPoint *v1,
                                          const AdvancedPoint *v2,
                                          const AdvancedPoint *v3,
                                          const double& o1,
                                          const double& o2);

    /**
     * @brief length
     * Distance from origin
     * @return
     */
    double length() const
    {
        return sqrt(double((x * x + y * y + z * z)));
    }

    /**
     * @brief normalize
     * Divides the vector by its length.
     * If isNull() the application exits with an error.
     */
    void normalize();

    /**
     * @brief rotate
     * Rotates the vector around 'axis' by 'ang' radians ccw.
     * @param axis
     * @param ang
     */
    void rotate(const AdvancedPoint& axis, const double& ang);

    /**
     * @brief distance
     *  Distance from 'b'
     * @param b
     * @return
     */
    double distance(const AdvancedPoint& b) const
    {
        return (((*(this)) - (b)).length());
    }

    /**
     * @brief distance
     * Distance from '*b'
     * @param b
     * @return
     */
    double distance(const AdvancedPoint *b) const
    {
        return (((*(this)) - (*b)).length());
    }

    /**
     * @brief distanceFromLine
     * Distance from straight line through 'a' and 'b'
     * @param a
     * @param b
     * @return
     */
    double distanceFromLine(const AdvancedPoint *a, const AdvancedPoint *b) const;

    /**
     * @brief distanceFromLine
     * Distance from straight line through 'a' and 'b'. *cc is set to the closest line point.
     * @param a
     * @param b
     * @param cc
     * @return
     */
    double distanceFromLine(const AdvancedPoint *a, const AdvancedPoint *b, AdvancedPoint *cc) const;

    /**
     * @brief distanceFromEdge
     * Distance from segment a-b
     * @param a
     * @param b
     * @return
     */
    double distanceFromEdge(const AdvancedPoint *a, const AdvancedPoint *b) const;

    /**
     * @brief distanceFromEdge
     * Distance from segment a-b. *cc is set to the closest edge point.
     * @param a
     * @param b
     * @param cc
     * @return
     */
    double distanceFromEdge(const AdvancedPoint *a, const AdvancedPoint *b, AdvancedPoint *cc) const;

    /**
     * @brief distanceLineLine
     * Distance between the straight lines through (this) - l1_p2 and l2_p1 - l2_p2.
     * @param l1_p2
     * @param l2_p1
     * @param l2_p2
     * @return
     */
    double distanceLineLine(const AdvancedPoint *l1_p2,
                            const AdvancedPoint *l2_p1,
                            const AdvancedPoint *l2_p2) const;

    /**
     * @brief getAngle
     * Angle between this vector and 'v' in radians.
     * @param v
     * @return
     */
    double getAngle(const AdvancedPoint& v) const;

    /**
     * @brief getAngle
     * Angle defined by <a, *this, b> in radians.
     * @param a
     * @param b
     * @return
     */
    double getAngle(const AdvancedPoint& a, const AdvancedPoint& b) const
    {
        return (a - (*this)).getAngle(b - (*this));
    }

    /**
     * @brief getAngle
     * Angle defined by <*a, *this, *b> in radians.
     * @param a
     * @param b
     * @return
     */
    double getAngle(const AdvancedPoint *a, const AdvancedPoint *b) const
    {
        return ((*a) - (*this)).getAngle((*b) - (*this));
    }

    /**
     * @brief closestPoints
     * Line-line closest point computation.
     * Computes the closest points of the line passing through this and this2,
     * and the line passing through p1 and p2. The computed points are used to
     * initialize the  coordinates  of  cpOnThis  and  cpOnOther.  The  method
     * returns 0 if the lines are parallel, 1 otherwise.
     * @param this2
     * @param p1
     * @param p2
     * @param cpOnThis
     * @param cpOnOther
     * @return
     */
    int closestPoints(const AdvancedPoint *this2,
                      const AdvancedPoint *p1, const AdvancedPoint *p2,
                      AdvancedPoint *cpOnThis, AdvancedPoint *cpOnOther) const;
};

/**
 * @brief xyzCompare
 * Lexycographic comparison to be used with jqsort() or abstractHeap.
 * @param p1
 * @param p2
 * @return
 */
int xyzCompare(const void *p1, const void *p2);

/**
 * @brief INFINITE_POINT
 * Static point with DBL_MAX coordinates.
 */
extern const AdvancedPoint INFINITE_POINT;

/// Checks whether a point is INFINITE_POINT.
#define IS_FINITE_POINT(p) ((p).x < DBL_MAX)

}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_POINT_H

