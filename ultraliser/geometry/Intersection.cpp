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

#include <common/Common.h>
#include <math/Macros.h>
#include <geometry/Intersection.h>

#define X 0
#define Y 1
#define Z 2

namespace Ultraliser
{

int planeBoxOverlap(double normal[3], double vertices[3], double maxBox[3])
{
    double vMin[3], vMax[3], v;
    for (int q = X; q <= Z; q++)
    {
        v = vertices[q];
        if (normal[q] > 0.0)
        {
            vMin[q] = -maxBox[q] - v;
            vMax[q] =  maxBox[q] - v;
        }
        else
        {
            vMin[q] =  maxBox[q] - v;
            vMax[q] = -maxBox[q] - v;
        }
    }

    if (DOT(normal, vMin) > 0.0)
        return 0;

    if (DOT(normal, vMax) >= 0.0)
        return 1;

    return 0;
}

int checkTriangleBoxIntersection(double boxCenter[3],
                                 double boxHalfSize[3],
                                 double triverts[3][3])

{
    double v0[3], v1[3], v2[3];
    double min, max, p0, p1, p2, rad, fex, fey, fez;
    double normal[3], e0[3], e1[3], e2[3];

    // Move everything so that the boxcenter is in (0,0,0)
    SUB(v0, triverts[0], boxCenter);
    SUB(v1, triverts[1], boxCenter);
    SUB(v2, triverts[2], boxCenter);

    // Compute triangle edges
    SUB(e0, v1, v0);
    SUB(e1, v2, v1);
    SUB(e2, v0, v2);

    // Test the 9 tests first
    fex = std::abs(e0[X]);
    fey = std::abs(e0[Y]);
    fez = std::abs(e0[Z]);

    AXIS_TEST_X01(e0[Z], e0[Y], fez, fey);
    AXIS_TEST_Y02(e0[Z], e0[X], fez, fex);
    AXIS_TEST_Z12(e0[Y], e0[X], fey, fex);

    fex = std::abs(e1[X]);
    fey = std::abs(e1[Y]);
    fez = std::abs(e1[Z]);

    AXIS_TEST_X01(e1[Z], e1[Y], fez, fey);
    AXIS_TEST_Y02(e1[Z], e1[X], fez, fex);
    AXIS_TEST_Z0(e1[Y], e1[X], fey, fex);

    fex = std::abs(e2[X]);
    fey = std::abs(e2[Y]);
    fez = std::abs(e2[Z]);

    AXIS_TEST_X2(e2[Z], e2[Y], fez, fey);
    AXIS_TEST_Y1(e2[Z], e2[X], fez, fex);
    AXIS_TEST_Z12(e2[Y], e2[X], fey, fex);

    // Bullet 1
    // First test overlap in the {x,y,z}-directions
    // Find min, max of the triangle each direction, and test for overlap in
    // that direction -- this is equivalent to testing a minimal AABB around
    // the triangle against the AABB

    // Test in X-direction
    FIND_MIN_MAX(v0[X], v1[X], v2[X], min, max);
    if (min > boxHalfSize[X] || max < -boxHalfSize[X])
        return 0;   // No intersection

    // Test in Y-direction
    FIND_MIN_MAX(v0[Y], v1[Y], v2[Y], min, max);
    if (min > boxHalfSize[Y] || max < -boxHalfSize[Y])
        return 0;   // No intersection

    // Test in Z-direction
    FIND_MIN_MAX(v0[Z], v1[Z], v2[Z], min, max);
    if (min > boxHalfSize[Z] || max < -boxHalfSize[Z])
        return 0;   // No intersection

    // Bullet 2
    // Test if the box intersects the plane of the triangle
    // Compute plane equation of triangle: normal * x + d = 0
    CROSS(normal, e0, e1);

    if (!planeBoxOverlap(normal, v0, boxHalfSize)) return 0;
    return 1;   // Box and triangle overlaps
}

}

namespace tri_tri_isct
{
    /* function prototype */

    /* coplanar returns whether the triangles are coplanar
    *  source and target are the endpoints of the segment of
    *  intersection if it exists)
    */


#define SCALAR(dest,alpha,v) dest[0] = alpha * v[0]; \
    dest[1] = alpha * v[1]; \
    dest[2] = alpha * v[2];

#define CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) {\
    SUB(v1,p2,q1)\
    SUB(v2,p1,q1)\
    CROSS(N1,v1,v2)\
    SUB(v1,q2,q1)\
    if (DOT(v1,N1) > 0.0f) return 0;\
    SUB(v1,p2,p1)\
    SUB(v2,r1,p1)\
    CROSS(N1,v1,v2)\
    SUB(v1,r2,p1) \
    if (DOT(v1,N1) > 0.0f) return 0;\
                else return 1; }

    /* Permutation in a canonical form of T2's vertices */

#define TRI_TRI_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
    if (dp2 > 0.0f) { \
    if (dq2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2) \
                                  else if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
                                  else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2) }\
                else if (dp2 < 0.0f) { \
  if (dq2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
                                else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
                                else CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
    } else { \
  if (dq2 < 0.0f) { \
  if (dr2 >= 0.0f)  CHECK_MIN_MAX(p1,r1,q1,q2,r2,p2)\
                                                else CHECK_MIN_MAX(p1,q1,r1,p2,q2,r2)\
        } \
                                else if (dq2 > 0.0f) { \
    if (dr2 > 0.0f) CHECK_MIN_MAX(p1,r1,q1,p2,q2,r2)\
                                                else  CHECK_MIN_MAX(p1,q1,r1,q2,r2,p2)\
        } \
                                else  { \
    if (dr2 > 0.0f) CHECK_MIN_MAX(p1,q1,r1,r2,p2,q2)\
                                                else if (dr2 < 0.0f) CHECK_MIN_MAX(p1,r1,q1,r2,p2,q2)\
                                                else return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
          }}}

    /*
    *
    *  Three-dimensional Triangle-Triangle Overlap Test
    *
    */

    int tri_tri_overlap_test_3d(
        const Scalar p1[3], const Scalar q1[3], const Scalar r1[3],
        const Scalar p2[3], const Scalar q2[3], const Scalar r2[3])
    {
        Scalar dp1, dq1, dr1, dp2, dq2, dr2;
        Scalar v1[3], v2[3];
        Scalar N1[3], N2[3];

        /* Compute distance signs  of p1, q1 and r1 to the plane of
        triangle(p2,q2,r2) */

        SUB(v1, p2, r2)
            SUB(v2, q2, r2)
            CROSS(N2, v1, v2)

            SUB(v1, p1, r2)
            dp1 = DOT(v1, N2);
        SUB(v1, q1, r2)
            dq1 = DOT(v1, N2);
        SUB(v1, r1, r2)
            dr1 = DOT(v1, N2);

        if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))  return 0;

        /* Compute distance signs  of p2, q2 and r2 to the plane of
        triangle(p1,q1,r1) */

        SUB(v1, q1, p1)
            SUB(v2, r1, p1)
            CROSS(N1, v1, v2)

            SUB(v1, p2, r1)
            dp2 = DOT(v1, N1);
        SUB(v1, q2, r1)
            dq2 = DOT(v1, N1);
        SUB(v1, r2, r1)
            dr2 = DOT(v1, N1);

        if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f)) return 0;

        /* Permutation in a canonical form of T1's vertices */

        if (dp1 > 0.0f) {
            if (dq1 > 0.0f) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
            else if (dr1 > 0.0f) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
            else TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
        }
        else if (dp1 < 0.0f) {
            if (dq1 < 0.0f) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
            else if (dr1 < 0.0f) TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
            else TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
        }
        else {
            if (dq1 < 0.0f) {
                if (dr1 >= 0.0f) TRI_TRI_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
                else TRI_TRI_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
            }
            else if (dq1 > 0.0f) {
                if (dr1 > 0.0f) TRI_TRI_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
                else TRI_TRI_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
            }
            else  {
                if (dr1 > 0.0f) TRI_TRI_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
                else if (dr1 < 0.0f) TRI_TRI_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
                else return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1, N2);
            }
        }
    };

    int coplanar_tri_tri3d(
        const Scalar p1[3], const Scalar q1[3], const Scalar r1[3],
        const Scalar p2[3], const Scalar q2[3], const Scalar r2[3],
        const Scalar normal_1[3], const Scalar normal_2[3])
    {
        Scalar P1[2], Q1[2], R1[2];
        Scalar P2[2], Q2[2], R2[2];

        Scalar n_x, n_y, n_z;

        n_x = ((normal_1[0] < 0) ? -normal_1[0] : normal_1[0]);
        n_y = ((normal_1[1] < 0) ? -normal_1[1] : normal_1[1]);
        n_z = ((normal_1[2] < 0) ? -normal_1[2] : normal_1[2]);

        /* Projection of the triangles in 3D onto 2D such that the area of
        the projection is maximized. */

        if ((n_x > n_z) && (n_x >= n_y)) {
            // Project onto plane YZ

            P1[0] = q1[2]; P1[1] = q1[1];
            Q1[0] = p1[2]; Q1[1] = p1[1];
            R1[0] = r1[2]; R1[1] = r1[1];

            P2[0] = q2[2]; P2[1] = q2[1];
            Q2[0] = p2[2]; Q2[1] = p2[1];
            R2[0] = r2[2]; R2[1] = r2[1];
        }
        else if ((n_y > n_z) && (n_y >= n_x)) {
            // Project onto plane XZ

            P1[0] = q1[0]; P1[1] = q1[2];
            Q1[0] = p1[0]; Q1[1] = p1[2];
            R1[0] = r1[0]; R1[1] = r1[2];

            P2[0] = q2[0]; P2[1] = q2[2];
            Q2[0] = p2[0]; Q2[1] = p2[2];
            R2[0] = r2[0]; R2[1] = r2[2];
        }
        else {
            // Project onto plane XY

            P1[0] = p1[0]; P1[1] = p1[1];
            Q1[0] = q1[0]; Q1[1] = q1[1];
            R1[0] = r1[0]; R1[1] = r1[1];

            P2[0] = p2[0]; P2[1] = p2[1];
            Q2[0] = q2[0]; Q2[1] = q2[1];
            R2[0] = r2[0]; R2[1] = r2[1];
        }

        return tri_tri_overlap_test_2d(P1, Q1, R1, P2, Q2, R2);
    };

    /*
    *
    *  Three-dimensional Triangle-Triangle Intersection
    *
    */

    /*
    This macro is called when the triangles surely intersect
    It constructs the segment of intersection of the two triangles
    if they are not coplanar.
    */

#define CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) { \
    SUB(v1,q1,p1) \
    SUB(v2,r2,p1) \
    CROSS(N,v1,v2) \
    SUB(v,p2,p1) \
    if (DOT(v,N) > 0.0f) {\
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) <= 0.0f) { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) > 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
                } else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p2,p1) \
    SUB(v2,p2,r2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p2,v1) \
    return 1; \
                } \
                } else { \
    return 0; \
                } \
                } else { \
    SUB(v2,q2,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) < 0.0f) { \
    return 0; \
                } else { \
    SUB(v1,r1,p1) \
    CROSS(N,v1,v2) \
    if (DOT(v,N) >= 0.0f) { \
    SUB(v1,p1,p2) \
    SUB(v2,p1,r1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p1,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
                } else { \
    SUB(v1,p2,p1) \
    SUB(v2,p2,q2) \
    alpha = DOT(v1,N1) / DOT(v2,N1); \
    SCALAR(v1,alpha,v2) \
    SUB(source,p2,v1) \
    SUB(v1,p1,p2) \
    SUB(v2,p1,q1) \
    alpha = DOT(v1,N2) / DOT(v2,N2); \
    SCALAR(v1,alpha,v2) \
    SUB(target,p1,v1) \
    return 1; \
                }}}}

#define TRI_TRI_INTER_3D(p1,q1,r1,p2,q2,r2,dp2,dq2,dr2) { \
    if (dp2 > 0.0f) { \
    if (dq2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2) \
                                  else if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
                                  else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2) }\
                else if (dp2 < 0.0f) { \
  if (dq2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
                                else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
                                else CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
        } else { \
    if (dq2 < 0.0f) { \
    if (dr2 >= 0.0f)  CONSTRUCT_INTERSECTION(p1,r1,q1,q2,r2,p2)\
                                                else CONSTRUCT_INTERSECTION(p1,q1,r1,p2,q2,r2)\
                } \
                                else if (dq2 > 0.0f) { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,p2,q2,r2)\
                                                else  CONSTRUCT_INTERSECTION(p1,q1,r1,q2,r2,p2)\
        } \
                                else  { \
    if (dr2 > 0.0f) CONSTRUCT_INTERSECTION(p1,q1,r1,r2,p2,q2)\
                                                else if (dr2 < 0.0f) CONSTRUCT_INTERSECTION(p1,r1,q1,r2,p2,q2)\
                                                else { \
      coplanar = true; \
      return coplanar_tri_tri3d(p1,q1,r1,p2,q2,r2,N1,N2);\
        } \
        }} }

    /*
    The following version computes the segment of intersection of the
    two triangles if it exists.
    coplanar returns whether the triangles are coplanar
    source and target are the endpoints of the line segment of intersection
    */

    int tri_tri_intersection_test_3d(
        const Scalar p1[3], const Scalar q1[3], const Scalar r1[3],
        const Scalar p2[3], const Scalar q2[3], const Scalar r2[3],
        bool & coplanar,
        Scalar source[3], Scalar target[3])
    {
        coplanar = false;

        Scalar dp1, dq1, dr1, dp2, dq2, dr2;
        Scalar v1[3], v2[3], v[3];
        Scalar N1[3], N2[3], N[3];
        Scalar alpha;

        // Compute distance signs  of p1, q1 and r1
        // to the plane of triangle(p2,q2,r2)

        SUB(v1, p2, r2)
            SUB(v2, q2, r2)
            CROSS(N2, v1, v2)

            SUB(v1, p1, r2)
            dp1 = DOT(v1, N2);
        SUB(v1, q1, r2)
            dq1 = DOT(v1, N2);
        SUB(v1, r1, r2)
            dr1 = DOT(v1, N2);

        if (((dp1 * dq1) > 0.0f) && ((dp1 * dr1) > 0.0f))
            return 0;

        // Compute distance signs  of p2, q2 and r2
        // to the plane of triangle(p1,q1,r1)

        SUB(v1, q1, p1)
            SUB(v2, r1, p1)
            CROSS(N1, v1, v2)

            SUB(v1, p2, r1)
            dp2 = DOT(v1, N1);
        SUB(v1, q2, r1)
            dq2 = DOT(v1, N1);
        SUB(v1, r2, r1)
            dr2 = DOT(v1, N1);

        if (((dp2 * dq2) > 0.0f) && ((dp2 * dr2) > 0.0f))
            return 0;

        // Permutation in a canonical form of T1's vertices

        if (dp1 > 0.0f) {
            if (dq1 > 0.0f)
                TRI_TRI_INTER_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
            else if (dr1 > 0.0f)
            TRI_TRI_INTER_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
            else
            TRI_TRI_INTER_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
        }
        else if (dp1 < 0.0f) {
            if (dq1 < 0.0f)
                TRI_TRI_INTER_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
            else if (dr1 < 0.0f)
            TRI_TRI_INTER_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
            else
            TRI_TRI_INTER_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
        }
        else {
            if (dq1 < 0.0f) {
                if (dr1 >= 0.0f)
                    TRI_TRI_INTER_3D(q1, r1, p1, p2, r2, q2, dp2, dr2, dq2)
                else
                TRI_TRI_INTER_3D(p1, q1, r1, p2, q2, r2, dp2, dq2, dr2)
            }
            else if (dq1 > 0.0f) {
                if (dr1 > 0.0f)
                    TRI_TRI_INTER_3D(p1, q1, r1, p2, r2, q2, dp2, dr2, dq2)
                else
                TRI_TRI_INTER_3D(q1, r1, p1, p2, q2, r2, dp2, dq2, dr2)
            }
            else  {
                if (dr1 > 0.0f)
                    TRI_TRI_INTER_3D(r1, p1, q1, p2, q2, r2, dp2, dq2, dr2)
                else if (dr1 < 0.0f)
                TRI_TRI_INTER_3D(r1, p1, q1, p2, r2, q2, dp2, dr2, dq2)
                else {
                    // triangles are co-planar

                    coplanar = true;
                    return coplanar_tri_tri3d(p1, q1, r1, p2, q2, r2, N1, N2);
                }
            }
        }
    };

    /*
    *
    *  Two dimensional Triangle-Triangle Overlap Test
    *
    */

    /* some 2D macros */

#define ORIENT_2D(a, b, c)  ((a[0]-c[0])*(b[1]-c[1])-(a[1]-c[1])*(b[0]-c[0]))

#define INTERSECTION_TEST_VERTEX(P1, Q1, R1, P2, Q2, R2) {\
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f)\
    if (ORIENT_2D(R2,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(P1,P2,Q1) > 0.0f) {\
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f) return 1; \
                                else return 0;} else {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f)\
    if (ORIENT_2D(Q1,R1,P2) >= 0.0f) return 1; \
                                                else return 0;\
                                else return 0;}\
                                else \
    if (ORIENT_2D(P1,Q2,Q1) <= 0.0f)\
    if (ORIENT_2D(R2,Q2,R1) <= 0.0f)\
    if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) return 1; \
                                                else return 0;\
                                else return 0;\
                                                else return 0;\
                else\
  if (ORIENT_2D(R2,P2,R1) >= 0.0f) \
  if (ORIENT_2D(Q1,R1,R2) >= 0.0f)\
  if (ORIENT_2D(P1,P2,R1) >= 0.0f) return 1;\
                                else return 0;\
                                                else \
      if (ORIENT_2D(Q1,R1,Q2) >= 0.0f) {\
      if (ORIENT_2D(R2,R1,Q2) >= 0.0f) return 1; \
                                                else return 0; }\
                                else return 0; \
                                else  return 0; \
        };

#define INTERSECTION_TEST_EDGE(P1, Q1, R1, P2, Q2, R2) { \
    if (ORIENT_2D(R2,P2,Q1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,Q1) >= 0.0f) { \
    if (ORIENT_2D(P1,Q1,R2) >= 0.0f) return 1; \
                                                                else return 0;} else { \
        if (ORIENT_2D(Q1,R1,P2) >= 0.0f){ \
        if (ORIENT_2D(R1,P1,P2) >= 0.0f) return 1; else return 0;} \
                                                else return 0; } \
                } else {\
    if (ORIENT_2D(R2,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,P2,R1) >= 0.0f) {\
    if (ORIENT_2D(P1,R1,R2) >= 0.0f) return 1;  \
                                else {\
    if (ORIENT_2D(Q1,R1,R2) >= 0.0f) return 1; else return 0;}}\
                                                else  return 0; }\
                                else return 0; }}

    int ccw_tri_tri_intersection_2d(
        const Scalar p1[2], const Scalar q1[2], const Scalar r1[2],
        const Scalar p2[2], const Scalar q2[2], const Scalar r2[2])
    {
        if (ORIENT_2D(p2, q2, p1) >= 0.0f) {
            if (ORIENT_2D(q2, r2, p1) >= 0.0f) {
                if (ORIENT_2D(r2, p2, p1) >= 0.0f) return 1;
                else INTERSECTION_TEST_EDGE(p1, q1, r1, p2, q2, r2)
            }
            else {
                if (ORIENT_2D(r2, p2, p1) >= 0.0f)
                    INTERSECTION_TEST_EDGE(p1, q1, r1, r2, p2, q2)
                else INTERSECTION_TEST_VERTEX(p1, q1, r1, p2, q2, r2)
            }
        }
        else {
            if (ORIENT_2D(q2, r2, p1) >= 0.0f) {
                if (ORIENT_2D(r2, p2, p1) >= 0.0f)
                    INTERSECTION_TEST_EDGE(p1, q1, r1, q2, r2, p2)
                else  INTERSECTION_TEST_VERTEX(p1, q1, r1, q2, r2, p2)
            }
            else INTERSECTION_TEST_VERTEX(p1, q1, r1, r2, p2, q2)
        }
    };

    int tri_tri_overlap_test_2d(
        const Scalar p1[2], const Scalar q1[2], const Scalar r1[2],
        const Scalar p2[2], const Scalar q2[2], const Scalar r2[2])
    {
        if (ORIENT_2D(p1, q1, r1) < 0.0f)
            if (ORIENT_2D(p2, q2, r2) < 0.0f)
                return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, r2, q2);
            else
                return ccw_tri_tri_intersection_2d(p1, r1, q1, p2, q2, r2);
        else
            if (ORIENT_2D(p2, q2, r2) < 0.0f)
                return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, r2, q2);
            else
                return ccw_tri_tri_intersection_2d(p1, q1, r1, p2, q2, r2);
    };
};

