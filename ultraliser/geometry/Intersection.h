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

#ifndef ULTRALISER_GEOMETRY_INTERSECTION_H
#define ULTRALISER_GEOMETRY_INTERSECTION_H

namespace Ultraliser
{

/**
 * @brief checkTriangleBoxIntersection
 * Use separating axis theorem to test overlap between triangle and box (AABB).
 * Need to test for overlap in these directions:
 *      1) The {x,y,z}-directions (actually, since we use the AABB of the
 *         triangle we do not even need to test these).
 *      2) Normal of the triangle.
 *      3) Cross-product(edge from tri, {x,y,z}-direction).
 * @param boxcenter
 * @param boxHalfSize
 * @param triverts
 * @note The original AABB-triangle overlap test code was written by
 * Tomas Akenine-Möller.
 * @return
 */
int checkTriangleBoxIntersection(double boxCenter[3],
                                 double boxHalfSize[3],
                                 double triverts[3][3]);

}

namespace tri_tri_isct
{
    typedef float Scalar;

    /* function prototype */

    int tri_tri_overlap_test_3d(
        const Scalar p1[3], const Scalar q1[3], const Scalar r1[3],
        const Scalar p2[3], const Scalar q2[3], const Scalar r2[3]);

    int coplanar_tri_tri3d(
        const Scalar  p1[3], const Scalar  q1[3], const Scalar  r1[3],
        const Scalar  p2[3], const Scalar  q2[3], const Scalar  r2[3],
        const Scalar  N1[3], const Scalar  N2[3]);

    int tri_tri_overlap_test_2d(
        const Scalar p1[2], const Scalar q1[2], const Scalar r1[2],
        const Scalar p2[2], const Scalar q2[2], const Scalar r2[2]);

    int tri_tri_intersection_test_3d(
        const Scalar p1[3], const Scalar q1[3], const Scalar r1[3],
        const Scalar p2[3], const Scalar q2[3], const Scalar r2[3],
        bool & coplanar,
        Scalar source[3], Scalar target[3]);
}

#endif // ULTRALISER_GEOMETRY_INTERSECTION_H
