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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_TRIANGLE_OPERATIONS_H
#define ULTRALISER_DATA_MESH_SIMPLE_TRIANGLE_OPERATIONS_H

#include <common/Common.h>
#include <data/common/CommonData.h>

namespace Ultraliser
{

bool isPointInSphere(const Vector3f& p, const Vector3f& center, const float& radius);

/**
 * @brief isTriangleInSphere
 * @param p0
 * @param p1
 * @param p2
 * @param center
 * @param radius
 * @return
 */
bool isTriangleInSphere(Vector3f p0, Vector3f p1, Vector3f p2,
                        const Vector3f& center, const float& radius);

/**
 * @brief computeTriangleSurfaceArea
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleSurfaceArea(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleSignedVolume
 * @param p1
 * @param p2
 * @param p3
 * @return
 */
float computeTriangleSignedVolume(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleAspectRatio
 * Computes the aspect ratio of the triangle.
 * The aspect ratio of a triangle is the ratio of its circumradius to the
 * diameter of its incircle. In other words, the aspect ratio of a triangle is
 * the ratio of the longest edge to shortest edge (so equilateral triangle has
 * aspect ratio 1).
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleAspectRatio(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleAspectFrobenius
 * The aspect Frobenius is the sum of the edge lengths squared divided by the
 * area and normalized so that a unit equilateral triangle has a value of 1.
 *      Dimension:                           1
 *      Acceptable Range:                   (1, 1.3]
 *      Normal Range:                       [1, FLOAT_MAXIMUM]
 *      Full Range:                         [1, FLOAT_MAXIMUM]
 *      Value for equilateral unit triangle: 1
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleAspectFrobenius(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleInRadius
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleInRadius(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleCircumradius
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleCircumradius(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleRadiusRatio
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleRadiusRatio(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleScaledJacobian
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleScaledJacobian(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleDistrotion
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleDistrotion(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleConditionNumber
 * The condition number of the weighted Jacobian matrix.
 * In theory the condition number is invari- ant to which node it is computed
 * at, but floating point truncation error can contribute to differences
 * between values computed for each node.
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleConditionNumber(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleEdgeRatio
 * Computes the triangle edge ratio. This ratio is defined as the ratio between
 * the longest and shortest edges of the triangle.
 * Dimension: 1
 * Acceptable Range: Between 1 and 1.3
 * Normal Range: 1 to infinity
 * Full Range: 1 to infinity
 * Value for equilateral unit triangle: 1
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleEdgeRatio(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleAngles
 * @param p0
 * @param p1
 * @param p2
 * @param angle1
 * @param angle2
 * @param angle3
 * @return
 */
void computeTriangleAngles(Vector3f p0, Vector3f p1, Vector3f p2,
                            float& angle1, float &angle2, float &angle3);

/**
 * @brief computeMinTriangleAngle
 * Computes the minimum angle of the given triangle. The result should be in
 * degrees. Note that if any of the edge vectors has zero length, the return
 * value will be zero.
 * Acceptable range is between 30 and 60 degrees.
 * Normal range is between 0 and 60 degrees.
 * Full range is between 0 and 360 degrees.
 * Value for equilateral unit triangle is 60 degrees.
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleMinAngle(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeMaxTriangleAngle
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleMaxAngle(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleRelativeSizeSquared
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleRelativeSizeSquared(Vector3f p0, Vector3f p1, Vector3f p2,
                                         const float averageArea);

/**
 * @brief computeTriangleShape
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleShape(Vector3f p0, Vector3f p1, Vector3f p2);

/**
 * @brief computeTriangleShapeAndSize
 * @param p0
 * @param p1
 * @param p2
 * @return
 */
float computeTriangleShapeAndSize(Vector3f p0, Vector3f p1, Vector3f p2,
                                  const float averageArea);
}
#endif // ULTRALISER_DATA_MESH_SIMPLE_TRIANGLE_OPERATIONS_H
