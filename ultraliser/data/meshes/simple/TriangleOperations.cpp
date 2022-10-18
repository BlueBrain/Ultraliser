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

#include "TriangleOperations.h"

#define ULTRALISER_FLOAT_MAX    1e30f
#define ULTRALISER_FLOAT_MIN   -1e30f

namespace Ultraliser
{

bool isPointInSphere(const Vector3f& p, const Vector3f& center, const float& radius)
{
    return (p - center).abs() < radius;
}

bool isTriangleInSphere(Vector3f p0, Vector3f p1, Vector3f p2,
                        const Vector3f& center, const float& radius)
{
    // If any of the points is not in the sphere, return false
    if (!isPointInSphere(p0, center, radius))
        return false;
    if (!isPointInSphere(p1, center, radius))
        return false;
    if (!isPointInSphere(p2, center, radius))
        return false;

    // Otherwise, return true
    return true;
}

bool isTriangleIntersectingSphere(Vector3f p0, Vector3f p1, Vector3f p2,
                        const Vector3f& center, const float& radius)
{
    // If any of the points is not in the sphere, return false
    if (isPointInSphere(p0, center, radius))
        return true;
    if (isPointInSphere(p1, center, radius))
        return true;
    if (isPointInSphere(p2, center, radius))
        return true;

    // Otherwise, return false
    return false;
}

Vector3f computeNormal(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Construct two vectors
    const Vector3f a = p1 - p0;
    const Vector3f b = p2 - p0;

    // Compute the normal, cross product and then normalize
    Vector3f normal = Vector3f::cross(a, b);
    normal.normalize();
    return normal;
}

float computeTriangleSurfaceArea(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Construct two vectors
    Vector3f a = p1 - p0;
    Vector3f b = p2 - p0;

    // Compute the area
    float vx = (a.y() * b.z()) - (a.z() * b.y());
    float vy = (a.z() * b.x()) - (a.x() * b.z());
    float vz = (a.x() * b.y()) - (a.y() * b.x());

    const float area = 0.5f * std::sqrt(vx * vx + vy * vy + vz * vz);
    return  area;
}

float computeTriangleSignedVolume(Vector3f p0, Vector3f p1, Vector3f p2)
{
    return Vector3f::dot(p0, Vector3f::cross(p1, p2)) / 6.0f;
}

float computeTriangleInRadius(Vector3f p0, Vector3f p1, Vector3f p2)
{
    const Vector3f edge0 = p1 - p0;
    const Vector3f edge1 = p2 - p1;
    const Vector3f edge2 = p0 - p2;

    const float area = computeTriangleSurfaceArea(p0, p1, p2);

    return (2.f * area) / (edge0.abs() + edge1.abs() + edge2.abs());
}

float computeTriangleCircumradius(Vector3f p0, Vector3f p1, Vector3f p2)
{
    const Vector3f edge0 = p1 - p0;
    const Vector3f edge1 = p2 - p1;
    const Vector3f edge2 = p0 - p2;

    const float numerator = edge0.abs() * edge1.abs() * edge2.abs();
    const float denominator = 2.f * computeTriangleInRadius(p0, p1, p2) *
            (edge0.abs() + edge1.abs() + edge2.abs());

    return numerator / denominator;
}

float computeTriangleRadiusRatio(Vector3f p0, Vector3f p1, Vector3f p2)
{
    return computeTriangleCircumradius(p0, p1, p2) /
            (2 * computeTriangleInRadius(p0, p1, p2));
}

float computeTriangleAspectRatio(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Get the edges from the points
    float a = (p0 - p1).abs();
    float b = (p1 - p2).abs();
    float c = (p2 - p0).abs();

    // Compute s
    float s = 0.5f * (a + b + c);

    // Return the aspect ratio
    return (a * b * c) / (8 * (s - a) * (s - b) * (s - c));
}

float computeTriangleAspectFrobenius(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Construct edge vectors
    Vector3f edge1 = p1 - p0;
    Vector3f edge2 = p2 - p1;
    Vector3f edge3 = p0 - p2;

    float sumRMS = edge1.absSquared() + edge2.absSquared() + edge3.absSquared();

    // Find two times the area of the triangle by cross product
    float twiceArea = ((edge1 * (-edge3)).abs());

    if(twiceArea == 0.0)
        return 1e10;

    float aspect = (sumRMS / (2.f * sqrt (3.f) * (twiceArea)));

    if (aspect > 0)
        return std::min(aspect, 1e10f);
    return std::max(aspect, -1e10f);
}

float computeTriangleConditionNumber(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Construct edge vectors
    const Vector3f e1 = p2 - p1;
    const Vector3f e2 = p0 - p2;

    const float numerator = Vector3f::dot(e2, e2) +
                            Vector3f::dot(e1, e1) +
                            Vector3f::dot(e1, e2);

    const float surfaceArea = computeTriangleSurfaceArea(p0, p1, p2);
    const float denominator = 2.f * sqrt(3.f) * surfaceArea;

    return numerator / denominator;
}

float computeTriangleScaledJacobian(Vector3f p0, Vector3f p1, Vector3f p2)
{
    static const double detw = 2./sqrt(3.0);

    Vector3f edge[3];
    edge[0] = Vector3f(p1.x() - p0.x(), p1.y() - p0.y(), p1.z() - p0.z());
    edge[1] = Vector3f(p2.x() - p0.x(), p2.y() - p0.y(), p2.z() - p0.z());
    edge[2] = Vector3f(p2.x() - p1.x(), p2.y() - p1.y(), p2.z() - p1.z());

    Vector3f first = edge[1] - edge[0];
    Vector3f second = edge[2] - edge[0];

    Vector3f cross = first * second;
    float jacobian = cross.abs();

    float maxEdgeLengthProduct = std::max(edge[0].abs() * edge[1].abs(),
            std::max(edge[1].abs() * edge[2].abs(), edge[0].abs() * edge[2].abs()));

    if (maxEdgeLengthProduct < ULTRALISER_FLOAT_MIN)
        return 0.f;

    jacobian *= detw;
    jacobian /= maxEdgeLengthProduct;

    if( jacobian > 0 )
        return std::min(jacobian, ULTRALISER_FLOAT_MAX);
    return std::max(jacobian, -ULTRALISER_FLOAT_MAX);
}

float computeTriangleDistrotion(Vector3f p0, Vector3f p1, Vector3f p2)
{
    const Vector3f p00(-1.f, -sqrt(3.f) / 3.f, 0.f);
    const Vector3f p11(1.f, -sqrt(3.f) / 3.f, 0.f);
    const Vector3f p22(0.f, 2.f * sqrt(3.f) / 3.f, 0.f);

    const float areaMaster = sqrt(3.f);
    return (areaMaster * std::abs(computeTriangleScaledJacobian(p0, p1, p2))) /
            computeTriangleSurfaceArea(p0, p1, p2);
}

float computeTriangleEdgeRatio(Vector3f p0, Vector3f p1, Vector3f p2)
{
    // Get the edges from the points
    float a = (p0 - p1).abs();
    float b = (p1 - p2).abs();
    float c = (p2 - p0).abs();

    const float largestEdge = std::max(a, (std::max(b, c)));
    const float smallestEdge = std::min(a, (std::min(b, c)));

    // Return the edge ratio
    const float ratio = largestEdge / smallestEdge;
    return  ratio;
}

void computeTriangleAngles(Vector3f p0, Vector3f p1, Vector3f p2,
                           float& angle1, float &angle2, float &angle3)
{
    Vector3f a, b;

    a = p1 - p0;
    b = p2 - p0;
    angle1 = std::acos(Vector3f::dot(a, b) / (a.abs() * b.abs()));

    a = p0 - p1;
    b = p2 - p1;
    angle2 = std::acos(Vector3f::dot(a, b) / (a.abs() * b.abs()));

    a = p0 - p2;
    b = p1 - p2;
    angle3 = std::acos(Vector3f::dot(a, b) / (a.abs() * b.abs()));
}

float computeTriangleMinAngle(Vector3f p0, Vector3f p1, Vector3f p2)
{
    float angle1, angle2, angle3;
    computeTriangleAngles(p0, p1, p2, angle1, angle2, angle3);
    return RAD2DEG(std::min(angle1, std::min(angle2, angle3)));
}

float computeTriangleMaxAngle(Vector3f p0, Vector3f p1, Vector3f p2)
{
    float angle1, angle2, angle3;
    computeTriangleAngles(p0, p1, p2, angle1, angle2, angle3);
    return RAD2DEG(std::max(angle1, std::min(angle2, angle3)));
}

float computeTriangleRelativeSizeSquared(Vector3f p0, Vector3f p1, Vector3f p2,
                                         const float averageArea)
{
    const float r = computeTriangleSurfaceArea(p0, p1, p2) / averageArea;

    if (r == 0.f)
    {
        LOG_WARNING("r = 0");
        return 0.f;
    }

    const float minValue = std::min(r, 1.f / r);
    return minValue * minValue;
}

float computeTriangleShape(Vector3f p0, Vector3f p1, Vector3f p2)
{
    return 1.f / computeTriangleConditionNumber(p0, p1, p2);
}

float computeTriangleShapeAndSize(Vector3f p0, Vector3f p1, Vector3f p2,
                                  const float averageArea)
{
    return computeTriangleRelativeSizeSquared(p0, p1, p2, averageArea) *
            computeTriangleShape(p0, p1, p2);
}

}
