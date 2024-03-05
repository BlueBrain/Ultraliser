/***************************************************************************************************
 * Copyright (c) 2010
 * Department of Mathematics, UC San Diego
 *
 * Copyright (c) 2018 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Zeyun Yu < zeyu@math.ucsd.edu >
 *      Michael Holst < mholst@ccom.ucsd.edu >
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
 * This file has been adapted from GAMer under the terms of the GNU General Public License as
 * published by the Free Software Foundation of version 3. GAMer is published under the GNU General
 * Public License as published by the Free Software Foundation; either version 2 of the License, or
 * any later version. Further information are available at < http://fetk.org/codes/gamer >.
 **************************************************************************************************/

#include "Mesh.h"
#include "MeshOperations.h"
#include "MeshStatistics.h"
#include <algorithms/utilities/KdTree.h>
#include "TriangleOperations.h"
#include <utilities/Utilities.h>
#include <math/Math.h>
#include <utilities/Utilities.h>
#include <geometry/Intersection.h>
#include <algorithms/geometry/Geometry.h>

#define ANGLE_ERROR                         0.123456789f
#define MAXIMUM_FLOAT_VALUE                 99999.0
#define MINIMUM_FLOAT_VALUE                 -99999.0
#define MAXIMUM_DOUBLE_VALUE                99999.0
#define MINIMUM_DOUBLE_VALUE                -99999.0
#define VERTEX_DELETION_VALUE               -998.3824223883588f

namespace Ultraliser
{

void Mesh::_destroyVertexMarkers()
{
    if (_vertexMarkers.size() > 0)
    {
        _vertexMarkers.clear();
        _vertexMarkers.shrink_to_fit();
    }
}

void Mesh::_resetVertexMarkers()
{
    if (_vertexMarkers.size() == 0 || _vertexMarkers.size() < _numberVertices)
    {
        _destroyVertexMarkers();
        _vertexMarkers.resize(_numberVertices);
    }

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        _vertexMarkers[i] = 0;
    }
}

void Mesh::_destroyNeighborlist()
{
    if (_neighborList != nullptr)
    {
        // Release the single neighbors
        // OMP_PARALLEL_FOR
        for (size_t i = 0; i < _numberVertices; ++i)
        {
            NeighborTriangle* firstNGR = nullptr;
            NeighborTriangle* auxNGR = nullptr;

            firstNGR = _neighborList[i];
            while (firstNGR != nullptr)
            {
                auxNGR = firstNGR->next;
                delete firstNGR;
                firstNGR = auxNGR;
            }
        }

        // Free the array of pointers
        delete[] _neighborList;
        _neighborList = nullptr;
    }
}

void Mesh::_createNeighbourList()
{
    // Destroy exsisting neigborlist, if found
    _destroyNeighborlist();

    // Create an array of NeighborTriangle, used to store
    NeighborTriangle **neighborList = new NeighborTriangle*[_numberVertices];

    // Initialize the neighbor list
    for (size_t i = 0; i < _numberVertices; ++i)
        neighborList[i] = nullptr;

    // Start the timer
    TIMER_SET;

    LOOP_STARTS("Creating Neighbour List");
    PROGRESS_SET;
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        PROGRESS_UPDATE;
        LOOP_PROGRESS(PROGRESS, _numberTriangles);

        // Iterate over the triangless and collect line segments (a, b) and
        // its connection to a triangles (c).
        // Save the line segment so it forms a counter clockwise triangle with the
        // origin vertices
        const int64_t a = _triangles[i][0];
        const int64_t b = _triangles[i][1];
        const int64_t c = _triangles[i][2];

        NeighborTriangle *firstNGR = new NeighborTriangle();
        firstNGR->a = b;
        firstNGR->b = c;
        firstNGR->c = i;
        firstNGR->next = neighborList[a];
        neighborList[a] = firstNGR;

        firstNGR = new NeighborTriangle();
        firstNGR->a = c;
        firstNGR->b = a;
        firstNGR->c = i;
        firstNGR->next = neighborList[b];
        neighborList[b] = firstNGR;

        firstNGR = new NeighborTriangle();
        firstNGR->a = a;
        firstNGR->b = b;
        firstNGR->c = i;
        firstNGR->next = neighborList[c];
        neighborList[c] = firstNGR;
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    NeighborTriangle *firstNGR, *secondNGR, *auxNGR;

    // Order the neighbors so they are connected counter clockwise
    TIMER_RESET;
    LOOP_STARTS("Ordering Vertices");
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        LOOP_PROGRESS(i, _numberVertices);

        firstNGR = neighborList[i];
        const int64_t c = firstNGR->a;

        while (firstNGR != nullptr)
        {
            const int64_t a = firstNGR->a;
            const int64_t b = firstNGR->b;

            secondNGR = firstNGR->next;

            while (secondNGR != nullptr)
            {
                const int64_t a0 = secondNGR->a;
                const int64_t b0 = secondNGR->b;

                if (a0 == b && b0 != a)
                {
                    auxNGR = firstNGR;

                    while (auxNGR != nullptr)
                    {
                        if (auxNGR->next == secondNGR)
                        {
                            auxNGR->next = secondNGR->next;
                            break;
                        }
                        auxNGR = auxNGR->next;
                    }

                    auxNGR = firstNGR->next;
                    firstNGR->next = secondNGR;
                    secondNGR->next = auxNGR;
                    break;
                }

                secondNGR = secondNGR->next;
            }
            if (firstNGR->next == nullptr)
            {
                if (firstNGR->b != c)
                {
                    LOG_WARNING("Some polygons are not closed, Vertices: [%d - %d]", firstNGR->b, c);
                    LOG_WARNING("[%f, %f, %f]", F2D(_vertices[firstNGR->b].x()),
                                                F2D(_vertices[firstNGR->b].y()),
                                                F2D(_vertices[firstNGR->b].z()));
                    LOG_WARNING("[%f, %f, %f]", F2D(_vertices[c].x()),
                                                F2D(_vertices[c].y()),
                                                F2D(_vertices[c].z()));
                }
            }

            firstNGR = firstNGR->next;
        }
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Attach the neigborlist to the surfmesh
    this->_neighborList = neighborList;
}

void Mesh::_selectVerticesInROI(const ROIs& regions)
{
    // Start the timer
    TIMER_SET;
    LOG_STATUS("Selecting ROI Vertices: Total [ %d ]", _numberVertices);

    // Reset the vertex markers
    _resetVertexMarkers();

    // Select the triangles that are located within the ROI
    LOOP_STARTS("Labeling Vertices");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        // Get the vertex
        const auto& vertex = _vertices[i];

        // Make sure that you cover all the regions
        for (size_t j = 0; j < regions.size(); ++j)
        {
            const auto region = regions[j];

            // If the vertex is in the sphere
            if (isPointInSphere(vertex, region->center, region->radius))
            {
                // Mark the vertices to avoid losing them
                _vertexMarkers[i] = 1;
            }
        }

        // Update the progress bar
        LOOP_PROGRESS_FRACTION(PROGRESS, _numberVertices);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    size_t numberSelectedVertices = 0;
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        if (_vertexMarkers[i] == 1)
        {
            numberSelectedVertices++;
        }
    }

    LOG_SUCCESS("[ %d ] vertices are selected for [ %d ] ROIs",
                numberSelectedVertices, regions.size());
}

void Mesh::optimizeAdapttivelyWithROI(const size_t &optimizationIterations,
                                      const size_t &smoothingIterations,
                                      const float &flatCoarseFactor,
                                      const float &denseFactor,
                                      const ROIs &regions)
{
    LOG_TITLE("Adaptive Mesh Optimization (ROI)");

    // Starting the timer
    TIMER_SET;

    // Select the vertices in the ROI
    _selectVerticesInROI(regions);

    // Coarse flat
    coarseFlat(flatCoarseFactor, optimizationIterations);

    // Coarse dense
    coarseDense(denseFactor, optimizationIterations);

    // Destroy the vertex markers
    _destroyVertexMarkers();

    // Smooth again
    smooth(15, 150, smoothingIterations);
    smoothNormals();

    // Statistics
    LOG_STATUS_IMPORTANT("Total Optimization");
    LOG_STATS(GET_TIME_SECONDS);
}

bool Mesh::isValidVertex(const Vector3f& v)
{
    if(isEqual(v.x(), VERTEX_DELETION_VALUE) ||
       isEqual(v.y(), VERTEX_DELETION_VALUE) ||
       isEqual(v.z(), VERTEX_DELETION_VALUE))
        return false;
    return true;
}

float Mesh::getAngleSurfaceOnly(const int64_t &a, const int64_t &b, const int64_t &c,
                                bool &angleError)
{
    const float ax = _vertices[a].x();
    const float ay = _vertices[a].y();
    const float az = _vertices[a].z();
    const float bx = _vertices[b].x();
    const float by = _vertices[b].y();
    const float bz = _vertices[b].z();
    const float cx = _vertices[c].x();
    const float cy = _vertices[c].y();
    const float cz = _vertices[c].z();

    const float length1 = (ax - bx) * (ax - bx) + (ay - by) * (ay - by) + (az - bz) * (az - bz);
    const float length2 = (ax - cx) * (ax - cx) + (ay - cy) * (ay - cy) + (az - cz) * (az - cz);
    const float length3 = (bx - cx) * (bx - cx) + (by - cy) * (by - cy) + (bz - cz) * (bz - cz);

    float angle;
    if (isZero(length1) || isZero(length2))
    {
        angle = ANGLE_ERROR;
        angleError = true;
    }
    else
    {
        angle = 0.5f * (length1 + length2 - length3) / D2F(sqrt(length1 * length2));
        angle = std::acos(angle) * 180.0 / 3.14159265358979f;
    }

    // Return the angle
    return angle;
}

float Mesh::computeDotProduct(const int64_t &a, const int64_t &b, const int64_t &c)
{
    // Construct the vector
    float bx = _vertices[b].x() - _vertices[a].x();
    float by = _vertices[b].y() - _vertices[a].y();
    float bz = _vertices[b].z() - _vertices[a].z();

    // Compute the length
    float length = sqrt(bx * bx + by * by + bz * bz);

    // Normalize
    if (length > 0.0)
    {
        bx /= length;
        by /= length;
        bz /= length;
    }

    // Construct the vector
    float cx = _vertices[c].x()- _vertices[a].x();
    float cy = _vertices[c].y()- _vertices[a].y();
    float cz = _vertices[c].z()- _vertices[a].z();

    // Compute the length
    length = sqrt(cx * cx + cy * cy + cz * cz);

    // Normalize
    if (length > 0.0)
    {
        cx /= length;
        cy /= length;
        cz /= length;
    }

    // Compute the dot product and return the value
    return bx * cx + by * cy + bz * cz;
}

Vector3f Mesh::computeCrossProduct(const int64_t &a, const int64_t &b, const int64_t &c)
{
    float bx = _vertices[b].x() - _vertices[a].x();
    float by = _vertices[b].y() - _vertices[a].y();
    float bz = _vertices[b].z() - _vertices[a].z();

    // Normalize
    float length = std::sqrt(bx * bx + by * by + bz * bz);
    if (length > 0.0)
    {
        bx /= length;
        by /= length;
        bz /= length;
    }

    float cx = _vertices[c].x() - _vertices[a].x();
    float cy = _vertices[c].y() - _vertices[a].y();
    float cz = _vertices[c].z() - _vertices[a].z();

    // Normalize
    length = std::sqrt(cx * cx + cy * cy + cz * cz);
    if (length > 0.0)
    {
        cx /= length;
        cy /= length;
        cz /= length;
    }

    float gx = cy * bz - cz * by;
    float gy = cz * bx - cx * bz;
    float gz = cx * by - cy * bx;

    // Normalize
    length = sqrt(gx * gx + gy * gy + gz * gz);
    if (length > 0.0)
    {
        gx /= length;
        gy /= length;
        gz /= length;
    }

    return Vector3f(gx, gy, gz);
}

bool Mesh::checkFlipAction(const int64_t &a, const int64_t &b, const int64_t &c, const int64_t &d)
{

    // Smaller angle criterion
    float minAngle1 = MINIMUM_FLOAT_VALUE;

    float angle = computeDotProduct(a, b, c);
    if (angle > minAngle1)
        minAngle1 = angle;

    angle = computeDotProduct(a, b, d);
    if (angle > minAngle1)
        minAngle1 = angle;

    angle = computeDotProduct(b, a, c);
    if (angle > minAngle1)
        minAngle1 = angle;

    angle = computeDotProduct(b, a, d);
    if (angle > minAngle1)
        minAngle1 = angle;

    float minAngle2 = MINIMUM_FLOAT_VALUE;
    angle = computeDotProduct(c, a, d);
    if (angle > minAngle2)
        minAngle2 = angle;

    angle = computeDotProduct(c, b, d);
    if (angle > minAngle2)
        minAngle2 = angle;

    angle = computeDotProduct(d, a, c);
    if (angle > minAngle2)
        minAngle2 = angle;

    angle = computeDotProduct(d, b, c);
    if (angle > minAngle2)
        minAngle2 = angle;

    // Return the check result based on the angles
    if (minAngle1 > minAngle2)
        return true;
    else
        return false;
}

void Mesh::edgeFlipping(int64_t index)
{
    int64_t a,b,c;

    NeighborTriangle* auxNGR1 = nullptr;
    NeighborTriangle* auxNGR2 = nullptr;
    NeighborTriangle* auxNGR = nullptr;

    NeighborTriangle* firstNGR = _neighborList[index];

    int64_t number = 0;
    while (firstNGR != nullptr)
    {
        number = 0;
        auxNGR = _neighborList[index];

        while (auxNGR != nullptr)
        {
            number++;
            auxNGR = auxNGR->next;
        }

        if (number <= 3)
        {
            if (number > 0)
            {
                float ax = 0.0;
                float ay = 0.0;
                float az = 0.0;

                auxNGR = _neighborList[index];
                while (auxNGR != nullptr)
                {
                    a = auxNGR->a;
                    ax += _vertices[a].x();
                    ay += _vertices[a].y();
                    az += _vertices[a].z();
                    auxNGR = auxNGR->next;
                }

                _vertices[index].x()= ax / (1.0 * number);
                _vertices[index].y()= ay / (1.0 * number);
                _vertices[index].z()= az / (1.0 * number);
            }
            return;
        }

        a = firstNGR->a;
        b = firstNGR->b;

        NeighborTriangle* secondNGR = firstNGR->next;
        if (secondNGR == nullptr)
            secondNGR = _neighborList[index];

        c = secondNGR->b;

        char flipFlag = 1;
        number = 0;
        auxNGR = _neighborList[b];
        while (auxNGR != nullptr)
        {
            number++;
            auxNGR = auxNGR->next;
        }
        if (number <= 3)
            flipFlag = 0;

        auxNGR = _neighborList[a];
        while (auxNGR != nullptr)
        {
            if (auxNGR->a == c)
                flipFlag = 0;
            auxNGR = auxNGR->next;
        }
        auxNGR = _neighborList[c];
        while (auxNGR != nullptr)
        {
            if (auxNGR->a == a)
                flipFlag = 0;
            auxNGR = auxNGR->next;
        }

        if (flipFlag)
        {
            if (checkFlipAction(index, b, a, c))
            {
                int64_t f1 = firstNGR->c;
                int64_t f2 = secondNGR->c;

                // Update triangles info
                _triangles[f1][0] = index;
                _triangles[f1][1] = a;
                _triangles[f1][2] = c;

                // Switch a and c here to make the triangles
                // normal point outward
                _triangles[f2][0] = b;
                _triangles[f2][1] = c;
                _triangles[f2][2] = a;

                // Delete the entries in neighbor lists
                firstNGR->b = c;
                if (firstNGR->next == nullptr)
                    _neighborList[index] = _neighborList[index]->next;
                else
                    firstNGR->next = secondNGR->next;
                auxNGR1 = secondNGR;

                auxNGR = _neighborList[b];
                while (auxNGR != nullptr)
                {
                    if (auxNGR->b == index)
                        break;
                    auxNGR = auxNGR->next;
                }

                if (auxNGR == nullptr)
                    LOG_WARNING("auxNGR == nullptr @ [%d]\n", index);

                if (auxNGR->a == c)
                {
                    auxNGR->b = a;
                    auxNGR->c = f2;

                    if (auxNGR->next == nullptr)
                    {
                        secondNGR = _neighborList[b];
                        _neighborList[b] = secondNGR->next;
                    }
                    else
                    {
                        secondNGR = auxNGR->next;
                        auxNGR->next = secondNGR->next;
                    }
                    auxNGR2 = secondNGR;
                }
                else
                {
                    LOG_WARNING("Deletion Error [%d : %d %d %d]",
                                index, a, b, c);
                    LOG_WARNING("[%f, %f, %f]", F2D(_vertices[index].x()),
                                F2D(_vertices[index].y()),
                                F2D(_vertices[index].z()));
                }

                // Add the entries in neighbor lists
                auxNGR = _neighborList[a];

                while (auxNGR != nullptr)
                {
                    if ((auxNGR->a == index && auxNGR->b == b) ||
                            (auxNGR->a == b && auxNGR->b == index))
                        break;
                    auxNGR = auxNGR->next;
                }

                // Assume neigbors are stored counter clockwise
                if (auxNGR->a == b && auxNGR->b == index)
                {
                    auxNGR->b = c;
                    auxNGR->c = f2;
                    auxNGR1->a = c;
                    auxNGR1->b = index;
                    auxNGR1->c = f1;
                    auxNGR1->next = auxNGR->next;
                    auxNGR->next = auxNGR1;
                }
                else
                    LOG_ERROR("Neighbour issue 1");

                auxNGR = _neighborList[c];
                while (auxNGR != nullptr)
                {
                    if ((auxNGR->a == index && auxNGR->b == b) ||
                            (auxNGR->a == b && auxNGR->b == index))
                        break;
                    auxNGR = auxNGR->next;
                }

                // Assume neigbors are stored counter clockwise
                if (auxNGR->a == index && auxNGR->b == b)
                {
                    auxNGR->b = a;
                    auxNGR->c = f1;
                    auxNGR2->a = a;
                    auxNGR2->b = b;
                    auxNGR2->c = f2;
                    auxNGR2->next = auxNGR->next;
                    auxNGR->next = auxNGR2;
                }
                else
                    LOG_ERROR("Neighbour issue 2");
            }
        }

        firstNGR = firstNGR->next;
    }
}

void Mesh::computeAngles(float* computedMinAngle, float* computedMaxAngle,
                         int64_t *numSmall, int64_t *numLarge,
                         int64_t maxMinAngle, int64_t minMaxAngle)
{
    // Initialize the angles
    float minAngle = std::numeric_limits< float >::max();
    float maxAngle = std::numeric_limits< float >::lowest();

    // Initialize the numbers
    int64_t smallNumber = 0;
    int64_t largeNumber = 0;

    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        const int64_t a = _triangles[i][0];
        const int64_t b = _triangles[i][1];
        const int64_t c = _triangles[i][2];

        // Compute the angle ab <-> bc
        bool angleError = false;
        float angle = getAngleSurfaceOnly(a, b, c, angleError);

        if (!angleError)
        {
            if (angle < minAngle)
                minAngle = angle;
            if (angle > maxAngle)
                maxAngle = angle;
            if (angle < maxMinAngle)
                smallNumber ++;
            if (angle > minMaxAngle)
                largeNumber++;
        }

        // Compute the angle ba <-> ac
        angle = getAngleSurfaceOnly(b, a, c, angleError);

        if (!angleError)
        {
            if (angle < minAngle)
                minAngle = angle;
            if (angle > maxAngle)
                maxAngle = angle;
            if (angle < maxMinAngle)
                smallNumber ++;
            if (angle > minMaxAngle)
                largeNumber++;
        }

        // Compute the angle ca <-> ab
        angle = getAngleSurfaceOnly(c, a, b, angleError);

        if (!angleError)
        {
            if (angle < minAngle)
                minAngle = angle;
            if (angle > maxAngle)
                maxAngle = angle;
            if (angle < maxMinAngle)
                smallNumber ++;
            if (angle > minMaxAngle)
                largeNumber++;
        }
    }

    // Return the computed values
    *computedMinAngle = minAngle;
    *computedMaxAngle = maxAngle;
    *numSmall = smallNumber ;
    *numLarge = largeNumber;
}

Vector3f Mesh::getPositionSurfaceOnly(const float &x, const float &y, const float &z,
                                      const int64_t &a, const int64_t &b, const int64_t &c)
{
    // Access the coordinates of the three vertices of the triangle
    // First vertex
    float ax = _vertices[a].x();
    float ay = _vertices[a].y();
    float az = _vertices[a].z();

    // Second vertex
    float bx = _vertices[b].x();
    float by = _vertices[b].y();
    float bz = _vertices[b].z();

    // Third vertex
    float cx = _vertices[c].x();
    float cy = _vertices[c].y();
    float cz = _vertices[c].z();

    // Subtraction
    bx -= ax; by -= ay; bz -= az;

    // Computing the distance
    float distance = sqrt(bx * bx + by * by + bz * bz);

    // Normalize
    if (distance > 0.0)
    {
        bx /= distance;
        by /= distance;
        bz /= distance;
    }

    // Subtraction
    cx -= ax; cy -= ay; cz -= az;

    // Compute the distance
    distance = sqrt(cx*cx+cy*cy+cz*cz);

    // Normalize
    if (distance > 0.0)
    {
        cx /= distance;
        cy /= distance;
        cz /= distance;
    }

    // Center
    float tx = 0.5f * (cx + bx);
    float ty = 0.5f * (cy + by);
    float tz = 0.5f * (cz + bz);

    // Compute the distance
    distance = sqrt(tx * tx + ty * ty + tz * tz);

    // Normalize
    if (distance > 0.0)
    {
        tx /= distance;
        ty /= distance;
        tz /= distance;
    }

    // Compute the distance
    float xx = by * cz - bz * cy;
    float yy = bz * cx - bx * cz;
    float zz = bx * cy - by * cx;
    distance = sqrt(xx * xx + yy * yy + zz * zz);

    // Normalize
    if (distance > 0)
    {
        xx /= distance;
        yy /= distance;
        zz /= distance;
    }

    // Second point
    bx = xx;
    by = yy;
    bz = zz;

    // Compute the distance
    distance = tx * (x - ax) + ty * (y - ay) + tz * (z - az);
    xx = distance * tx + ax;
    yy = distance * ty + ay;
    zz = distance * tz + az;

    // Compute the distance
    distance = bx * (x - xx) + by * (y - yy) + bz * (z - zz);

    // Compute the position
    Vector3f position;
    position.x()= distance * bx + xx;
    position.y()= distance * by + yy;
    position.z()= distance * bz + zz;

    // Return the position
    return position;
}

Vector3f Mesh::getVertexNormal(const int64_t &index)
{
    // Get vertex coordinates
    const float x = _vertices[index].x();
    const float y = _vertices[index].y();
    const float z = _vertices[index].z();

    int64_t number = 0;

    // Computed normal
    Vector3f normal;

    NeighborTriangle* firstNGR = _neighborList[index];
    while (firstNGR != nullptr)
    {
        int64_t a = firstNGR->a;
        int64_t b = firstNGR->b;

        float ax = _vertices[a].x()- x;
        float ay = _vertices[a].y()- y;
        float az = _vertices[a].z()- z;

        // Compute the length
        float length = sqrt(ax * ax + ay * ay + az * az);

        // Normalize
        if (length > 0.0)
        {
            ax /= length;
            ay /= length;
            az /= length;
        }

        float bx = _vertices[b].x()-x;
        float by = _vertices[b].y()-y;
        float bz = _vertices[b].z()-z;

        // Compute the length
        length = sqrt(bx * bx + by * by + bz * bz);

        // Normalize
        if (length > 0.0)
        {
            bx /= length;
            by /= length;
            bz /= length;
        }

        float gx = ay*bz-az*by;
        float gy = az*bx-ax*bz;
        float gz = ax*by-ay*bx;

        // Compute the length
        length = sqrt(gx * gx + gy * gy + gz * gz);

        // Normalize
        if (length > 0.0)
        {
            gx /= length;
            gy /= length;
            gz /= length;
        }

        // Compute the length
        length = normal.x()* gx + normal.y()* gy + normal.z()* gz;

        // Shift
        if (length < 0.0)
        {
            gx = -gx;
            gy = -gy;
            gz = -gz;
        }

        // Update the normals
        normal.x()+= gx;
        normal.y()+= gy;
        normal.z()+= gz;

        number++;
        firstNGR = firstNGR->next;
    }

    if (number > 0)
    {
        normal.x()/= (1.0 * number);
        normal.y()/= (1.0 * number);
        normal.z()/= (1.0 * number);

        // Compute the length
        float length = sqrt(normal.x()* normal.x()+
                            normal.y()* normal.y()+
                            normal.z()* normal.z());

        // Normalize
        if (length > 0.0)
        {
            normal.x()/= length;
            normal.y()/= length;
            normal.z()/= length;
        }
    }
    else
    {
        normal.x()= 0;
        normal.y()= 0;
        normal.z()= 0;
    }

    return normal;
}

EigenVector Mesh::computeEigenVector(const int64_t &index0,
                                     Vector3f *eigenValue,
                                     float *maxAngleComputed)
{
    // Compute the normal
    Vector3f normal = getVertexNormal(index0);

    double A[3][3];
    A[0][0] = F2D(normal.x() * normal.x());
    A[0][1] = F2D(normal.x() * normal.y());
    A[0][2] = F2D(normal.x() * normal.z());
    A[1][1] = F2D(normal.y() * normal.y());
    A[1][2] = F2D(normal.y() * normal.z());
    A[2][2] = F2D(normal.z() * normal.z());

    int64_t startPointer = 0;
    int64_t endPointer = 1;

    int64_t indexArray[333];
    indexArray[startPointer] = index0;

    int64_t distArray[333];
    distArray[startPointer] = 0;

    // Maximum angle
    float maxAngle = MAXIMUM_FLOAT_VALUE;

    // Copy the normal
    Vector3f normal0(normal.x(), normal.y(), normal.z());

    int64_t index, dist, visited;
    while (startPointer < endPointer)
    {
        index = indexArray[startPointer];
        dist = distArray[startPointer];
        startPointer++;

        // NOTE: The comparison against 2.0 is not well understood.
        if (dist < 2.0)
        {
            NeighborTriangle *firstNGR = _neighborList[index];
            while (firstNGR != nullptr)
            {
                int64_t m = firstNGR->a;
                visited = 0;
                for (int64_t n = 0; n < endPointer; n++)
                {
                    if (indexArray[n] == m)
                    {
                        visited = 1;
                        break;
                    }
                }
                if (visited == 0)
                {
                    normal = getVertexNormal(m);
                    float angle = normal0.x() * normal.x() +
                            normal0.y() * normal.y() +
                            normal0.z() * normal.z();

                    if (angle < 0.0)
                        angle = -angle;

                    if (angle < maxAngle)
                        maxAngle = angle;

                    A[0][0] += F2D(normal.x() * normal.x());
                    A[0][1] += F2D(normal.x() * normal.y());
                    A[0][2] += F2D(normal.x() * normal.z());
                    A[1][1] += F2D(normal.y() * normal.y());
                    A[1][2] += F2D(normal.y() * normal.z());
                    A[2][2] += F2D(normal.z() * normal.z());

                    indexArray[endPointer] = m;
                    distArray[endPointer] = dist+1;
                    endPointer ++;
                }
                firstNGR = firstNGR->next;
            }
        }
    }

    // Assign to the return value
    *maxAngleComputed = maxAngle;

    A[1][0] = A[0][1];
    A[2][0] = A[0][2];
    A[2][1] = A[1][2];

    double c0 = A[0][0] * A[1][1] * A[2][2] +
                A[0][1] * A[0][2] * A[1][2] * 2 -
                A[0][0] * A[1][2] * A[1][2] -
                A[1][1] * A[0][2] * A[0][2] -
                A[2][2] * A[0][1] * A[0][1];
    double c1 = A[0][0] * A[1][1] -  A[0][1] * A[0][1] +
                A[0][0] * A[2][2]- A[0][2] * A[0][2] +
                A[1][1] * A[2][2] - A[1][2] * A[1][2];
    double c2 = A[0][0] + A[1][1] + A[2][2];

    double a = (3.0 * c1 - c2 * c2) / 3.0;
    double b = (-2.0 * c2 * c2 * c2 + 9.0 * c1 * c2 - 27.0 * c0) / 27.0;
    double Q = b * b / 4.0 + a * a * a / 27.0;

    double theta = atan2(sqrt(-Q), -0.5 * b);
    double p = sqrt(0.25 * b * b - Q);

    double x1 = c2 / 3.0 + 2.0 * pow(p, 1.0 / 3.0) * cos(theta/3.0);
    double x2 = c2 / 3.0 - pow(p, 1.0 / 3.0) * (cos(theta / 3.0) + sqrt(3.0) * sin(theta / 3.0));
    double x3 = c2 / 3.0 - pow(p, 1.0 / 3.0) * (cos(theta / 3.0) - sqrt(3.0) * sin(theta / 3.0));

    double tx = 0.0;
    double ty = 0.0;
    double tz = 0.0;

    EigenVector computedEigenVector;
    if (std::isnan(x1) || std::isnan(x2) || std::isnan(x3))
    {
        computedEigenVector.isValid = false;
        return computedEigenVector;
    }
    else
    {
        tx = std::max(x1, std::max(x2,x3));
        if (isEqual(tx, x1))
        {
            if (x2 >= x3)
            {
                ty = x2;
                tz = x3;
            }
            else
            {
                ty = x3;
                tz = x2;
            }
        }
        else if (isEqual(tx, x2))
        {
            if (x1 >= x3)
            {
                ty = x1;
                tz = x3;
            }
            else
            {
                ty = x3;
                tz = x1;
            }
        }
        else if (isEqual(tx, x3))
        {
            if (x1 >= x2)
            {
                ty = x1;
                tz = x2;
            }
            else
            {
                ty = x2;
                tz = x1;
            }
        }

        x1 = tx;
        x2 = ty;
        x3 = tz;

        // Update the eigen value
        eigenValue->x() = static_cast< float >(tx);
        eigenValue->y() = static_cast< float >(ty);
        eigenValue->z() = static_cast< float >(tz);

        if (x1 > MAXIMUM_DOUBLE_VALUE || x1 < MINIMUM_DOUBLE_VALUE ||
                x2 > MAXIMUM_DOUBLE_VALUE || x2 < MINIMUM_DOUBLE_VALUE ||
                x3 > MAXIMUM_DOUBLE_VALUE || x3 < MINIMUM_DOUBLE_VALUE)
        {
            LOG_ERROR("ERROR: [%f %f %f]", x1, x2, x3);
            exit(0);
        }

        A[0][0] -= x1;
        A[1][1] -= x1;
        A[2][2] -= x1;

        double B[6];
        B[0] = A[1][1] * A[2][2] - A[1][2] * A[1][2];
        B[1] = A[0][2] * A[1][2] - A[0][1] * A[2][2];
        B[2] = A[0][0] * A[2][2] - A[0][2] * A[0][2];
        B[3] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
        B[4] = A[0][1] * A[0][2] - A[1][2] * A[0][0];
        B[5] = A[0][0] * A[1][1] - A[0][1] * A[0][1];

        c0 = B[0] * B[0]+B[1] * B[1]+B[3] * B[3];
        c1 = B[1] * B[1]+B[2] * B[2]+B[4] * B[4];
        c2 = B[3] * B[3]+B[4] * B[4]+B[5] * B[5];

        if (c0 >= c1 && c0 >= c2)
        {
            tx = B[0];
            ty = B[1];
            tz = B[3];
        }
        else if (c1 >= c0 && c1 >= c2)
        {
            tx = B[1];
            ty = B[2];
            tz = B[4];
        }
        else if (c2 >= c0 && c2 >= c1)
        {
            tx = B[3];
            ty = B[4];
            tz = B[5];
        }

        double p = sqrt(tx * tx + ty * ty + tz * tz);
        if (p > 0.0)
        {
            tx /= p;
            ty /= p;
            tz /= p;
        }

        // Update the eigen vector
        computedEigenVector.x1 = D2F(tx);
        computedEigenVector.y1 = D2F(ty);
        computedEigenVector.z1 = D2F(tz);

        A[0][0] += x1;
        A[1][1] += x1;
        A[2][2] += x1;

        A[0][0] -= x2;
        A[1][1] -= x2;
        A[2][2] -= x2;

        B[0] = A[1][1] * A[2][2] - A[1][2] * A[1][2];
        B[1] = A[0][2] * A[1][2] - A[0][1] * A[2][2];
        B[2] = A[0][0] * A[2][2] - A[0][2] * A[0][2];
        B[3] = A[0][1] * A[1][2] - A[0][2] * A[1][1];
        B[4] = A[0][1] * A[0][2] - A[1][2] * A[0][0];
        B[5] = A[0][0] * A[1][1] - A[0][1] * A[0][1];

        c0 = B[0] * B[0] + B[1] * B[1] +B[3] * B[3];
        c1 = B[1] * B[1] + B[2] * B[2] +B[4] * B[4];
        c2 = B[3] * B[3] + B[4] * B[4] +B[5] * B[5];

        if (c0 >= c1 && c0 >= c2)
        {
            tx = B[0];
            ty = B[1];
            tz = B[3];
        }
        else if (c1 >= c0 && c1 >= c2)
        {
            tx = B[1];
            ty = B[2];
            tz = B[4];
        }
        else if (c2 >= c0 && c2 >= c1)
        {
            tx = B[3];
            ty = B[4];
            tz = B[5];
        }

        p = sqrt(tx*tx+ty*ty+tz*tz);

        if (p > 0.0)
        {
            tx /= p;
            ty /= p;
            tz /= p;
        }

        computedEigenVector.x2 = D2F(tx);
        computedEigenVector.y2 = D2F(ty);
        computedEigenVector.z2 = D2F(tz);

        computedEigenVector.x3 = computedEigenVector.y1 * D2F(tz) - computedEigenVector.z1 * D2F(ty);
        computedEigenVector.y3 = computedEigenVector.z1 * D2F(tx) - computedEigenVector.x1 * D2F(tz);
        computedEigenVector.z3 = computedEigenVector.x1 * D2F(ty) - computedEigenVector.y1 * D2F(tx);

        computedEigenVector.isValid = true;
    }

    return computedEigenVector;
}

Vector3f Mesh::rotate(const float& sx, const float& sy, const float& sz,
                      const float& theta, const float& phi, const float& angle)
{

    Vector3f rotatedVector;

    float a[3][3], b[3][3];
    a[0][0] = cos(0.5f * ULTRALISER_PIF - phi) * cos(theta);
    a[0][1] = cos(0.5f * ULTRALISER_PIF - phi) * sin(theta);
    a[0][2] = -sin(0.5f * ULTRALISER_PIF - phi);
    a[1][0] = -sin(theta);
    a[1][1] =  cos(theta);
    a[1][2] = 0.0;
    a[2][0] = sin(0.5f * ULTRALISER_PIF - phi) * cos(theta);
    a[2][1] = sin(0.5f * ULTRALISER_PIF - phi) * sin(theta);
    a[2][2] = cos(0.5f * ULTRALISER_PIF - phi);

    b[0][0] = cos(0.5f * ULTRALISER_PIF - phi) * cos(theta);
    b[0][1] = -sin(theta);
    b[0][2] = sin(0.5f * ULTRALISER_PIF - phi) * cos(theta);
    b[1][0] = cos(0.5f * ULTRALISER_PIF - phi) * sin(theta);
    b[1][1] = cos(theta);
    b[1][2] = sin(0.5f * ULTRALISER_PIF - phi) * sin(theta);
    b[2][0] = -sin(0.5f * ULTRALISER_PIF - phi);
    b[2][1] = 0.0;
    b[2][2] = cos(0.5f * ULTRALISER_PIF - phi);

    const float x = a[0][0] * sx + a[0][1] * sy + a[0][2] * sz;
    const float y = a[1][0] * sx + a[1][1] * sy + a[1][2] * sz;
    const float z = a[2][0] * sx + a[2][1] * sy + a[2][2] * sz;

    const float xx = cos(angle) * x - sin(angle) * y;
    const float yy = sin(angle) * x + cos(angle) * y;
    const float zz = z;

    rotatedVector.x() = b[0][0] * xx + b[0][1] * yy + b[0][2] * zz;
    rotatedVector.y() = b[1][0] * xx + b[1][1] * yy + b[1][2] * zz;
    rotatedVector.z() = b[2][0] * xx + b[2][1] * yy + b[2][2] * zz;

    return(rotatedVector);

}

void Mesh::subdividePolygin(NeighborTriangle *starNGR,
                            int64_t *triangleAvailableList,
                            int64_t *triangleAvailableIndex)
{
    NeighborTriangle *firstNGR = nullptr, *secondNGR = nullptr;
    NeighborTriangle *auxNGR, *firstCopyNGR,*secondCopyNGR;

    int64_t minNum, degree;
    int64_t faceIndex, number;
    int64_t a,b,c;

    number = 1;
    auxNGR = starNGR;
    while (auxNGR->next != starNGR)
    {
        number++;
        auxNGR = auxNGR->next;
    }

    if (number < 3)
    {
        LOG_ERROR("Number of nodes is less than 3");
    }
    if (number == 3)
    {
        a = starNGR->a;
        auxNGR = starNGR->next;
        delete starNGR;
        starNGR = auxNGR;

        b = starNGR->a;
        auxNGR = starNGR->next;
        delete starNGR;
        starNGR = auxNGR;

        c = starNGR->a;
        auxNGR = starNGR->next;
        delete starNGR;
        starNGR = auxNGR;

        faceIndex = triangleAvailableList[*triangleAvailableIndex];
        _triangles[faceIndex][0] = a;
        _triangles[faceIndex][1] = b;
        _triangles[faceIndex][2] = c;
        *triangleAvailableIndex += 1;

        firstNGR = new NeighborTriangle();
        firstNGR->a = b;
        firstNGR->b = c;
        firstNGR->c = faceIndex;
        firstNGR->next = _neighborList[a];
        _neighborList[a] = firstNGR;

        firstNGR = new NeighborTriangle();
        firstNGR->a = c;
        firstNGR->b = a;
        firstNGR->c = faceIndex;
        firstNGR->next = _neighborList[b];
        _neighborList[b] = firstNGR;

        firstNGR = new NeighborTriangle();
        firstNGR->a = a;
        firstNGR->b = b;
        firstNGR->c = faceIndex;
        firstNGR->next = _neighborList[c];
        _neighborList[c] = firstNGR;
    }
    else
    {
        auxNGR = starNGR;
        minNum = auxNGR->b;
        firstNGR = auxNGR;
        auxNGR = auxNGR->next;
        while (auxNGR != starNGR)
        {
            degree = auxNGR->b;
            if (degree < minNum)
            {
                minNum = degree;
                firstNGR = auxNGR;
            }
            auxNGR = auxNGR->next;
        }

        minNum = 99999;
        auxNGR = starNGR;
        if (auxNGR != firstNGR &&
                auxNGR != firstNGR->next &&
                auxNGR->next != firstNGR)
        {
            minNum = auxNGR->b;
            secondNGR = auxNGR;
        }
        auxNGR = auxNGR->next;
        while (auxNGR != starNGR)
        {
            degree = auxNGR->b;
            if (auxNGR != firstNGR && auxNGR != firstNGR->next &&
                    auxNGR->next != firstNGR && degree < minNum)
            {
                minNum = degree;
                secondNGR = auxNGR;
            }
            auxNGR = auxNGR->next;
        }

        firstNGR->b += 1;
        secondNGR->b += 1;

        firstCopyNGR = new NeighborTriangle();
        firstCopyNGR->a = firstNGR->a;
        firstCopyNGR->b = firstNGR->b;

        secondCopyNGR = new NeighborTriangle();
        secondCopyNGR->a = secondNGR->a;
        secondCopyNGR->b = secondNGR->b;

        auxNGR = firstNGR;

        while (auxNGR->next != firstNGR)
            auxNGR = auxNGR->next;

        auxNGR->next = firstCopyNGR;
        firstCopyNGR->next = secondCopyNGR;
        secondCopyNGR->next = secondNGR->next;
        secondNGR->next = firstNGR;

        subdividePolygin(firstNGR, triangleAvailableList,
                         triangleAvailableIndex);
        subdividePolygin(firstCopyNGR, triangleAvailableList,
                         triangleAvailableIndex);
    }
}

void Mesh::moveVertexAlongSurface(int64_t index)
{
    // Get the coordinates of the vertex
    float x = _vertices[index].x();
    float y = _vertices[index].y();
    float z = _vertices[index].z();

    // New vertex position
    float nx = 0.0;
    float ny = 0.0;
    float nz = 0.0;

    // Vertex weight, initially set to zero
    float weight = 0;

    NeighborTriangle *firstNGR = _neighborList[index];

    // If there is a neighbour
    while (firstNGR != nullptr)
    {
        // Get the first points
        int64_t a = firstNGR->a;
        int64_t b = firstNGR->b;
        int64_t c;

        // Link to the next
        NeighborTriangle* secondNGR = firstNGR->next;

        // If there is neighbour
        if (secondNGR == nullptr)
            secondNGR = _neighborList[index];

        // Connect
        c = secondNGR->b;

        // Vertex position
        Vector3f position = getPositionSurfaceOnly(x, y, z, b, a, c);

        // Compute the angle between the triangle vectors
        float angle = computeDotProduct(b, a, c);

        // TODO: Understand why do we have to increment by 1.0
        angle += 1.0;

        // Update the vertex position
        nx += angle * position.x();
        ny += angle * position.y();
        nz += angle * position.z();

        // Update the weight
        weight += angle;

        // Connect
        firstNGR = firstNGR->next;
    }

    // If the weight is greater than zero
    if (weight > 0.0)
    {
        // Update the vertex position by the weight
        nx /= weight;
        ny /= weight;
        nz /= weight;

        // Compute the eigen vector
        Vector3f eigenValue;
        float maxAngle;
        EigenVector eigenVector =
                computeEigenVector(index, &eigenValue, &maxAngle);

        if (!eigenVector.isValid)
        {
            // Compute the new vertex position
            _vertices[index].x()= nx;
            _vertices[index].y()= ny;
            _vertices[index].z()= nz;
        }
        else
        {
            // Shift the new point
            nx -= x;
            ny -= y;
            nz -= z;

            // Compute XYZ weights
            const float w1 = (nx * eigenVector.x1 +
                              ny * eigenVector.y1 +
                              nz * eigenVector.z1) / (1.0f + eigenValue.x());
            const float w2 = (nx * eigenVector.x2 +
                              ny * eigenVector.y2 +
                              nz * eigenVector.z2) / (1.0 + eigenValue.y());
            const float w3 = (nx * eigenVector.x3 +
                              ny * eigenVector.y3 +
                              nz * eigenVector.z3) / (1.0 + eigenValue.z());

            // Compute the new values
            nx = w1 * eigenVector.x1 + w2 * eigenVector.x2 + w3 * eigenVector.x3 + x;
            ny = w1 * eigenVector.y1 + w2 * eigenVector.y2 + w3 * eigenVector.y3 + y;
            nz = w1 * eigenVector.z1 + w2 * eigenVector.z2 + w3 * eigenVector.z3 + z;

            // Update the vertex position
            _vertices[index].x()= nx;
            _vertices[index].y()= ny;
            _vertices[index].z()= nz;
        }
    }
}

void Mesh::smoothNormal(const int64_t n)
{
    // The updated position of the vertex after smoothing the normal
    Vector3f position;

    // Counter
    int64_t number = 0;

    NeighborTriangle* firstNGR = _neighborList[n];
    while (firstNGR != nullptr)
    {
        int64_t a = firstNGR->a;
        int64_t b = firstNGR->b;

        NeighborTriangle* secondNGR = firstNGR->next;
        if (secondNGR == nullptr)
            secondNGR = _neighborList[n];
        int64_t c = secondNGR->b;

        NeighborTriangle* thirdNGR = secondNGR->next;
        if (thirdNGR == nullptr)
            thirdNGR = _neighborList[n];
        int64_t d = thirdNGR->b;

        NeighborTriangle* auxNGR = _neighborList[b];
        while (auxNGR != nullptr)
        {
            if ((auxNGR->a == c && auxNGR->b != n) ||
                    (auxNGR->b == c && auxNGR->a != n))
                break;
            auxNGR = auxNGR->next;
        }

        int64_t e;
        if (auxNGR->a == c && auxNGR->b != n)
            e = auxNGR->b;
        else if (auxNGR->b == c && auxNGR->a != n)
            e = auxNGR->a;
        else
            LOG_WARNING("Normal computation error");

        // Computing the cross product
        Vector3f normal = computeCrossProduct(n, b, c);

        // Store the initially computed normal
        Vector3f iNormal = normal;

        float dx = 0.0;
        float dy = 0.0;
        float dz = 0.0;

        int64_t num = 0;
        normal = computeCrossProduct(n, a, b);

        float length = normal.x() * iNormal.x() +
                normal.y() * iNormal.y() +
                normal.z() * iNormal.z();
        if (length > 0.0)
        {
            num++;
            dx += length * normal.x();
            dy += length * normal.y();
            dz += length * normal.z();
        }

        normal = computeCrossProduct(n, c, d);
        length = normal.x() * iNormal.x() +
                normal.y() * iNormal.y() +
                normal.z() *  iNormal.z();
        if (length > 0.0)
        {
            num++;
            dx += length * normal.x();
            dy += length * normal.y();
            dz += length * normal.z();
        }

        normal = computeCrossProduct(b, e, c);
        length = normal.x() * iNormal.x() +
                normal.y() * iNormal.y() +
                normal.z() * iNormal.z();
        if (length > 0.0)
        {
            num++;
            dx += length * normal.x();
            dy += length * normal.y();
            dz += length * normal.z();
        }

        length = sqrt(dx * dx + dy * dy + dz * dz);
        if (length > 0.0)
        {
            dx /= length;
            dy /= length;
            dz /= length;

            const float fx = iNormal.y() * dz - iNormal.z() * dy;
            const float fy = iNormal.z() * dx - iNormal.x() * dz;
            const float fz = iNormal.x() * dy - iNormal.y() * dx;

            const float cx = _vertices[c].x();
            const float cy = _vertices[c].y();
            const float cz = _vertices[c].z();

            const float bx = _vertices[b].x();
            const float by = _vertices[b].y();
            const float bz = _vertices[b].z();

            float theta, phi;
            length = fx * (bx - cx) + fy * (by - cy) + fz * (bz - cz);
            if (length >= 0)
            {
                theta = D2F(atan2(by - cy, bx - cx));
                phi = D2F(atan2(bz - cz,
                                sqrt((bx - cx) * (bx - cx) +
                                     (by - cy) * (by - cy))));
            }
            else
            {
                theta = D2F(atan2(cy - by, cx - bx));
                phi = D2F(atan2(cz - bz, sqrt((bx - cx) * (bx - cx) +
                                              (by - cy)*(by - cy))));
            }

            float alpha = acos(dx * iNormal.x() +
                               dy * iNormal.y() +
                               dz * iNormal.z()) / D2F(4.0 - num);

            Vector3f sv = rotate(_vertices[n].x() - cx,
                                 _vertices[n].y() - cy,
                                 _vertices[n].z() - cz,
                                 theta, phi, alpha);

            position += sv + Vector3f(cx, cy, cz);
            number++;
        }

        firstNGR = firstNGR->next;
    }

    if (number > 0 && !position.isNan())
    {
        _vertices[n] = position / I2F(number);
    }

}

void Mesh::smoothNormals()
{
    // Starting timer
    TIMER_SET;

    LOG_STATUS("Smoothing Mesh Normals");

    // Check if neighborlist is created
    if (_neighborList == nullptr)
        _createNeighbourList();

    // Normal smooth all vertices
    LOOP_STARTS("Smoothing Normals");
    for (size_t n = 0; n < _numberVertices; n++)
    {
        LOOP_PROGRESS(n, _numberVertices);

        smoothNormal(n);
    }
    LOOP_DONE;

    float minAngle, maxAngle;
    int64_t numSmall, numLarge;
    computeAngles(&minAngle, &maxAngle, &numSmall, &numLarge, 15, 150);

    LOG_WARNING("Min. Angle: [%f],  Max. Angle : [%f], Smaller than 15: [%d], Larger than 150: [%d]",
                F2D(minAngle), F2D(maxAngle), numSmall, numLarge);

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

bool Mesh::smooth(const int64_t &maxMinAngle, const int64_t &minMaxAngle,
                  const int64_t &maxIteration, const bool &flipEdges)
{
    // Starting the timer
    TIMER_SET;

    LOG_STATUS("Smoothing Mesh Surface");

    float minAngle, maxAngle;
    int64_t numSmall, numLarge;
    bool smoothed = false;

    // Check if neighborlist is created
    if (_neighborList == nullptr)
        _createNeighbourList();

    // Global index
    int64_t i = 0;

    // Print the initial quality only when doing 1 or more iterations
    if (maxIteration > 1)
    {
        computeAngles(&minAngle, &maxAngle, &numSmall,
                      &numLarge, maxMinAngle, minMaxAngle);

        LOG_DEBUG("[%2d]: Min. angle: %f, Max. angle: %f, "
                  "Smaller than %d: %d, Larger than %d: %d\n",
                  i, F2D(minAngle), F2D(maxAngle), maxMinAngle, numSmall,
                  minMaxAngle, numLarge);
    }

    LOOP_STARTS("Smoothing Faces");
    while (!smoothed && i < maxIteration)
    {
        LOOP_PROGRESS(i, maxIteration);

        // Increment the index
        ++i;

        // Smooth all vertices
        if (flipEdges)
        {
            for (size_t index = 0; index < _numberVertices; index++)
            {
                // Move the vertex along the surface and flip the edge
                moveVertexAlongSurface(index);
                edgeFlipping(index);
            }
        }
        else
        {
            for (size_t index = 0; index < _numberVertices; index++)
            {
                // Move only the vertex along the edge
                moveVertexAlongSurface(index);
            }
        }

        // Calculate and print quality after surtriangles smooth
        computeAngles(&minAngle, &maxAngle, &numSmall,
                      &numLarge, maxMinAngle, minMaxAngle);

        // Print the iteration number only when doing 1 or more iterations
        if (maxIteration != 1)
            LOG_DEBUG("[%2d ]: Min. angle: %f, Max. angle: %f, "
                      "Smaler than %d: %d, Larger than %d: %d \n",
                      i, F2D(minAngle), F2D(maxAngle), maxMinAngle, numSmall,
                      minMaxAngle, numLarge);
        else
        {
            LOG_DEBUG("Min. angle: %f - max_angle: %f - "
                      "Smaler than %d: %d, Larger than %d: %d \n",
                      F2D(minAngle), F2D(maxAngle), maxMinAngle, numSmall,
                      minMaxAngle, numLarge);
        }

        // Check if the mesh is smoothed or not
        smoothed = (minAngle > maxMinAngle) && (maxAngle < minMaxAngle);
    }
    LOOP_DONE;

    LOG_WARNING("Min. angle: %f, Max. angle: %f, "
                "Smaler than %d: %d, Larger than %d: %d \n",
                F2D(minAngle), F2D(maxAngle), maxMinAngle, numSmall, minMaxAngle, numLarge);

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    return smoothed;
}

void Mesh::subdivideTriangleAtCentroid(const size_t &triangleIndex,
                                       std::vector< Vector3f >& vertexList,
                                       std::vector< Triangle >& triangleList)
{
    // Get the triangle
    Triangle t = _triangles[triangleIndex];

    // Geth the vertices of the triangles
    const Vector3f& v0 = _vertices[t[0]];
    const Vector3f& v1 = _vertices[t[1]];
    const Vector3f& v2 = _vertices[t[2]];

    // Compute the centroid
    const Vector3f centroid = (v0 + v1 + v2) / 3.0;

    // Add the centroid to the vertex list
    vertexList.push_back(centroid);

    // Get the new vertex index
    const size_t vertexIndex = _numberVertices + vertexList.size() - 1;

    // Create three new triangles
    Triangle t0, t1, t2;
    t0[0] = t[0]; t0[1] = t[1]; t0[2] = vertexIndex;
    t1[0] = t[1]; t1[1] = t[2]; t1[2] = vertexIndex;
    t2[0] = t[2]; t2[1] = t[0]; t2[2] = vertexIndex;

    // Add the triangles to the list
    triangleList.push_back(t0);
    triangleList.push_back(t1);
    triangleList.push_back(t2);
}

void Mesh::subdivideTriangleAtMidPoints(const size_t &triangleIndex,
                                        std::vector< Vector3f >& vertexList,
                                        std::vector< Triangle >& triangleList)
{
    // Get the triangle
    const Triangle& t = _triangles[triangleIndex];

    const size_t& v0Index = t[0];
    const size_t& v1Index = t[1];
    const size_t& v2Index = t[2];

    // Geth the vertices of the triangles
    const Vector3f& v0 = _vertices[v0Index];
    const Vector3f& v1 = _vertices[v1Index];
    const Vector3f& v2 = _vertices[v2Index];

    // Compute the midpoints along each edge
    const Vector3f v3 = 0.5f * (v0 + v1); // v01
    const Vector3f v4 = 0.5f * (v1 + v2); // v12
    const Vector3f v5 = 0.5f * (v2 + v0); // v20

    // Add the midpoint vertices the vertex list
    vertexList.push_back(v3);
    vertexList.push_back(v4);
    vertexList.push_back(v5);

    // Get the new vertex index
    const size_t v3Index = _numberVertices + vertexList.size() - 3;
    const size_t v4Index = _numberVertices + vertexList.size() - 2;
    const size_t v5Index = _numberVertices + vertexList.size() - 1;

    // Create three new triangles
    Triangle t0, t1, t2, t3;

    // v0 v3 v5
    t0[0] = v0Index;
    t0[1] = v3Index;
    t0[2] = v5Index;

    // v3 v1 v4
    t1[0] = v3Index;
    t1[1] = v1Index;
    t1[2] = v4Index;

    // v3 v4 v5
    t2[0] = v3Index;
    t2[1] = v4Index;
    t2[2] = v5Index;

    // v4 v2 v5
    t3[0] = v4Index;
    t3[1] = v2Index;
    t3[2] = v5Index ;

    // Add the triangles to the list
    triangleList.push_back(t0);
    triangleList.push_back(t1);
    triangleList.push_back(t2);
    triangleList.push_back(t3);
}

void Mesh::subdivideTrianglseAtMidPoints()
{
    // New lists of vertices and triangles
    std::vector< Vector3f > createdVertices;
    std::vector< Triangle > createdTriangles;

    // Subdivide all the triangles
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        subdivideTriangleAtMidPoints(i, createdVertices, createdTriangles);
    }

    // Update the _vertices and _triangles lists
    const size_t totalNumberVertices = _numberVertices + createdVertices.size();
    const size_t totalNumberTriangles = createdTriangles.size();

    // Allocate the new arrays
    Vector3f* newVertices = new Vector3f[totalNumberVertices];
    Triangle* newTriangles = new Triangle[totalNumberTriangles];

    // Copy the old vertices to the new vertices list
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        newVertices[i] = _vertices[i];
    }

    // Release the old vertices list
    delete [] _vertices; _vertices = nullptr;

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < createdVertices.size(); ++i)
    {
        newVertices[_numberVertices + i] = createdVertices[i];
    }

    createdVertices.clear();
    createdVertices.shrink_to_fit();

    delete [] _triangles; _triangles = nullptr;

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < createdTriangles.size(); ++i)
    {
        newTriangles[i] = createdTriangles[i];
    }

    createdTriangles.clear();
    createdTriangles.shrink_to_fit();

    _numberVertices = totalNumberVertices;
    _numberTriangles = totalNumberTriangles;

    _vertices = newVertices;
    _triangles = newTriangles;
}

void Mesh::refineSelectedTriangles(const std::vector< size_t > &trianglesIndices)
{
    // Starting the timer
    TIMER_SET;

    LOG_STATUS("Refining Triangles");

    // Destroy the neighbouring list becuase the number of vertices is going to change later
    _destroyNeighborlist();

    // Create the triangle markers
    std::vector< bool > triangleMarkers;
    triangleMarkers.resize(_numberTriangles);

    // Initialize them to false
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < triangleMarkers.size(); ++i)
        triangleMarkers[i] = false;

    // Create the vectors
    std::vector< Vector3f > vertices;
    std::vector< Triangle > triangles;

    // Keeps track on the number of divided triangles
    size_t numberDividedTriangles = 0;

    // For the moment, it refines all the triangles, but we will use it to only refine the
    // selected triangles.
    LOOP_STARTS("Refine Selected Triangles")
    LOOP_PROGRESS(0, 5);
    for (size_t i = 0; i < trianglesIndices.size(); ++i)
    {
        // If the triangle has been visited before, then skip it
        if (triangleMarkers[trianglesIndices[i]])
            continue;

        // Subdivide at the centroid
        subdivideTriangleAtCentroid(trianglesIndices[i], vertices, triangles);

        // Update the count
        ++numberDividedTriangles;

        // Update the marker
        triangleMarkers[trianglesIndices[i]] = true;
    }

    // Compute the total number of vertices in the subdivided mesh
    const size_t totalNumberVertices = _numberVertices + vertices.size();

    // Compute the total number of MANIFOLD triangles in the subdivided mesh
    const size_t numberManifoldTriangles = _numberTriangles - numberDividedTriangles;
    const size_t totalNumberTriangles = numberManifoldTriangles + triangles.size();

    // Allocate the new arrays
    Vector3f* newVertices = new Vector3f[totalNumberVertices];
    Triangle* newTriangles = new Triangle[totalNumberTriangles];

    LOOP_PROGRESS(1, 5);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        newVertices[i] = _vertices[i];
    }

    delete [] _vertices; _vertices = nullptr;

    LOOP_PROGRESS(2, 5);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        newVertices[_numberVertices + i] = vertices[i];
    }

    vertices.clear();
    vertices.shrink_to_fit();

    LOOP_PROGRESS(3, 5);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        newTriangles[i] = triangles[i];
    }

    LOOP_PROGRESS(4, 5);
    size_t tIndex = 0, nIndex = 0;
    while (true)
    {
        tIndex++;

        // If we reach the maximum break
        if (tIndex > _numberTriangles)
            break;

        // If the triangles is in the markers list, then it will not be added
        if (triangleMarkers[tIndex - 1])
            continue;
        else
        {
            newTriangles[triangles.size() + nIndex] = _triangles[tIndex - 1];
            nIndex++;
        }

    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    delete [] _triangles; _triangles = nullptr;

    triangles.clear();
    triangles.shrink_to_fit();

    _numberVertices = totalNumberVertices;
    _numberTriangles = totalNumberTriangles;

    _vertices = newVertices;
    _triangles = newTriangles;
}

void Mesh::refineROIs(const ROIs& regions)
{
    // Starting the timer
    TIMER_SET;

    LOG_STATUS("Refining Mesh Surface - ROI");

    // Selected triangles
    std::vector< size_t > trianglesIndices;

    // Select the triangles that are located within the ROI
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        // Get the triangle
        const auto& triangle = _triangles[i];

        // Get the vertices of the triangles
        const auto v0 = _vertices[triangle[0]];
        const auto v1 = _vertices[triangle[1]];
        const auto v2 = _vertices[triangle[2]];

        for (size_t j = 0; j < regions.size(); ++j)
        {
            const auto region = regions[j];
            if (isTriangleInSphere(v0, v1, v2, region->center, region->radius))
            {
                trianglesIndices.push_back(i);
            }
        }
    }

    // The triangles are already selected, now time to refine them
    refineSelectedTriangles(trianglesIndices);

    // Clear the selected triangles list
    trianglesIndices.clear();
    trianglesIndices.shrink_to_fit();

    // Statistics
    LOG_STATUS("Refinement");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::kdTreeMapping(std::vector< Vector3f >& pointCloud, const bool& showProgress)
{
    auto kdtree = KdTree::from(pointCloud);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        auto nearestPoint = kdtree.findNearestPoint(_vertices[i]);
        _vertices[i].x() = nearestPoint.position.x();
        _vertices[i].y() = nearestPoint.position.y();
        _vertices[i].z() = nearestPoint.position.z();
    }
}

void Mesh::kdTreeMapping(const KdTree& kdTree, const bool& showProgress)
{
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        auto nearestPoint = kdTree.findNearestPoint(_vertices[i]);
        _vertices[i].x() = nearestPoint.position.x();
        _vertices[i].y() = nearestPoint.position.y();
        _vertices[i].z() = nearestPoint.position.z();
    }
}

void Mesh::map(std::vector< Vector3f >& pointCloud, const bool& showProgress)
{
    if (showProgress)
    {
        // Starting the timer
        TIMER_SET;

        LOG_STATUS("Mapping Vertices");

        LOOP_STARTS("Mapping Loop");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
                for (size_t i = 0; i < _numberVertices; ++i)
        {
            PROGRESS_UPDATE;
            LOOP_PROGRESS(PROGRESS, _numberVertices);

            float minDistance = std::numeric_limits<float>::max();
            int64_t minIndex = -1;

            for (size_t j = 0; j < pointCloud.size(); ++j)
            {
                float distance = _vertices[i].distance(pointCloud[j]);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    minIndex = j;
                }
            }

            _vertices[i] = pointCloud[minIndex];
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);

        LOG_STATUS("Mapping");
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for (size_t i = 0; i < _numberVertices; ++i)
        {
            float minDistance = std::numeric_limits<float>::max();
            int64_t minIndex = -1;

            for (size_t j = 0; j < pointCloud.size(); ++j)
            {
                float distance = _vertices[i].distance(pointCloud[j]);

                if (distance < minDistance)
                {
                    minDistance = distance;
                    minIndex = j;
                }
            }

            _vertices[i] = pointCloud[minIndex];
        }
    }
}

void Mesh::map(Mesh* toMesh)
{
    // Starting the timer
    TIMER_SET;

    LOG_STATUS("Mapping Vertices");

    LOOP_STARTS("Mapping Loop");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        PROGRESS_UPDATE;
        LOOP_PROGRESS(PROGRESS, _numberVertices);

        float minDistance = std::numeric_limits<float>::max();
        int64_t minIndex = -1;

        for (size_t j = 0; j < toMesh->getNumberVertices(); ++j)
        {
            float distance = _vertices[i].distance(toMesh->getVertices()[j]);

            if (distance < minDistance)
            {
                minDistance = distance;
                minIndex = j;
            }
        }

        _vertices[i].x() = toMesh->getVertices()[minIndex].x();
        _vertices[i].y() = toMesh->getVertices()[minIndex].y();
        _vertices[i].z() = toMesh->getVertices()[minIndex].z();

    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOG_STATUS("Mapping");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::refine()
{
    // Starting the timer
    TIMER_SET;
    LOG_STATUS("Refining Mesh Surface");

    // Check if neighborlist is created
    if (_neighborList == nullptr)
        _createNeighbourList();

    // Create an array with the number of edges associated with each vertex
    int64_t* numEdges = new int64_t[_numberVertices];

    // Create an array with the offsets into the vertex2edge array for each vertex
    int64_t* offsets = new int64_t[_numberVertices];

    // Iterate over all vertices and collect edges
    int64_t totalNumberEdges = 0, localNumberEdges = 0;
    for (int64_t i = 0; i < UI2I64(_numberVertices); ++i)
    {
        // Offsets
        offsets[i] = totalNumberEdges;
        localNumberEdges = 0;

        NeighborTriangle* ngr = _neighborList[i];
        while (ngr != nullptr)
        {
            // If i is smaller than ngr->a we have an edge
            if (i < ngr->a)
            {
                totalNumberEdges++;
                localNumberEdges++;
            }

            ngr = ngr->next;
        }

        numEdges[i] = localNumberEdges;
    }

    // Create the refined mesh structure
    Mesh* refinedMesh = instanciate(_numberVertices + totalNumberEdges, _numberTriangles * 4);

    // For the moment, just keep the old numbers.
    refinedMesh->_numberVertices = _numberVertices;
    refinedMesh->_numberTriangles = _numberTriangles;

    int64_t edgeNumber = 0;

    // Copy the vertices of the original mesh to the new mesh
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        refinedMesh->_vertices[i].x()= _vertices[i].x();
        refinedMesh->_vertices[i].y()= _vertices[i].y();
        refinedMesh->_vertices[i].z()= _vertices[i].z();
    }

    // Copy the triangles of the original mesh to the new mesh
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        refinedMesh->_triangles[i]= _triangles[i];
    }

    // Create the map from vertices to edges
    int64_t* vertex2edge = new int64_t[totalNumberEdges];

    // Iterate over all vertices and split edges
    LOOP_STARTS("Splitting Edges");
    for (int64_t i = 0; i < UI2I64(_numberVertices); ++i)
    {
        LOOP_PROGRESS(i, _numberVertices);

        // Get the coordinates of vertex i
        const float nx = refinedMesh->_vertices[i].x();
        const float ny = refinedMesh->_vertices[i].y();
        const float nz = refinedMesh->_vertices[i].z();

        // Get the neighbours list at vertex i
        NeighborTriangle* ngr = _neighborList[i];
        while (ngr != nullptr)
        {
            // If n is smaller than ngr->a we have an edge
            if (i < ngr->a)
            {
                // Add the value of the opposite vertex to the map
                vertex2edge[edgeNumber] = ngr->a;

                // Get the coordinates of vertex ngr->a
                const float ax = refinedMesh->_vertices[ngr->a].x();
                const float ay = refinedMesh->_vertices[ngr->a].y();
                const float az = refinedMesh->_vertices[ngr->a].z();

                // Add the new vertex coordinates of the splitted edge
                refinedMesh->_vertices[_numberVertices + edgeNumber].x() = 0.5f * (ax + nx);
                refinedMesh->_vertices[_numberVertices + edgeNumber].y() = 0.5f * (ay + ny);
                refinedMesh->_vertices[_numberVertices + edgeNumber].z() = 0.5f * (az + nz);

                // Increase the edge number
                edgeNumber++;
            }

            ngr = ngr->next;
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // A counter for adding new faces or triangles
    int64_t triangleNumber = refinedMesh->_numberTriangles;

    // Iterate over faces and add information of the refined face
    int64_t localVertices[3], additionalLocalVertices[3];

    for (size_t i = 0; i < refinedMesh->_numberTriangles; ++i)
    {
        localVertices[0] = refinedMesh->_triangles[i][0];
        localVertices[1] = refinedMesh->_triangles[i][1];
        localVertices[2] = refinedMesh->_triangles[i][2];

        // Iterate over the vertices and find the edges
        for (int64_t j = 0; j < 3; j++)
        {
            int64_t minVertexNum =
                    std::min(localVertices[j], localVertices[(j + 1) % 3]);
            int64_t maxVertexNum =
                    std::max(localVertices[j], localVertices[(j + 1) % 3]);

            // Find the edge number that fit the pair of vertices
            int64_t k = 0;
            for (k = 0; k < numEdges[minVertexNum]; k++)
            {
                if (vertex2edge[offsets[minVertexNum] + k] == maxVertexNum)
                    break;
            }

            // The edge number represents the number of the added vertex plus
            // the number of original vertices
            additionalLocalVertices[j] = _numberVertices +
                    offsets[minVertexNum] + k;
        }

        // Add information of the four new faces

        // First the mid face
        refinedMesh->_triangles[i][0] = additionalLocalVertices[0];
        refinedMesh->_triangles[i][1] = additionalLocalVertices[1];
        refinedMesh->_triangles[i][2] = additionalLocalVertices[2];

        // Then the three corner faces
        for (int64_t m = 0; m < 3; m++)
        {
            refinedMesh->_triangles[triangleNumber][0] = localVertices[m];
            refinedMesh->_triangles[triangleNumber][1] = additionalLocalVertices[m];
            refinedMesh->_triangles[triangleNumber][2] = additionalLocalVertices[(m + 2) % 3];
            triangleNumber++;
        }
    }

    // Release memory
    delete [] numEdges;
    delete [] offsets;
    delete [] vertex2edge;

    // Update number information
    refinedMesh->_numberVertices += totalNumberEdges;
    refinedMesh->_numberTriangles *= 4;

    // Release old data
    _releaseData();

    // Assign the refined mesh to the passed
    _numberVertices = refinedMesh->_numberVertices;
    _numberTriangles = refinedMesh->_numberTriangles;
    _vertices = refinedMesh->_vertices;
    _triangles = refinedMesh->_triangles;

    // Recreate the neigborlist
    _createNeighbourList();

    // Statistics
    LOG_STATUS("Refinement");
    LOG_STATS(GET_TIME_SECONDS);
}

bool Mesh::coarse(const float& coarseRate,
                  const float& flatnessRate,
                  const float& densenessWeight,
                  const float& maxNormalAngle,
                  const int64_t &iteration)
{
    // If the mesh has less than 1k polygons, then don't coarse it
    if (_numberTriangles < 1024)
        return false;

    // Starting the timer
    TIMER_SET;

    LOG_STATUS("Coarsing Mesh [%d]", iteration + 1);

    const size_t initialNumberVertices = _numberVertices;
    const size_t initialNumberTriangles = _numberTriangles;

    // Check if neighborlist is created, otherwise create it
    if (_neighborList == nullptr)
    {
        _createNeighbourList();
    }

    // Allocate the vertex and triangle index arrays
    int64_t* vertexIndex = new int64_t[_numberVertices];
    int64_t* triangleIndex = new int64_t[_numberTriangles];

    // If using sparseness weight, calculate the average segment length of mesh
    float averageEdgeLength = 0.0;
    if (densenessWeight > 0.0)
    {
        // Average edge length of the triangle
        float averageLength = 0.0;

        LOOP_STARTS("Computing Edges");
        TIMER_RESET;
        for (size_t i = 0; i < _numberTriangles; ++i)
        {
            LOOP_PROGRESS(i, _numberTriangles);

            // Indices
            Triangle& t = _triangles[i];
            const float nx = sqrt((_vertices[t[0]].x()- _vertices[t[1]].x()) *
                                  (_vertices[t[0]].x()- _vertices[t[1]].x()) +
                                  (_vertices[t[0]].y()- _vertices[t[1]].y()) *
                                  (_vertices[t[0]].y()- _vertices[t[1]].y()) +
                                  (_vertices[t[0]].z()- _vertices[t[1]].z()) *
                                  (_vertices[t[0]].z()- _vertices[t[1]].z()));
            const float ny = sqrt((_vertices[t[0]].x()- _vertices[t[2]].x()) *
                                  (_vertices[t[0]].x()- _vertices[t[2]].x()) +
                                  (_vertices[t[0]].y()- _vertices[t[2]].y()) *
                                  (_vertices[t[0]].y()- _vertices[t[2]].y()) +
                                  (_vertices[t[0]].z()- _vertices[t[2]].z()) *
                                  (_vertices[t[0]].z()- _vertices[t[2]].z()));
            const float nz = sqrt((_vertices[t[2]].x()- _vertices[t[1]].x()) *
                                  (_vertices[t[2]].x()- _vertices[t[1]].x()) +
                                  (_vertices[t[2]].y()- _vertices[t[1]].y()) *
                                  (_vertices[t[2]].y()- _vertices[t[1]].y()) +
                                  (_vertices[t[2]].z()- _vertices[t[1]].z()) *
                                  (_vertices[t[2]].z()- _vertices[t[1]].z()));

            averageLength += ((nx + ny + nz) * 0.33333333334);
        }
        LOOP_DONE;

        // Statistics
        LOG_STATS(GET_TIME_SECONDS);

        if (_numberTriangles == 0)
        {
            LOG_ERROR("Zero degree on a vertex !");
        }
        else
        {
            averageEdgeLength = averageLength / (1.0 * _numberTriangles);
        }
    }

    // The main loop over all vertices
    int64_t faceAvailableList[64];
    int64_t faceAvailableIndex;
    int64_t neighborAuxList[64];

    // Start the algorithm with the current number of vertices in the mesh
    // This number will be reduced later if any coarsening happens
    TIMER_RESET;
    LOOP_STARTS("Removing Vertices");
    size_t vertexNumber = _numberVertices;
    for (int64_t n = 0; n < UI2I64(_numberVertices); ++n)
    {
        LOOP_PROGRESS(n, _numberVertices);

        // If the vertex is labeled in the vertex markers, then ignore it
        if (_vertexMarkers.size() > 0)
        {
            if (_vertexMarkers[n] == 1)
            {
                continue;
            }
        }

        // Check if the vertex has enough neigborgs to be deleted
        char deleteFlag = 1;
        NeighborTriangle* firstNGR = _neighborList[n];

        float maxLength;
        while (firstNGR != nullptr)
        {
            const int64_t a = firstNGR->a;
            int64_t number = 0;
            int64_t someNumber = 0;

            NeighborTriangle* secondNGR = _neighborList[a];
            while (secondNGR != nullptr)
            {
                const int64_t b = secondNGR->a;
                NeighborTriangle* auxNGR = _neighborList[n];

                while (auxNGR != nullptr)
                {
                    if (auxNGR ->a == b)
                        someNumber++;
                    auxNGR  = auxNGR ->next;
                }

                number++;
                secondNGR = secondNGR->next;
            }

            if (number <= 3 || someNumber > 2)
                deleteFlag = 0;

            firstNGR = firstNGR->next;
        }

        if (deleteFlag)
        {
            const float x = _vertices[n].x();
            const float y = _vertices[n].y();
            const float z = _vertices[n].z();

            maxLength = -1;

            firstNGR = _neighborList[n];

            // If using sparseness as a criteria for coarsening
            // calculate the maximal segment length
            float ratio2 = 1.0;
            if (densenessWeight > 0.0)
            {
                while (firstNGR != nullptr)
                {
                    const int64_t a = firstNGR->a;
                    const int64_t b = firstNGR->b;

                    const float nx = sqrt((x - _vertices[a].x()) *
                                          (x - _vertices[a].x()) +
                                          (y - _vertices[a].y()) *
                                          (y - _vertices[a].y()) +
                                          (z - _vertices[a].z()) *
                                          (z - _vertices[a].z()));
                    const float ny = sqrt((x - _vertices[b].x()) *
                                          (x - _vertices[b].x()) +
                                          (y - _vertices[b].y()) *
                                          (y - _vertices[b].y()) +
                                          (z - _vertices[b].z()) *
                                          (z - _vertices[b].z()));
                    if (nx > maxLength)
                    {
                        maxLength = nx;
                    }

                    if (ny > maxLength)
                    {
                        maxLength = ny;
                    }

                    firstNGR = firstNGR->next;
                }

                // Max segment length over the average segment length of the mesh
                ratio2 = maxLength / averageEdgeLength;
                ratio2 = pow(ratio2, densenessWeight);
            }

            // If using curvatory as a coarsening criteria
            // calculate the local structure tensor
            float ratio1 = 1.0;
            float maxAngle = 0.0;
            if (flatnessRate > 0.0)
            {
                Vector3f eigenValue;
                EigenVector eigenVector = computeEigenVector(n, &eigenValue, &maxAngle);

                if (!eigenVector.isValid)
                {
                    ratio1 = 999999.0f;
                }
                else
                {
                    if (isZero(eigenValue.x()))
                    {
                        LOG_ERROR("ERROR: Maximum eigen value is zero! \n");
                    }
                    else

                    {
                        ratio1 = std::fabs((eigenValue.y())/(eigenValue.x()));
                        ratio1 = pow(ratio1, flatnessRate);
                    }
                }
            }

            // Compare the two coarseness criterias against the given
            // coarsing rate
            char deleteVertex = ratio1 * ratio2 < coarseRate;

            // Use maximal angle between vertex normal as a complementary
            // coarse criteria
            if (maxNormalAngle > 0)
                deleteVertex = deleteVertex && maxAngle > maxNormalAngle;

            // Deleting a vertex and retrianglulate the hole
            if (deleteVertex)
            {
                vertexNumber--;

                // Delete the vertex by labeling it here and removing it later.
                _vertices[n].x()= VERTEX_DELETION_VALUE;
                _vertices[n].y()= VERTEX_DELETION_VALUE;
                _vertices[n].z()= VERTEX_DELETION_VALUE;

                int64_t neighborNumber = 0;
                firstNGR = _neighborList[n];

                while (firstNGR != nullptr)
                {
                    const int64_t a = firstNGR->a;
                    const int64_t c = firstNGR->c;
                    faceAvailableList[neighborNumber] = c;
                    neighborAuxList[neighborNumber] = a;
                    neighborNumber++;

                    // Delete faces associated with vertex n
                    _triangles[c][0] = -1;
                    _triangles[c][1] = -1;
                    _triangles[c][2] = -1;

                    /* delete neighbors associated with vertex n */
                    NeighborTriangle* secondNGR = _neighborList[a];
                    NeighborTriangle* auxNGR  = secondNGR;

                    while (secondNGR != nullptr)
                    {
                        if (secondNGR->a == n || secondNGR->b == n)
                        {
                            if (secondNGR == _neighborList[a]) {
                                _neighborList[a] = secondNGR->next;
                                delete secondNGR;
                                secondNGR = _neighborList[a];
                                auxNGR  = secondNGR;
                            }
                            else
                            {
                                auxNGR ->next = secondNGR->next;
                                delete secondNGR;
                                secondNGR = auxNGR ->next;
                            }
                        }
                        else {
                            if (secondNGR == _neighborList[a])
                            {
                                secondNGR = secondNGR->next;
                            }
                            else
                            {
                                auxNGR  = secondNGR;
                                secondNGR = secondNGR->next;
                            }
                        }
                    }

                    int64_t number = 0;
                    secondNGR = _neighborList[a];
                    while (secondNGR != nullptr)
                    {
                        number++;
                        secondNGR = secondNGR->next;
                    }

                    firstNGR->b = number;
                    firstNGR = firstNGR->next;
                }

                firstNGR = _neighborList[n];

                while (firstNGR->next != nullptr)
                    firstNGR = firstNGR->next;

                firstNGR->next = _neighborList[n];

                faceAvailableIndex = 0;

                subdividePolygin(_neighborList[n], faceAvailableList,
                                 &faceAvailableIndex);

                // Ordering the neibghours
                for (int64_t m = 0; m < neighborNumber; m++)
                {
                    firstNGR = _neighborList[neighborAuxList[m]];
                    const int64_t c = firstNGR->a;
                    while (firstNGR != nullptr)
                    {
                        const int64_t a = firstNGR->a;
                        const int64_t b = firstNGR->b;
                        NeighborTriangle* secondNGR = firstNGR->next;
                        while (secondNGR != nullptr)
                        {
                            const int64_t a0 = secondNGR->a;
                            const int64_t b0 = secondNGR->b;

                            // Assume counter clockwise orientation
                            if (a0==b && b0!=a)
                            {
                                NeighborTriangle* auxNGR  = firstNGR;
                                while (auxNGR  != nullptr)
                                {
                                    if (auxNGR ->next == secondNGR)
                                    {
                                        auxNGR ->next = secondNGR->next;
                                        break;
                                    }
                                    auxNGR  = auxNGR ->next;
                                }
                                auxNGR  = firstNGR->next;
                                firstNGR->next = secondNGR;
                                secondNGR->next = auxNGR ;
                                break;
                            }

                            secondNGR = secondNGR->next;
                        }
                        if (firstNGR->next == nullptr) {
                            if (firstNGR->b != c) {
                                LOG_ERROR("Some polygons are "
                                          "NOT closed: @[%d]", n);
                            }
                        }

                        firstNGR = firstNGR->next;
                    }
                }

                // Smooth the neighbors
                for (int64_t m = 0; m < neighborNumber; m++)
                {
                    const int64_t someNumber = neighborAuxList[m];
                    const float x = _vertices[someNumber].x();
                    const float y = _vertices[someNumber].y();
                    const float z = _vertices[someNumber].z();
                    float nx = 0.0;
                    float ny = 0.0;
                    float nz = 0.0;

                    float weight = 0.0;
                    firstNGR = _neighborList[someNumber];
                    while (firstNGR != nullptr)
                    {
                        const int64_t a = firstNGR->a;
                        const int64_t b = firstNGR->b;
                        NeighborTriangle* secondNGR = firstNGR->next;
                        if (secondNGR == nullptr)
                            secondNGR = _neighborList[someNumber];
                        const int64_t c = secondNGR->b;
                        Vector3f position = getPositionSurfaceOnly(x,y,z,b,a,c);
                        float angle = computeDotProduct(b, a, c);
                        angle += 1.0;
                        nx += angle*position.x();
                        ny += angle*position.y();
                        nz += angle*position.z();

                        weight += angle;
                        firstNGR = firstNGR->next;
                    }

                    if (weight > 0)
                    {
                        nx /= weight;
                        ny /= weight;
                        nz /= weight;

                        Vector3f eigenValue;
                        EigenVector eigenVector =
                                computeEigenVector(someNumber,
                                                   &eigenValue, &maxAngle);
                        if(!eigenVector.isValid)
                        {
                            _vertices[someNumber].x()= nx;
                            _vertices[someNumber].y()= ny;
                            _vertices[someNumber].z()= nz;
                        }
                        else
                        {
                            nx -= x; ny -= y; nz -= z;

                            const float w1 = (nx * eigenVector.x1 +
                                              ny * eigenVector.y1 +
                                              nz * eigenVector.z1) /
                                    (1.0  + eigenValue.x());
                            const float w2 = (nx * eigenVector.x2 +
                                              ny * eigenVector.y2 +
                                              nz * eigenVector.z2) /
                                    (1.0 + eigenValue.y());
                            const float w3 = (nx * eigenVector.x3 +
                                              ny * eigenVector.y3 +
                                              nz * eigenVector.z3) /
                                    (1.0 + eigenValue.z());

                            _vertices[someNumber].x()= w1 * eigenVector.x1 +
                                    w2 * eigenVector.x2 +
                                    w3 * eigenVector.x3 + x;
                            _vertices[someNumber].y()= w1 * eigenVector.y1 +
                                    w2 * eigenVector.y2 +
                                    w3 * eigenVector.y3 + y;
                            _vertices[someNumber].z()= w1 * eigenVector.z1 +
                                    w2 * eigenVector.z2 +
                                    w3 * eigenVector.z3 + z;
                        }
                    }
                }
            }
        }
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Clean the lists of nodes and faces
    int64_t startIndex = 0;
    for (int64_t n = 0; n < UI2I64(_numberVertices); ++n)
    {
        // If the vertex is valid
        if (isValidVertex(_vertices[n]))
        {
            if (startIndex != n)
            {
                _vertices[startIndex].x() = _vertices[n].x();
                _vertices[startIndex].y() = _vertices[n].y();
                _vertices[startIndex].z() = _vertices[n].z();
                _neighborList[startIndex] = _neighborList[n];
            }

            if (_vertexMarkers.size() > 1)
            {
                _vertexMarkers[startIndex] = _vertexMarkers[n];
            }

            vertexIndex[n] = startIndex;
            startIndex++;
        }

        // Otherwise, delete it
        else
        {
            vertexIndex[n] = -1;
        }
    }

    _numberVertices = startIndex;

    startIndex = 0;
    for (size_t n = 0; n < _numberTriangles; n++)
    {
        const int64_t a = _triangles[n][0];
        const int64_t b = _triangles[n][1];
        const int64_t c = _triangles[n][2];

        if (a >= 0 && vertexIndex[a] >= 0 &&
            b >= 0 && vertexIndex[b] >= 0 &&
            c >= 0 && vertexIndex[c] >= 0)
        {
            _triangles[startIndex][0] = vertexIndex[a];
            _triangles[startIndex][1] = vertexIndex[b];
            _triangles[startIndex][2] = vertexIndex[c];

            triangleIndex[n] = startIndex;
            startIndex++;
        }
        else
        {
            triangleIndex[n] = -1;
        }
    }
    _numberTriangles = startIndex;
    for (size_t i = 0; i < _numberVertices; ++i)
    {
        NeighborTriangle* firstNGR = _neighborList[i];

        while (firstNGR != nullptr)
        {
            const int64_t a = firstNGR->a;
            const int64_t b = firstNGR->b;
            const int64_t c = firstNGR->c;

            firstNGR->a = vertexIndex[a];
            firstNGR->b = vertexIndex[b];
            firstNGR->c = triangleIndex[c];

            firstNGR = firstNGR->next;
        }
    }

    size_t numberRemovedVertices = initialNumberVertices - _numberVertices;
    size_t numberRemovedTriangles = initialNumberTriangles - _numberTriangles;
    LOG_DETAIL("Removed Vertices: [ %s ], Removed Triangles: [ %s ]",
               FORMAT(numberRemovedVertices),
               FORMAT(numberRemovedTriangles));

    delete[] vertexIndex;
    delete[] triangleIndex;

    if (numberRemovedVertices > 0)
        return true;

    return false;
}

void Mesh::coarseDense(const float& denseRate, const int64_t &iterations)
{
    // Starting the timer
    TIMER_SET;

    for (int64_t i = 0; i < iterations; ++i)
        if (!coarse(denseRate, 0, 10, -1, i)) break;

    // Statistics
    LOG_STATUS_IMPORTANT("Dense Coarsing (Decimation) Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::coarseFlat(const float& flatnessRate,
                      const int64_t &iterations)
{
    // Starting the timer
    TIMER_SET;

    for (int64_t i = 0; i < iterations; ++i)
        if (!coarse(flatnessRate, 1, 0, -1, i)) break;

    // Statistics
    LOG_STATUS("Flat Coarsing (Decimation) Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::optimizeAdaptively(const size_t &optimizationIterations,
                              const size_t &smoothingIterations,
                              const float &flatFactor,
                              const float &denseFactor)
{
    LOG_TITLE("Adaptive Mesh Optimization");

    // Starting the timer
    TIMER_SET;

    // Coarse flat
    coarseFlat(flatFactor, optimizationIterations);

    // Smooth normals
    smoothNormals();

    // Smooth
    smooth(15, 150, smoothingIterations);

    // Smooth normals
    smoothNormals();

    // Statistics
    _optimizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Total Optimization");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::optimize(const size_t &optimizationIterations,
                    const int64_t &smoothingIterations,
                    const float& denseFactor)
{
    LOG_TITLE("Mesh Optimization");

    // Starting the timer
    TIMER_SET;

    // Remove the unnecessary vertices in multiple iterations
    coarseDense(denseFactor, optimizationIterations);

    // Smooth the normals
    smoothNormals();

    // Smoothing the surface
    smooth(15, 150, smoothingIterations);

    // Smooth the normals again
    smoothNormals();

    // Statistics
    _optimizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Total Optimization");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::optimizeUsingDefaultParameters()
{
    LOG_TITLE("Mesh Optimization");

    // Starting the timer
    TIMER_SET;

    // Refine
    refine();

    // Smooth
    smooth(15, 150, 4);

    // Coarse dense
    coarseDense(2.0, 4);

    // Smooth
    smooth(15, 150, 2);

    // Coarse dense
    coarseDense(2.0, 1);

    // Smooth
    smooth(20, 140, 6);

    // Normal smooth
    smoothNormals();

    // Coarse flat
    coarseFlat(0.035f);

    // Smooth
    smooth(20, 140, 6);

    // Coarse dense
    coarseDense(1.5, 4);

    // Smooth
    smooth(20, 140, 6);

    // Normal smooth
    smoothNormals();

    // Statistics
    _optimizationTime = GET_TIME_SECONDS;

    // Statistics
    LOG_STATUS_IMPORTANT("Total Optimization");
    LOG_STATS(GET_TIME_SECONDS);
}

void Mesh::removeFloatingFaces()
{
    // For all the triangles, if they have zero surface area, then they must
    // be removed
    Triangles triangles;
    std::vector< size_t > vertexIndex;
    for (size_t i = 0; i < _numberTriangles; ++i)
    {
        Triangle t = _triangles[i];
        const float triangleArea = computeTriangleSurfaceArea(
                    _vertices[t[0]], _vertices[t[1]], _vertices[t[2]]);

        if (triangleArea > 0.0)
            triangles.push_back(t);
        else
            vertexIndex.push_back(i);
    }

    // Delete the old triangles list
    delete _triangles;

    const size_t removedFaces = _numberTriangles - triangles.size();
    if (removedFaces > 0)
    {
        LOG_WARNING("Removing [%d] floating faces", removedFaces);

        // Create the new list
        _triangles = new Triangle[triangles.size()];

        for (size_t i = 0; i < triangles.size(); ++i)
            _triangles[i] = triangles[i];
        _numberTriangles = triangles.size();
    }
}

}
