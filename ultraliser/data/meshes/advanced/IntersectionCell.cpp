/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#include "IntersectionCell.h"
#include <algorithms/Sorting.h>
#include <utilities/Timer.h>
#include <data/meshes/advanced/Defines.hh>
#include <data/meshes/advanced/Utilities.h>

namespace Ultraliser
{

IntersectionCell::IntersectionCell(AdvancedMesh *input, bool useAll)
{
    Node *iNode;
    AdvancedVertex *vertex;
    AdvancedTriangle *triangle;

    // Initially
    pMax.x = -DBL_MAX, pMin.x = DBL_MAX;
    pMax.y = -DBL_MAX, pMin.y = DBL_MAX;
    pMax.z = -DBL_MAX, pMin.z = DBL_MAX;

    FOR_EACH_VV_VERTEX((&(input->_vertices)), vertex, iNode)
    {
        if (useAll || IS_BIT(vertex,5))
        {
            if (vertex->x < pMin.x) pMin.x = vertex->x;
            if (vertex->x > pMax.x) pMax.x = vertex->x;
            if (vertex->y < pMin.y) pMin.y = vertex->y;
            if (vertex->y > pMax.y) pMax.y = vertex->y;
            if (vertex->z < pMin.z) pMin.z = vertex->z;
            if (vertex->z > pMax.z) pMax.z = vertex->z;
        }
    }

    pMin -= DI_EPSILON_POINT;
    pMax += DI_EPSILON_POINT;

    FOR_EACH_VT_TRIANGLE((&(input->_triangles)), triangle, iNode)
    {
        if (useAll || IS_VISITED(triangle))
        {
            triangles.appendTail(triangle);
        }
    }
}

bool IntersectionCell::isTriangleBoundingBoxInCell(AdvancedTriangle *triangle) const
{
    AdvancedVertex *v1 = triangle->v1();
    AdvancedVertex *v2 = triangle->v2();
    AdvancedVertex *v3 = triangle->v3();

    const double minX = MIN(v1->x, MIN(v2->x, v3->x));
    const double maxX = MAX(v1->x, MAX(v2->x, v3->x));
    const double minY = MIN(v1->y, MIN(v2->y, v3->y));
    const double maxY = MAX(v1->y, MAX(v2->y, v3->y));
    const double minZ = MIN(v1->z, MIN(v2->z, v3->z));
    const double maxZ = MAX(v1->z, MAX(v2->z, v3->z));

    // Triangle BB is not entirely out of cell
    if (!(maxX < pMin.x || minX > pMax.x ||
          maxY < pMin.y || minY > pMax.y ||
          maxZ < pMin.z || minZ > pMax.z))
    {
        return true;
    }
    else
    {
        return false;
    }
}

IntersectionCell *IntersectionCell::fork()
{
    Node *iNode;
    AdvancedTriangle *triangle;
    AdvancedPoint diagonal = pMax - pMin;
    IntersectionCell *intersectionCell = new IntersectionCell;
    int8_t whichCoordinate = 2;

    if (diagonal.x >= diagonal.y && diagonal.x >= diagonal.z)
    {
        whichCoordinate = 0;
    }
    else if (diagonal.y >= diagonal.x && diagonal.y >= diagonal.z)
    {
        whichCoordinate = 1;
    }

    intersectionCell->pMin = pMin;
    intersectionCell->pMax = pMax;
    intersectionCell->pMax[whichCoordinate] -= (diagonal[whichCoordinate] / 2.0f);

    pMin[whichCoordinate] = intersectionCell->pMax[whichCoordinate];

    iNode = triangles.head();
    while (iNode != nullptr)
    {
        triangle = (AdvancedTriangle *)iNode->data;
        iNode = iNode->next();

        if (!isTriangleBoundingBoxInCell(triangle))
        {
            triangles.moveNodeTo((iNode != nullptr) ?
                    (iNode->prev()) : triangles.tail(), &(intersectionCell->triangles));
        }
        else if (intersectionCell->isTriangleBoundingBoxInCell(triangle))
        {
            intersectionCell->triangles.appendHead(triangle);
        }
    }

    return intersectionCell;
}

void IntersectionCell::selectIntersections(bool justProper)
{
    Node *iNode, *jNode;
    AdvancedTriangle *iTriangle, *jTriangle;
    List *triangleList;

    for (iNode = triangles.head(); iNode != nullptr; iNode = iNode->next())
    {
        for (jNode = iNode->next(); jNode != nullptr; jNode = jNode->next())
        {
            // For any pair of triangles in the cell, the same triangle pair can be in different
            // cells. The following avoids redoing the check.
            iTriangle = (AdvancedTriangle *) iNode->data;
            jTriangle = (AdvancedTriangle *) jNode->data;

            if (iTriangle->info == nullptr || jTriangle->info == nullptr ||
                    (((List *)iTriangle->info)->containsNode(jTriangle) == nullptr))
            {
                if (iTriangle->intersects(jTriangle, justProper))
                {
                    MARK_VISIT(iTriangle);
                    MARK_VISIT(jTriangle);

                    triangleList = ((iTriangle->info != nullptr) ?
                                    ((List *)iTriangle->info) : (new List));

                    triangleList->appendTail(jTriangle);
                    iTriangle->info = triangleList;

                    triangleList = ((jTriangle->info != nullptr) ?
                                    ((List *)jTriangle->info) : (new List));

                    triangleList->appendTail(iTriangle);
                    jTriangle->info = triangleList;
                }
            }
        }
    }
}

uint64_t AdvancedMesh::getNumberSelfIntersectingFaces()
{
    // Generic
    AdvancedTriangle* triangle;
    AdvancedVertex* vertex;
    Node* node;
    uint16_t trianglesPerCell = 50;

    List *selectedTriangles = new List;
    List *selectedVertices = new List;

    // At leaset there is a selection
    bool isSelection = 0;
    TIMER_SET;
    LOOP_COUNTER_SET;
    LOOP_STARTS("Selecting Intersecting Triangles");
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, _triangles.numberElements());

        if (IS_VISITED(triangle))
        {
            isSelection = 1;
            selectedTriangles->appendTail(triangle);

            vertex=triangle->v1();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }

            vertex=triangle->v2();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }

            vertex=triangle->v3();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // If any triangle is selected, delete it an its vertices
    if (!isSelection)
    {
        delete(selectedTriangles);
        delete(selectedVertices);

        selectedTriangles = &_triangles;
        selectedVertices = &_vertices;
    }

    IntersectionCell *c2, *cell = new IntersectionCell(this, !isSelection);
    List cells, toProess(cell);

    int i = 0;
    while ((cell = (IntersectionCell *)toProess.popHead()) != nullptr)
    {
        if (i > DI_MAX_NUMBER_OF_CELLS ||
            cell->triangles.numberElements() <= trianglesPerCell)
            cells.appendHead(cell);
        else
        {
            i++;
            c2 = cell->fork();
            toProess.appendTail(cell);
            toProess.appendTail(c2);
        }
    }

    // Deselect everything and select only intersecting triangles
    deselectTriangles();

    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Releasing Triangles");
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, _triangles.numberElements());
        triangle->info = nullptr;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Selecting Intersections");
    FOR_EACH_NODE(cells, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, cells.numberElements());

        (((IntersectionCell *) node->data)->selectIntersections(false));
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Dispose memory allocated for cells
    FOR_EACH_VT_TRIANGLE(selectedTriangles, triangle, node)
    {
        if (triangle->info!=nullptr)
            delete((List *)triangle->info);
        triangle->info = nullptr;
    }

    while (cells.numberElements())
    {
        delete((IntersectionCell *) cells.popHead());
    }

    // Count selected triangles for final report and delete stored normals
    int numberIntersectingTriangles = 0;
    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Counting Intersecting Triangles");
    FOR_EACH_VT_TRIANGLE(selectedTriangles, triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, selectedTriangles->numberElements());

        if (IS_VISITED(triangle))
            numberIntersectingTriangles++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return numberIntersectingTriangles;
}

int AdvancedMesh::selectIntersectingTriangles(uint16_t trianglesPerCell,
                                              bool justProper)
{
    // Generic
    AdvancedTriangle* triangle;
    AdvancedVertex* vertex;
    Node* node;

    List *selectedTriangles = new List;
    List *selectedVertices = new List;

    // At leaset there is a selection
    bool isSelection = 0;
    TIMER_SET;
    LOOP_COUNTER_SET;
    LOOP_STARTS("Selecting Intersecting Triangles");
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, _triangles.numberElements());

        if (IS_VISITED(triangle))
        {
            isSelection = 1;
            selectedTriangles->appendTail(triangle);

            vertex=triangle->v1();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }

            vertex=triangle->v2();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }

            vertex=triangle->v3();
            if (!IS_BIT(vertex, 5))
            {
                MARK_BIT(vertex, 5);
                selectedVertices->appendTail(vertex);
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // If any triangle is selected, delete it an its vertices
    if (!isSelection)
    {
        delete(selectedTriangles);
        delete(selectedVertices);

        selectedTriangles = &_triangles;
        selectedVertices = &_vertices;
    }

    IntersectionCell *c2, *cell = new IntersectionCell(this, !isSelection);
    List cells, toProess(cell);

    int i = 0;
    while ((cell = (IntersectionCell *)toProess.popHead()) != nullptr)
    {
        if (i > DI_MAX_NUMBER_OF_CELLS ||
            cell->triangles.numberElements() <= trianglesPerCell)
            cells.appendHead(cell);
        else
        {
            i++;
            c2 = cell->fork();
            toProess.appendTail(cell);
            toProess.appendTail(c2);
        }
    }

    // Deselect everything and select only intersecting triangles
    deselectTriangles();

    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Releasing Triangles");
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, _triangles.numberElements());
        triangle->info = nullptr;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Selecting Intersections");
    FOR_EACH_NODE(cells, node)
    {
        // LOOP_PROGRESS_FRACTION(++COUNTER, cells.numberElements());

        (((IntersectionCell *) node->data)->selectIntersections(justProper));
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Dispose memory allocated for cells
    FOR_EACH_VT_TRIANGLE(selectedTriangles, triangle, node)
    {
        if (triangle->info!=nullptr)
            delete((List *)triangle->info);
        triangle->info = nullptr;
    }

    while (cells.numberElements())
        delete((IntersectionCell *) cells.popHead());

    // Count selected triangles for final report and delete stored normals
    int numberIntersectingTriangles = 0;
    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Counting Intersecting Triangles");
    FOR_EACH_VT_TRIANGLE(selectedTriangles, triangle, node)
    {
        LOOP_PROGRESS_FRACTION(++COUNTER, selectedTriangles->numberElements());

        if (IS_VISITED(triangle))
            numberIntersectingTriangles++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    if (numberIntersectingTriangles)
        LOG_WARNING("The mesh has [%d] intersecting triangles. DIRTY MESH!",
                    numberIntersectingTriangles);
    else
        LOG_SUCCESS("No Self Intersections Detected");

    FOR_EACH_VV_VERTEX(selectedVertices, vertex, node)
    {
        UNMARK_BIT(vertex, 5);
    }

    if (isSelection)
    {
        delete(selectedTriangles);
        delete(selectedVertices);
    }

    return numberIntersectingTriangles;
}

bool AdvancedMesh::safeCoordBackApproximation()
{
    Node *node;
    AdvancedVertex *vertex;

    deselectTriangles();

    FOR_EACH_VERTEX(vertex, node)
    {
        jitterCoordinate(vertex->x, 0);
        jitterCoordinate(vertex->y, 0);
        jitterCoordinate(vertex->z, 0);
    }

    AdvancedEdge *e;
    AdvancedVertex *ov1, *ov2;

    int pnos = 0, nos;
    nos = 0;
    FOR_EACH_EDGE(e, node)
    {
        if (e->overlaps())
            nos++;
    }

    do
    {
        pnos = nos;
        FOR_EACH_EDGE(e, node)
        {
            if (e->overlaps())
            {
                ov1 = e->t1->oppositeVertex(e);
                ov2 = e->t2->oppositeVertex(e);
                vertex = (AdvancedPoint::squaredTriangleArea3D(e->v1, e->v2, ov1) <
                          AdvancedPoint::squaredTriangleArea3D(e->v1, e->v2, ov2)) ? (ov1) : (ov2);
                for (int a = -1; a <= 1; a++)
                {
                    for (int b = -1; b <= 1; b++)
                    {
                        for (int c = -1; c <= 1; c++)
                        {
                            jitterCoordinate(vertex->x, a);
                            jitterCoordinate(vertex->y, b);
                            jitterCoordinate(vertex->z, c);

                            if (e->overlaps())
                            {
                                jitterCoordinate(vertex->x, -a);
                                jitterCoordinate(vertex->y, -b);
                                jitterCoordinate(vertex->z, -c);
                            }
                            else
                            {
                                a = b = c = 2;
                            }
                        }
                    }
                }
            }
        }

        nos = 0;

        FOR_EACH_EDGE(e, node)
        {
            if (e->overlaps())
            {
                nos++;
            }
        }
    } while (nos < pnos);

    return (nos == 0);
}

bool AdvancedMesh::strongIntersectionRemoval(int maxIterations)
{
    int iterationCount = 0;

    LOG_STATUS("Removing Self Intersections: STRONG");
    while ((++iterationCount) <= maxIterations && selectIntersectingTriangles())
    {
        for (int n = 1; n < iterationCount; ++n)
        {
            growSelection();
        }

        removeSelectedTriangles();
        removeSmallestComponents();

        fillHoles(_edges.numberElements(), false);
        coordBackApproximation();
        remintsSelectTrianglesInCubes(this);
    }

    if (iterationCount > maxIterations)
        return false;

    return true;
}

} 
