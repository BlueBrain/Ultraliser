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

#include "AdvancedMesh.h"
#include <utilities/Timer.h>

namespace Ultraliser
{

int AdvancedMesh::starTriangulateHole(AdvancedEdge *edge)
{
    // Boundary edge, should be fixed
    if (!edge->isOnBoundary())
        return 0;

    List bvs;
    Node *iNode;
    AdvancedEdge *e1, *e2, *e3;
    AdvancedPoint np;
    AdvancedVertex *vertex;
    AdvancedVertex *nv, *v1, *v2;

    // Holes
    int numberTriangulatedHoles = 0;
    vertex = edge->v1;

    do {
        bvs.appendHead(vertex);
        vertex = vertex->nextOnBoundary();
        // MARK_BIT(v, 5);
    } while (vertex != edge->v1);

    FOR_EACH_VV_VERTEX((&(bvs)), vertex, iNode)
    {
        np = np + (*vertex);
    }

    np = np / bvs.numberElements();
    nv = newVertex(&np);
    _vertices.appendHead(nv);

    v1 = ((AdvancedVertex *)bvs.head()->data);
    AdvancedEdge *ep = v1->e0;
    e1 = createEdge(nv, v1);
    v1->e0 = ep;

    for (iNode = bvs.head()->next(); iNode != nullptr; iNode = iNode->next())
    {
        v2 = ((AdvancedVertex *)iNode->data);
        e2 = createEdge(nv, v2);
        e3 = v1->getEdge(v2);

        createTriangle(e1, e2, e3);
        numberTriangulatedHoles++;

        v1 = v2;
        e1 = e2;
    }

    v2 = ((AdvancedVertex *) bvs.head()->data);
    e2 = nv->getEdge(v2);
    e3 = v1->getEdge(v2);

    // Create the triangle
    createTriangle(e1, e2, e3);
    numberTriangulatedHoles++;

    // Create the number of triangles that were created to fill the holes
    return numberTriangulatedHoles;
}

int AdvancedMesh::TriangulateHole(AdvancedEdge *edge, AdvancedPoint *normal)
{
    if (!edge->isOnBoundary())
        return 0;

    List vertexList;
    Node *iNode, *currentNode = nullptr;
    AdvancedEdge *e1, *e2;
    AdvancedVertex *vertex, *v1, *v2;
    double angle, currentAngle;

    // Number of triangles created
    uint64_t numberTriangles = 0;

    vertex = edge->v1;
    do
    {
        vertexList.appendHead(vertex);
        vertex = vertex->nextOnBoundary();

    } while (vertex != edge->v1);

    while (vertexList.numberElements() > 2)
    {
        currentAngle = DBL_MAX;

        FOR_EACH_VV_VERTEX((&(vertexList)), vertex, iNode)
        {
            if (!IS_BIT(vertex, 5) && vertex->e0 &&
                    (angle = vertex->getAngleOnAveragePlane(normal)) < currentAngle)
            {
                currentAngle = angle;
                currentNode = iNode;
            }
        }
        if (currentAngle == DBL_MAX)
        {
            LOG_WARNING("TriangulateHole: CANNOT complete the triangulation!");
            FOR_EACH_VV_VERTEX((&(vertexList)), vertex, iNode)
            {
                UNMARK_BIT(vertex, 5);
            }
            return 0;
        }

        vertex = ((AdvancedVertex *)currentNode->data);

        v1 = (AdvancedVertex *)((currentNode->next() != nullptr) ?
                                    (currentNode->next()) : (vertexList.head()))->data;
        v2 = (AdvancedVertex *)((currentNode->prev() != nullptr) ?
                                    (currentNode->prev()) : (vertexList.tail()))->data;

        e1 = vertex->getEdge(v1);
        e2 = vertex->getEdge(v2);

        if (!eulerEdgeTriangle(e1, e2))
            MARK_BIT(vertex, 5);
        else
        {
            vertexList.removeCell(currentNode);

            UNMARK_BIT(v1, 5);
            UNMARK_BIT(v2, 5);

            numberTriangles++;
        }
    }

    int64_t i, skips;
    do
    {
        skips = 0;
        for (iNode=_edges.head(), i = 2 * numberTriangles * numberTriangles;
             i < UI2I64(numberTriangles); iNode = iNode->next(), i--)
        {
            edge = ((AdvancedEdge *)iNode->data);
            angle = edge->delaunayMinAngle();

            if (edge->swap())
            {
                if (edge->delaunayMinAngle() <= angle)
                    edge->swap(1);
                else
                    skips++;
            }
        }

        if (i < 0)
        {
            LOG_WARNING("Optimization is taking too long. Giving up!");
            break;
        }
    } while (skips);

    return numberTriangles;
}

int AdvancedMesh::TriangulateHole(AdvancedEdge *edge, List *vertexList)
{
    if (!edge->isOnBoundary())
        return 0;

    List boundaryVertexList, boundaryVertexLists;
    Node *node, *currentNode = nullptr;
    AdvancedEdge *e1, *e2;
    AdvancedVertex *vertex, *v1, *v2;
    double angle, currentAngle;

    uint64_t numberTriangles = 0, numberEdges;

    vertex = edge->v1;

    do {
        boundaryVertexList.appendHead(vertex);
        vertex = vertex->nextOnBoundary();
    } while (vertex != edge->v1);
    boundaryVertexLists.appendList(&boundaryVertexList);

    // While there are more than two boundary vertices
    List edgeList;
    while (boundaryVertexList.numberElements() > 2)
    {
        currentAngle = DBL_MAX;
        FOR_EACH_VV_VERTEX((&(boundaryVertexList)), vertex, node)
        {
            if (!IS_BIT(vertex, 5) && vertex->e0 &&
                    (angle = vertex->getAngleForTriangulation()) < currentAngle)
            {
                currentAngle = angle;
                currentNode = node;
            }
        }

        if (currentAngle == DBL_MAX)
        {
            LOG_WARNING("TriangulateHole: CANNOT complete the triangulation!");

            FOR_EACH_VV_VERTEX((&(boundaryVertexList)), vertex, node)
            {
                UNMARK_BIT(vertex, 5);
            }
            return 0;
        }

        vertex = ((AdvancedVertex *)currentNode->data);
        v1 = (AdvancedVertex *)((currentNode->next() != nullptr) ?
                                    (currentNode->next()) : (boundaryVertexList.head()))->data;
        v2 = (AdvancedVertex *)((currentNode->prev() != nullptr)
                                ? (currentNode->prev()) : (boundaryVertexList.tail()))->data;

        e1 = vertex->getEdge(v1);
        e2 = vertex->getEdge(v2);

        numberEdges = _edges.numberElements();
        if (!eulerEdgeTriangle(e1, e2)) MARK_BIT(vertex, 5);
        else
        {
            boundaryVertexList.removeCell(currentNode);

            UNMARK_BIT(v1, 5);
            UNMARK_BIT(v2, 5);

            numberTriangles++;

            if (_edges.numberElements() > numberEdges)
                edgeList.appendHead(_edges.head()->data);
        }
    }


    uint64_t i;
    AdvancedPoint triangleNormal;

    for (i = 0, node = _triangles.head(); i < numberTriangles; i++, node = node->next())
        triangleNormal = triangleNormal+((AdvancedTriangle *)node->data)->getNormal();

    if (triangleNormal.isNull())
    {
        LOG_WARNING("TriangulateHole: Unable to compute an average normal. CANNOT Optimize!");
        return numberTriangles;
    }
    triangleNormal.normalize();

    double *ovps = new double[3*(boundaryVertexLists.numberElements())];
    int j = 0;
    FOR_EACH_VV_VERTEX((&(boundaryVertexLists)), vertex, node)
    {
        ovps[j++] = vertex->x;
        ovps[j++] = vertex->y;
        ovps[j++] = vertex->z;

        vertex->project(&triangleNormal);
    }

    AdvancedPoint *p;
    double *ovpsi = new double[3*(vertexList->numberElements())];
    j = 0;
    FOR_EACH_NODE((*vertexList), node)
    {
        p = ((AdvancedPoint *)node->data);

        ovpsi[j++] = p->x;
        ovpsi[j++] = p->y;
        ovpsi[j++] = p->z;
    }
    FOR_EACH_NODE((*vertexList), node)
    {
        p = ((AdvancedPoint *)node->data);
        p->project(&triangleNormal);
    }

    int sw;
    do
    {
        sw = 0;
        FOR_EACH_VE_EDGE((&(edgeList)), e1, node)
        {
            angle = e1->delaunayMinAngle();
            if (e1->swap())
            {if (e1->delaunayMinAngle() <= angle) e1->swap(1); else sw++;}
        }
    } while (sw);

    int64_t totalNumberTriangles = _triangles.numberElements() - numberTriangles;
    List ivs;

    FOR_EACH_NODE((*vertexList), node)
    {
        p = ((AdvancedPoint *)node->data);
        ivs.appendTail(_watsonInsert(p, &_triangles,
                                     _triangles.numberElements() - totalNumberTriangles));
    }

    numberTriangles = (_triangles.numberElements() - totalNumberTriangles);

    j = 0;
    FOR_EACH_VV_VERTEX((&(ivs)), vertex, node)
    {
        if (vertex != nullptr)
        {
            vertex->x = ovpsi[j++];
            vertex->y = ovpsi[j++];
            vertex->z = ovpsi[j++];

        }
        else
            j+=3;
    }
    delete [] ovpsi;

    j = 0;
    FOR_EACH_VV_VERTEX((&(boundaryVertexLists)), vertex, node)
    {
        vertex->x = ovps[j++];
        vertex->y = ovps[j++];
        vertex->z = ovps[j++];
    }
    delete [] ovps;

    return numberTriangles;
}

AdvancedVertex *AdvancedMesh::_watsonInsert(AdvancedPoint *point,
                                            List *delaunayTriangulation,
                                            int numberTriangles)
{
    Node *iNode, *jNode;
    AdvancedEdge *edge;
    AdvancedTriangle *triangle;
    List vertexList, verticesList, toProcess, *incidentEdges;
    AdvancedVertex *v1, *v2, *v3;
    int i;

    for (i = 0, iNode = _triangles.head(); i < numberTriangles; iNode = iNode->next(), i++)
    {
        triangle = ((AdvancedTriangle *) iNode->data);
        if (triangle->edge1 != nullptr && triangle->inSphere(point))
        {
            v1 = triangle->v1();
            v2 = triangle->v2();
            v3 = triangle->v3();

            if (!IS_BIT(v1, 5))
                vertexList.appendHead(v1);

            if (!IS_BIT(v2, 5))
                vertexList.appendHead(v2);

            if (!IS_BIT(v3, 5))
                vertexList.appendHead(v3);

            MARK_BIT(v1, 5);
            MARK_BIT(v2, 5);
            MARK_BIT(v3, 5);
            MARK_BIT(triangle, 6);
            toProcess.appendHead(triangle);
        }
    }

    if (vertexList.numberElements() == 0)
        return nullptr;

    FOR_EACH_VV_VERTEX((&(vertexList)), v1, iNode)
    {
        incidentEdges = v1->getIncidentEdges();
        FOR_EACH_VE_EDGE(incidentEdges, edge, jNode)
        {
            if (!IS_BIT(edge->t1, 6) || !IS_BIT(edge->t2, 6))
                v1->e0 = edge;
        }
        delete(incidentEdges);
    }

    while (toProcess.numberElements())
    {
        triangle = ((AdvancedTriangle *)toProcess.head()->data);
        toProcess.removeCell(toProcess.head());
        unlinkTriangleNoManifold(triangle);
    }

    Node *tmpNode;
    for (i = 0, iNode = _triangles.head(); i<numberTriangles; i++)
    {
        triangle = ((AdvancedTriangle *)iNode->data);

        if (triangle->edge1 == nullptr)
        {
            tmpNode = iNode;
            iNode = iNode->next();
            _triangles.freeCell(tmpNode);
        }

        else iNode = iNode->next();
    }

    for (iNode = vertexList.head(); iNode != nullptr;)
    {
        v1 = ((AdvancedVertex *)iNode->data);
        if (v1->e0 == nullptr)
        {
            tmpNode = iNode;
            iNode = iNode->next();
            vertexList.removeCell(tmpNode);
        }
        else
            iNode = iNode->next();
    }

    v1 = v2 = ((AdvancedVertex *)vertexList.head()->data);
    do
    {
        verticesList.appendHead(v1);
        v1 = v1->nextOnBoundary();
    } while (v1 != v2);

    AdvancedVertex *v = newVertex(point->x, point->y, point->z);
    _vertices.appendHead(v);

    v1 = ((AdvancedVertex *)verticesList.head()->data);
    v->e0 = edge = newEdge(v, v1);
    UNMARK_BIT(v1, 5);
    _edges.appendHead(edge);

    for (iNode = verticesList.head()->next(); iNode != nullptr; iNode = iNode->next())
    {
        v1 = ((AdvancedVertex *)iNode->data);
        UNMARK_BIT(v1, 5);
        v2 = ((AdvancedVertex *)iNode->prev()->data);
        edge = newEdge(v, v1);
        createTriangle(edge, v1->getEdge(v2), (AdvancedEdge *)_edges.head()->data);
        _edges.appendHead(edge);
    }
    eulerEdgeTriangle(v->e0, (AdvancedEdge *)_edges.head()->data);

    return v;
}

int AdvancedMesh::retriangulateVT(AdvancedVertex *vertex)
{
    AdvancedPoint vertexNormal;
    AdvancedEdge *edge, *e0 = vertex->e0->t1->oppositeEdge(vertex);
    List *incidentTriangles = vertex->VT();
    List oppositeEdges;
    AdvancedTriangle *triangle;
    Node *iNode, *node;

    int i, numberTriangles;

    FOR_EACH_VT_TRIANGLE(incidentTriangles, triangle, iNode)
    {
        edge = triangle->oppositeEdge(vertex);
        oppositeEdges.appendTail(triangle->prevEdge(edge));
        oppositeEdges.appendTail(edge);
        oppositeEdges.appendTail(triangle->nextEdge(edge));
        vertexNormal = vertexNormal+triangle->getNormal();
        unlinkTriangle(triangle);
    }

    removeUnlinkedElements();

    vertexNormal.normalize();

    numberTriangles = TriangulateHole(e0, &vertexNormal);

    for (iNode=_triangles.head(), i = 0; i < numberTriangles; i++, iNode=iNode->next())
    {
        triangle = ((AdvancedTriangle *)iNode->data);
        if (triangle->overlaps() || triangle->isExactlyDegenerate())
            break;
    }

    if (i < numberTriangles)
    {
        LOG_WARNING("Re-triangulation failed. Restoring...");

        for (iNode=_triangles.head(), i = 0; i < numberTriangles; i++, iNode=iNode->next())
            unlinkTriangle(((AdvancedTriangle *)iNode->data));

        node = oppositeEdges.head();
        FOR_EACH_VT_TRIANGLE(incidentTriangles, triangle, iNode)
        {
            triangle->edge1 = ((AdvancedEdge *)node->data); node=node->next();
            triangle->edge2 = ((AdvancedEdge *)node->data); node=node->next();
            triangle->edge3 = ((AdvancedEdge *)node->data); node=node->next();

            triangle->edge1->v1 = vertex;
            triangle->edge1->v2 = (triangle->edge2->t1 == nullptr) ?
                        (triangle->edge2->v1) : (triangle->edge2->v2);

            triangle->edge3->v1 = vertex;
            triangle->edge3->v2 = (triangle->edge2->t1 == nullptr) ?
                        (triangle->edge2->v2) : (triangle->edge2->v1);

            ((triangle->edge2->t1 == nullptr) ?
                        (triangle->edge2->t1) : (triangle->edge2->t2)) = triangle;

            triangle->edge1->t1 = triangle;
            triangle->edge3->t2 = triangle;
        }

        vertex->e0 = ((AdvancedTriangle *)incidentTriangles->head()->data)->edge1;
    }

    delete(incidentTriangles);

    return 1;
}

int AdvancedMesh::TriangulateHole(AdvancedEdge *edge)
{
    // Not boundary edge
    if (!edge->isOnBoundary())
        return 0;

    List vertexList;
    Node *node, *currentNode = nullptr;
    AdvancedEdge *e1, *e2;
    AdvancedVertex *vertex, *v1, *v2;
    AdvancedTriangle *triangle;

    vertex = edge->v1;
    triangle = (edge->t1!=nullptr) ? (edge->t1) : (edge->t2);
    if (triangle->nextEdge(edge)->isOnBoundary() && triangle->prevEdge(edge)->isOnBoundary())
        return 0;

    do
    {
        vertexList.appendHead(vertex);
        vertex = vertex->nextOnBoundary();
    } while (vertex != edge->v1);


    int numberTriangles = 0;
    double angle, currentAngle;
    while (vertexList.numberElements() > 2)
    {
        currentAngle = DBL_MAX;
        FOR_EACH_VV_VERTEX((&vertexList), vertex, node)
        {
            if (!IS_BIT(vertex, 5) && vertex->e0 &&
                    (angle = vertex->getAngleForTriangulation()) < currentAngle)
            {
                currentAngle = angle;
                currentNode = node;
            }
        }

        if (currentAngle == DBL_MAX)
        {
            LOG_WARNING("TriangulateHole: CANNOT complete the triangulation!");

            FOR_EACH_VV_VERTEX((&vertexList), vertex, node)
            {
                UNMARK_BIT(vertex, 5);
            }

            int i = 0;
            FOR_EACH_TRIANGLE(triangle, node)
            {
                if (i++ == numberTriangles)
                    break;
                else
                    unlinkTriangle(triangle);
            }

            removeUnlinkedElements();
            return 0;
        }

        vertex = ((AdvancedVertex *)currentNode->data);
        v1 = (AdvancedVertex *)((currentNode->next() != nullptr) ?
                                    (currentNode->next()) : (vertexList.head()))->data;
        v2 = (AdvancedVertex *)((currentNode->prev() != nullptr)
                                ? (currentNode->prev()) : (vertexList.tail()))->data;

        e1 = vertex->getEdge(v1);
        e2 = vertex->getEdge(v2);

        if ((triangle=eulerEdgeTriangle(e1, e2)) == nullptr)
            MARK_BIT(vertex, 5);
        else
        {
            vertexList.removeCell(currentNode);

            UNMARK_BIT(v1, 5);
            UNMARK_BIT(v2, 5);

            MARK_VISIT(triangle);

            numberTriangles++;
        }
    }

    return numberTriangles;
}

void AdvancedMesh::fillHole(AdvancedEdge *edge, bool refine)
{
    Node *node;
    AdvancedTriangle *triangle;
    AdvancedVertex *vertex;

    // Deslect all triangles
    deselectTriangles();

    // Unmark all triangles
    FOR_EACH_VERTEX(vertex, node)
    {
        UNMARK_BIT(vertex, 5);
    }
    int numberTriangles = TriangulateHole(edge);
    if (!numberTriangles)
        return;

    int i = 0;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        if (i++ == numberTriangles)
            break;
        else
            MARK_VISIT(triangle);
    }

    if (refine)
        refineSelectedHolePatches((AdvancedTriangle *)_triangles.head()->data);
}

uint64_t AdvancedMesh::fillHoles(const uint64_t minNumberBoundaryEdges,
                                 const bool refinePatches)
{
    TIMER_SET;

    uint64_t numberEdges = minNumberBoundaryEdges;
    if (numberEdges == 0)
        numberEdges = _edges.numberElements();

    // Processing data structures
    AdvancedVertex *vertex;
    AdvancedTriangle* triangle;
    Node* node;

    // At least a single triangle must be selected to proceed with the cleaning
    int isSelection = 0;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        if (IS_VISITED(triangle))
        {
            isSelection = 1;
            break;
        }
    }

    if (isSelection)
    {
        // Some triangles are selected
        FOR_EACH_TRIANGLE(triangle, node)
        {
            if (!IS_VISITED(triangle))
            {
                MARK_BIT(triangle->v1(), 6);
                MARK_BIT(triangle->v2(), 6);
                MARK_BIT(triangle->v3(), 6);
            }
        }
    }
    else
    {
        // Unmark all the vertices in the mesh
        FOR_EACH_VERTEX(vertex, node)
        {
            UNMARK_BIT(vertex, 6);
        }
    }

    uint64_t grd = 0;
    uint64_t tbds = 0;
    // uint64_t pct = 100;

    // A list of all the boundary edges in the mesh.
    List boundaries;

    // Collecting the boundary vertices
    AdvancedVertex* bounadryVertex;
    LOOP_COUNTER_SET;
    TIMER_RESET;
    LOOP_STARTS("Collecting Boundary Vertices");
    FOR_EACH_VERTEX(vertex, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _vertices.numberElements());

        grd = 0;

        // If this is a boundary vertex
        if (!IS_BIT(vertex, 6) && vertex->isOnBoundary())
        {
            tbds++;
            bounadryVertex = vertex;
            do
            {
                if (IS_BIT(bounadryVertex, 6))
                    grd=numberEdges+1;

                MARK_BIT(bounadryVertex, 6);
                grd++;
                bounadryVertex = bounadryVertex->nextOnBoundary();
            } while (bounadryVertex != vertex);

            if (grd <= numberEdges)
                boundaries.appendHead(bounadryVertex->nextBoundaryEdge());
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    FOR_EACH_VERTEX(vertex, node)
    {
        UNMARK_BIT(vertex, 5);
        UNMARK_BIT(vertex, 6);
    }

    deselectTriangles();

    // Once the holes are filled, we must triangulate them to proceed with a triangular valid mesh
    LOOP_COUNTER_RESET;
    TIMER_RESET;
    LOOP_STARTS("Triangulating Holes");
    FOR_EACH_NODE(boundaries, node)
    {
        if (TriangulateHole((AdvancedEdge *)node->data) && refinePatches)
        {
            triangle = (AdvancedTriangle *)_triangles.head()->data;
            refineSelectedHolePatches(triangle);
        }

        LOOP_PROGRESS(++COUNTER, boundaries.numberElements())
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Return the number of boundary edges (actual elements)
    return boundaries.numberElements();
}

int AdvancedMesh::refineSelectedHolePatches(AdvancedTriangle *inputTriangulation)
{
    Node *iNode;
    AdvancedTriangle *triangle, *t1, *t2;
    AdvancedEdge *edge, *nextEdge;
    AdvancedVertex *vertex;

    List regionList, edgedToSwap;
    if (inputTriangulation != nullptr)
    {
        if (!IS_VISITED(inputTriangulation))
        {
            LOG_ERROR("refineSelectedHolePatches: unexpected unselected inputTriangulation!");
        }

        UNMARK_VISIT(inputTriangulation);
        edgedToSwap.appendHead(inputTriangulation);

        while ((triangle = (AdvancedTriangle *) edgedToSwap.popHead()) != nullptr)
        {
            regionList.appendHead(triangle);

            t1 = triangle->t1();
            if (IS_VISITED(t1))
            {
                UNMARK_VISIT(t1);
                edgedToSwap.appendHead(t1);
            }

            t1 = triangle->t2();
            if (IS_VISITED(t1))
            {
                UNMARK_VISIT(t1);
                edgedToSwap.appendHead(t1);
            }

            t1 = triangle->t3();
            if (IS_VISITED(t1))
            {
                UNMARK_VISIT(t1);
                edgedToSwap.appendHead(t1);
            }
        }

        FOR_EACH_VT_TRIANGLE((&regionList), triangle, iNode)
        {
            MARK_VISIT(triangle);
        }
    }
    else
    {
        FOR_EACH_TRIANGLE(triangle, iNode)
        {
            if (IS_VISITED(triangle))
            {
                regionList.appendHead(triangle);
            }
        }
    }

    List allEdges;
    FOR_EACH_VT_TRIANGLE((&regionList), triangle, iNode)
    {
        edge = triangle->edge1;
        if (!IS_BIT(edge, 5))
        {
            MARK_BIT(edge, 5);
            allEdges.appendHead(edge);
        }
        else
            UNMARK_BIT(edge, 5);

        edge = triangle->edge2;
        if (!IS_BIT(edge, 5))
        {
            MARK_BIT(edge, 5);
            allEdges.appendHead(edge);
        }
        else
            UNMARK_BIT(edge, 5);

        edge = triangle->edge3;
        if (!IS_BIT(edge, 5))
        {
            MARK_BIT(edge, 5);
            allEdges.appendHead(edge);
        }
        else
            UNMARK_BIT(edge, 5);
    }

    List interiorEdges, boundaryEdges;
    while (allEdges.numberElements())
    {
        edge = (AdvancedEdge *)allEdges.popHead();
        if (IS_BIT(edge, 5))
        {
            boundaryEdges.appendHead(edge);
            UNMARK_BIT(edge, 5);
        }
        else
        {
            interiorEdges.appendHead(edge);
            MARK_BIT(edge, 5);
        }
    }

    List boundaryVertices;
    FOR_EACH_VE_EDGE((&boundaryEdges), edge, iNode)
    {
        vertex = edge->v1;
        if (!IS_BIT(vertex, 5))
        {
            MARK_BIT(vertex, 5);
            boundaryVertices.appendHead(vertex);
        }

        vertex = edge->v2;
        if (!IS_BIT(vertex, 5))
        {
            MARK_BIT(vertex, 5);
            boundaryVertices.appendHead(vertex);
        }
    }

    FOR_EACH_VV_VERTEX((&boundaryVertices), vertex, iNode)
    {
        UNMARK_BIT(vertex, 5);
    }

    Node *jNode;
    List *incidentEdges;

    double sigma;
    int  numberEdges;
    FOR_EACH_VV_VERTEX((&boundaryVertices), vertex, iNode)
    {
        incidentEdges = vertex->getIncidentEdges();
        sigma = 0;
        numberEdges = 0;

        FOR_EACH_VE_EDGE(incidentEdges, edge, jNode)
        {
            if (!IS_BIT(edge, 5))
            {
                numberEdges++;
                sigma += edge->length();
            }
        }

        sigma /= numberEdges;
        vertex->info = new double(sigma);
        delete(incidentEdges);
    }

    FOR_EACH_VE_EDGE((&interiorEdges), edge, iNode)
    {
        UNMARK_BIT(edge, 5);
    }

    FOR_EACH_VE_EDGE((&boundaryEdges), edge, iNode)
    {
        MARK_BIT(edge, 6);
    }

    const double alpha = sqrt(2.0);

    int swaps, totalIterations;
    int totalNumberTriangles = -1, currentNumberTriangles;
    int gits = 0;
    AdvancedPoint triangleCenter;
    List interiorVertices;
    do
    {
        currentNumberTriangles = totalNumberTriangles;
        totalNumberTriangles = 0;
        FOR_EACH_VT_TRIANGLE((&regionList), triangle, iNode)
        {
            triangleCenter = triangle->getCenter();

            const double v1Sigma = (*(double *)triangle->v1()->info);
            const double v2Sigma = (*(double *)triangle->v2()->info);
            const double v3Sigma = (*(double *)triangle->v3()->info);

            sigma = (v1Sigma + v2Sigma + v3Sigma) / 3.0;
            const double dv1 = alpha * (triangle->v1()->distance(&triangleCenter));
            const double dv2 = alpha * (triangle->v2()->distance(&triangleCenter));
            const double dv3 = alpha * (triangle->v3()->distance(&triangleCenter));

            if (dv1 > sigma && dv1 > v1Sigma &&
                dv2 > sigma && dv2 > v2Sigma &&
                dv3 > sigma && dv3 > v3Sigma)
            {
                const uint64_t numberTriangles = _triangles.numberElements();
                vertex = splitTriangle(triangle,&triangleCenter,1);
                totalNumberTriangles += (_triangles.numberElements() - numberTriangles);

                if (_triangles.numberElements() == numberTriangles + 2)
                {
                    vertex->info = new double(sigma);

                    interiorVertices.appendHead(vertex);
                    interiorEdges.appendHead(vertex->e0);
                    interiorEdges.appendHead(vertex->e0->leftTriangle(vertex)->prevEdge(vertex->e0));
                    interiorEdges.appendHead(vertex->e0->rightTriangle(vertex)->nextEdge(vertex->e0));

                    t1 = ((AdvancedTriangle *)_triangles.head()->data);
                    t2 = ((AdvancedTriangle *)_triangles.head()->next()->data);
                    t1->mask = t2->mask = triangle->mask;

                    regionList.appendHead(t1);
                    regionList.appendHead(t2);
                }
            }
        }

        FOR_EACH_VE_EDGE((&interiorEdges), edge, iNode)
        {
            MARK_BIT(edge, 5);
            edgedToSwap.appendHead(edge);
        }

        totalIterations = 0;
        swaps = 1;

        while (swaps && totalIterations++ < 10)
        {
            swaps = 0;
            while ((edge=(AdvancedEdge *)edgedToSwap.popHead()) != nullptr)
            {
                UNMARK_BIT(edge, 5);

                const double edgeLength = edge->squaredLength();
                if (edge->swap())
                {
                    if (edge->squaredLength() >= edgeLength * 0.999999)
                        edge->swap(1);
                    else
                    {
                        swaps++;
                        edgedToSwap.appendTail(edge);
                        nextEdge = edge->t1->nextEdge(edge);

                        if (!IS_BIT(nextEdge, 5) && !IS_BIT(nextEdge, 6))
                        {
                            MARK_BIT(nextEdge, 5);
                            edgedToSwap.appendTail(nextEdge);
                        }

                        nextEdge = edge->t1->prevEdge(edge);
                        if (!IS_BIT(nextEdge, 5) && !IS_BIT(nextEdge, 6))
                        {
                            MARK_BIT(nextEdge, 5);
                            edgedToSwap.appendTail(nextEdge);
                        }

                        nextEdge = edge->t2->nextEdge(edge);
                        if (!IS_BIT(nextEdge, 5) && !IS_BIT(nextEdge, 6))
                        {
                            MARK_BIT(nextEdge, 5);
                            edgedToSwap.appendTail(nextEdge);
                        }

                        nextEdge = edge->t2->prevEdge(edge);
                        if (!IS_BIT(nextEdge, 5) && !IS_BIT(nextEdge, 6))
                        {
                            MARK_BIT(nextEdge, 5);
                            edgedToSwap.appendTail(nextEdge);
                        }
                    }
                }
            }
        }

        if (currentNumberTriangles == totalNumberTriangles)
            gits++;
    } while (totalNumberTriangles && gits < 10);

    FOR_EACH_VV_VERTEX((&boundaryVertices), vertex, iNode)
    {
        delete((double *)vertex->info);
        vertex->info = nullptr;
        MARK_BIT(vertex, 5);
    }

    FOR_EACH_VV_VERTEX((&interiorVertices), vertex, iNode)
    {
        delete((double *)vertex->info);
        vertex->info = nullptr;
        MARK_BIT(vertex, 6);
    }

    if (gits >= 10)
    {
        LOG_WARNING("refineSelectedHolePatches: Refinement stage failed to converge. Breaking!");
        return 1;
    }

    return 0;
}

AdvancedEdge *AdvancedMesh::joinBoundaryLoops(AdvancedVertex *boundaryVertex1,
                                              AdvancedVertex *boundaryVertex2,
                                              const bool &justConnect,
                                              const bool &refine)
{
    AdvancedVertex *vertex, *nextVertex, *previousVertex;
    AdvancedEdge *edge, *boundaryVertexEdge1, *boundaryVertexEdge2;
    AdvancedTriangle *triangle;
    Node *iNode;


    if (boundaryVertex1 == nullptr || boundaryVertex2 == nullptr ||
       !boundaryVertex1->isOnBoundary() || !boundaryVertex2->isOnBoundary())
    {
        return nullptr;
    }

    FOR_EACH_VERTEX(vertex, iNode)
    {
        UNMARK_VISIT(vertex);
    }

    deselectTriangles();

    vertex = boundaryVertex1;
    if (!justConnect)
    {
        do {
            vertex = vertex->nextOnBoundary();
            if (vertex == boundaryVertex2)
                return nullptr;
        } while (vertex != boundaryVertex1);
    }
    else
    {
        nextVertex = boundaryVertex1->nextOnBoundary();
        previousVertex = boundaryVertex1->prevOnBoundary();

        if (boundaryVertex2 == nextVertex || boundaryVertex2 == previousVertex)
            return nullptr;

        if (boundaryVertex2 == nextVertex->nextOnBoundary())
        {
            triangle = eulerEdgeTriangle(nextVertex->prevBoundaryEdge(),
                                         nextVertex->nextBoundaryEdge());
            MARK_VISIT(triangle);

            return triangle->oppositeEdge(nextVertex);
        }

        if (boundaryVertex2 == previousVertex->prevOnBoundary())
        {
            triangle = eulerEdgeTriangle(previousVertex->prevBoundaryEdge(),
                                         previousVertex->nextBoundaryEdge());
            MARK_VISIT(triangle);

            return triangle->oppositeEdge(previousVertex);
        }
    }

    boundaryVertexEdge1 = boundaryVertex1->prevBoundaryEdge();
    nextVertex = boundaryVertexEdge1->oppositeVertex(boundaryVertex1);

    boundaryVertexEdge2 = boundaryVertex2->nextBoundaryEdge();
    previousVertex = boundaryVertexEdge2->oppositeVertex(boundaryVertex2);

    AdvancedEdge *justEdge = createEdge(boundaryVertex1, boundaryVertex2);
    AdvancedEdge *je1 = createEdge(boundaryVertex1, previousVertex);
    AdvancedEdge *je2 = createEdge(previousVertex, nextVertex);

    triangle = createTriangle(justEdge, boundaryVertexEdge2, je1);
    MARK_VISIT(triangle);

    triangle = createTriangle(je1, je2, boundaryVertexEdge1);
    MARK_VISIT(triangle);

    if (justConnect)
        return justEdge;

    // Compute the total edge length
    double totalEdgeLength1 = 0.0;
    vertex = boundaryVertex1;
    do {
        edge = vertex->nextBoundaryEdge();
        vertex = edge->oppositeVertex(vertex);
        totalEdgeLength1 += edge->length();
    } while (vertex != boundaryVertex1);

    double totalEdgeLength2 = 0.0;
    vertex = boundaryVertex2;
    do {
        edge = vertex->nextBoundaryEdge();
        vertex = edge->oppositeVertex(vertex);
        totalEdgeLength2 += edge->length();
    } while (vertex != boundaryVertex2);

    double pl1 = totalEdgeLength1;
    double pl2 = totalEdgeLength2;

    edge = justEdge;
    while (edge->isOnBoundary())
    {
        boundaryVertex1 = (edge->t2 != nullptr) ? (edge->v2) : (edge->v1);
        boundaryVertexEdge1 = boundaryVertex1->nextBoundaryEdge();

        boundaryVertex2 = (edge->t1 != nullptr) ? (edge->v2) : (edge->v1);
        boundaryVertexEdge2 = boundaryVertex2->prevBoundaryEdge();

        double const c1 = std::fabs((pl1 - boundaryVertexEdge1->length()) *
                                    totalEdgeLength2 - pl2 * totalEdgeLength1);
        double const c2 = std::fabs((pl2 - boundaryVertexEdge2->length()) *
                                    totalEdgeLength1 - pl1 * totalEdgeLength2);

        if (c1 < c2)
        {
            triangle = eulerEdgeTriangle(edge, boundaryVertexEdge1);
            MARK_VISIT(triangle);

            pl1 -= boundaryVertexEdge1->length();
            edge = triangle->nextEdge(boundaryVertexEdge1);
        }
        else
        {
            triangle = eulerEdgeTriangle(boundaryVertexEdge2, edge);
            MARK_VISIT(triangle);

            pl2 -= boundaryVertexEdge2->length();
            edge = triangle->prevEdge(boundaryVertexEdge2);
        }
    }

    if (refine)
        refineSelectedHolePatches();

    return justEdge;
}

}
