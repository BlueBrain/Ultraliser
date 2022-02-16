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

#include <common/Common.h>
#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/primitives/AdvancedEdge.h>
#include <data/meshes/advanced/primitives/AdvancedTriangle.h>
#include <data/meshes/advanced/AdvancedMesh.h>
#include <algorithms/Sorting.h>
#include <utilities/Timer.h>
#include <data/meshes/advanced/Defines.hh>

namespace Ultraliser
{

const char* AdvancedMesh::checkConnectivity()
{
    AdvancedVertex *vertex;
    AdvancedEdge *edge1,*edge2;
    AdvancedTriangle *triangle;
    Node *nNode,*mNode;
    List *vertexEdges;

    // Vertex connectivity
    FOR_EACH_VERTEX(vertex, nNode)
    {
        if (vertex == nullptr)
        {
            return "checkConnectivity: Detected nullptr element in the vertex list!";
        }

        if (vertex->e0 == nullptr)
        {
            return "checkConnectivity: Detected nullptr e0 pointer for a vertex!";
        }

        if (!vertex->e0->hasVertex(vertex))
        {
            return "checkConnectivity: Detected wrong e0 pointer for a vertex!";
        }
    }

    // Edge connectivity
    FOR_EACH_EDGE(edge1, nNode)
    {
        if (edge1 == nullptr)
        {
            return "checkConnectivity: Detected nullptr element in E list!";
        }

        if (edge1->v1 == nullptr || edge1->v2 == nullptr)
        {
            return "checkConnectivity: Detected edge with one or two nullptr end-points!";
        }

        if (edge1->v1 == edge1->v2)
        {
            return "checkConnectivity: Detected edge with two coincident end-points!";
        }

        if (edge1->t1 == nullptr && edge1->t2 == nullptr)
        {
            return "checkConnectivity: Detected edge with no incident triangles!";
        }

        if (edge1->t1 != nullptr)
        {
            if (!edge1->t1->hasEdge(edge1))
            {
                return "checkConnectivity: Detected wrong t1 triangle at an edge";
            }
            if (edge1->commonVertex(edge1->t1->nextEdge(edge1)) == edge1->v1)
            {
                return "checkConnectivity: Edge orientation does not match t1 normal";
            }
        }

        if (edge1->t2 != nullptr)
        {
            if (!edge1->t2->hasEdge(edge1))
            {
                return "checkConnectivity: Detected wrong t2 triangle at an edge";
            }
            if (edge1->commonVertex(edge1->t2->nextEdge(edge1)) == edge1->v2)
            {
                return "checkConnectivity: Edge orientation does not match t2 normal";
            }
        }
    }

    // Triangle connectivity
    FOR_EACH_TRIANGLE(triangle, nNode)
    {
        if (triangle == nullptr)
        {
            return "checkConnectivity: Detected nullptr element in T list!";
        }

        if (triangle->edge1 == nullptr ||
            triangle->edge2 == nullptr ||
            triangle->edge3 == nullptr)
        {
            return "checkConnectivity: Detected nullptr as a triangle edge!";
        }

        if (triangle->edge1 == triangle->edge2 ||
            triangle->edge1 == triangle->edge3 ||
            triangle->edge2 == triangle->edge3)
        {
            return "checkConnectivity: Detected triangle with two coincident edges!";
        }

        if (triangle->v1() == nullptr ||
            triangle->v2() == nullptr ||
            triangle->v3() == nullptr)
        {
            return "checkConnectivity: Triangle edges do not share vertices!";
        }

        if (triangle->edge1->t1 != triangle && triangle->edge1->t2 != triangle)
        {
            return "checkConnectivity: Detected triangle with 1st edge not pointing to the "
                   "triangle itself!";
        }

        if (triangle->edge2->t1 != triangle && triangle->edge2->t2 != triangle)
        {
            return "checkConnectivity: Detected triangle with 2nd edge not pointing to the "
                   "triangle itself!";
        }

        if (triangle->edge3->t1 != triangle && triangle->edge3->t2 != triangle)
        {
            return "checkConnectivity: Detected triangle with 3rd edge not pointing to the "
                   "triangle itself!";
        }
    }

    // Now go edge by edge
    FOR_EACH_EDGE(edge1, nNode)
    {
        vertexEdges = edge1->v1->getIncidentEdges();
        FOR_EACH_VE_EDGE(vertexEdges, edge2, mNode)
        {
            if (edge2 != edge1 && edge2->oppositeVertex(edge1->v1) == edge1->v2)
            {
                return "checkConnectivity: Detected duplicate edge!";
            }
        }

        if (vertexEdges->containsNode(edge1) == nullptr)
        {
            return "checkConnectivity: Detected non manifold vertex!";
        }

        delete(vertexEdges);

        vertexEdges = edge1->v2->getIncidentEdges();
        FOR_EACH_VE_EDGE(vertexEdges, edge2, mNode)
        {
            if (edge2 != edge1 && edge2->oppositeVertex(edge1->v2) == edge1->v1)
            {
                return "checkConnectivity: Detected duplicate edge!";
            }
        }

        if (vertexEdges->containsNode(edge1) == nullptr)
        {
            return "checkConnectivity: Detected non manifold vertex!";
        }

        delete(vertexEdges);
    }

    return nullptr;
}

int AdvancedMesh::duplicateNonManifoldVertices()
{
    // Generic
    AdvancedVertex *vertex;
    AdvancedEdge *edge1, *edge2;
    Node *node1, *node2;
    List *vertexEdges;

    // The number of duplicated vertices that indicate the number of singular vertices in the mesh
    int numberDuplicatedVertices = 0;

    FOR_EACH_EDGE(edge1, node1)
    {
        vertexEdges = edge1->v1->getIncidentEdges();

        if (vertexEdges->containsNode(edge1) == nullptr)
        {
            vertex = newVertex(edge1->v1);
            vertex->info = edge1->v1->info;
            vertex->mask = 0;
            _vertices.appendHead(vertex);

            FOR_EACH_VE_EDGE(vertexEdges, edge2, node2)
            {
                edge2->replaceVertex(edge1->v1, vertex);
            }

            vertex->e0 = edge1->v1->e0;
            edge1->v1->e0 = edge1;

            numberDuplicatedVertices++;
        }

        delete(vertexEdges);
    }

    FOR_EACH_EDGE(edge1, node1)
    {
        vertexEdges = edge1->v2->getIncidentEdges();

        if (vertexEdges->containsNode(edge1) == nullptr)
        {
            vertex = newVertex(edge1->v2);
            vertex->info = edge1->v2->info;
            vertex->mask = 0;
            _vertices.appendHead(vertex);

            FOR_EACH_VE_EDGE(vertexEdges, edge2, node2)
            {
                edge2->replaceVertex(edge1->v2, vertex);
            }

            vertex->e0 = edge1->v2->e0;
            edge1->v2->e0 = edge1;

            numberDuplicatedVertices++;
        }

        delete(vertexEdges);
    }

    // Update the state
    if (numberDuplicatedVertices)
        _dBoundaries = _dHandles = _dShells = 1;

    // Return the number of duplicated vertices
    return numberDuplicatedVertices;
}

AdvancedVertex *AdvancedMesh::checkGeometry()
{
    // The closest vertex that will be returned when there is an issue during the geometry check
    AdvancedVertex *closestVertexToIssue = nullptr;

    // Generic
    AdvancedTriangle *triangle;
    AdvancedEdge *edge;

    // If this return a nullptr, there is no enough memory to check
    AdvancedVertex **vertexArray = (AdvancedVertex **) _vertices.toArray();
    AdvancedEdge **edgeArray;
    AdvancedVertex *vertex1, *vertex2;
    Node *node;

    if (vertexArray == nullptr)
    {
        LOG_WARNING("checkGeometry: Not enough memory. Cannot check for coincident vertices.");
    }
    else
    {
        jqSort((void **) vertexArray, _vertices.numberElements(), xyzCompare);

        for (uint64_t i = 0; i < (_vertices.numberElements() - 1); ++i)
        {
            vertex1 = ((AdvancedVertex *) vertexArray[i]);
            vertex2 = ((AdvancedVertex *) vertexArray[i + 1]);

            if ((*vertex1) == (*vertex2))
            {
                closestVertexToIssue = vertex1;

                LOG_WARNING("checkGeometry: detected coincident vertices.");

                if (vertex1->getEdge(vertex2))
                {
                    LOG_WARNING("and there is an edge connecting them!");

                    free(vertexArray);

                    return vertex1;
                }
            }
        }
        free(vertexArray);
    }

    edgeArray = (AdvancedEdge **)_edges.toArray();

    if (edgeArray == nullptr)
    {
        LOG_WARNING("Not enough memory. Could NOT check for coincident edges.");
    }
    else
    {
        jqSort((void **)edgeArray, _edges.numberElements(), lexEdgeCompare);

        for (uint64_t i = 0; i < (_edges.numberElements() - 1); ++i)
        {
            if (!lexEdgeCompare(edgeArray[i], edgeArray[i + 1]))
            {
                closestVertexToIssue = ((AdvancedEdge *)edgeArray[i])->v1;

                LOG_WARNING("checkGeometry: detected coincident edges.");
            }
        }

        free(edgeArray);
    }

    double angle, minda = 0;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        angle = triangle->getAngle(triangle->v1());
        if (angle == 0 || angle == M_PI)
        {
            LOG_WARNING("Degenerate triangle detected.");
            return triangle->v1();
        }

        angle = triangle->getAngle(triangle->v2());
        if (angle == 0 || angle == M_PI)
        {
            LOG_WARNING("Degenerate triangle detected.");
            return triangle->v2();
        }

        angle = triangle->getAngle(triangle->v3());
        if (angle == 0 || angle == M_PI)
        {
            LOG_WARNING("Degenerate triangle detected.");
            return triangle->v3();
        }
    }

    angle = minda = 0;
    FOR_EACH_EDGE(edge, node)
    {
        if (edge->t1 != nullptr && edge->t2 != nullptr &&
           (angle = edge->t1->getDAngle(edge->t2)) == M_PI)
        {
            LOG_WARNING("checkGeometry: overlapping triangles detected.");
            return edge->v1;
        }
        else
        {
            minda = MAX(minda, angle);
        }
    }

    LOG_INFO("checkGeometry: minimum dihedral angle = %f (%f DEGs)",
                    M_PI - minda, ((M_PI - minda) * 360) / (2 * M_PI));
    return closestVertexToIssue;
}

int AdvancedMesh::mergeCoincidentEdges()
{
    _vertices.sort(&xyzCompare);
    Node *iNode;
    AdvancedVertex *vertex, *pVertex = (AdvancedVertex *)_vertices.head()->data;

    FOR_EACH_VERTEX(vertex, iNode)
    {
        UNMARK_BIT(vertex, 5);
    }

    AdvancedEdge *edge;
    FOR_EACH_EDGE(edge, iNode)
    {
        if (edge->isOnBoundary())
        {
            MARK_BIT(edge->v1, 5);
            MARK_BIT(edge->v2, 5);
        }
    }

    FOR_EACH_VERTEX(vertex, iNode)
    {
        if ((*vertex) != (*pVertex) || !IS_BIT(vertex,5))
            pVertex = vertex;

        vertex->info = pVertex;
        UNMARK_BIT(vertex, 5);
    }

    // At this point any vertex points (through 'info') to its unique representative (possibly itself)
    FOR_EACH_VERTEX(vertex, iNode)
    {
        vertex->e0 = nullptr;
    }

    FOR_EACH_EDGE(edge, iNode)
    {
        if (edge->v1->info != edge->v1)
            edge->v1 = (AdvancedVertex *)edge->v1->info;

        if (edge->v2->info != edge->v2)
            edge->v2 = (AdvancedVertex *)edge->v2->info;

        edge->v1->e0 = edge->v2->e0 = edge;
    }

    removeVertices();

    // At this point the mesh should no longer have duplicated vertices, but may have duplicated edges
    _edges.sort(&vtxEdgeCompare);
    AdvancedEdge *pEdge = (AdvancedEdge *)_edges.head()->data;
    FOR_EACH_EDGE(edge, iNode)
    {
        if (!edge->isOnBoundary() || vtxEdgeCompare(edge, pEdge))
            pEdge = edge;
        edge->info = pEdge;
    }

    FOR_EACH_EDGE(edge, iNode)
    {
        if (edge->info != edge)
        {
            AdvancedTriangle *boundaryTriangle = edge->getBoundaryTriangle();
            AdvancedEdge *keyEdge = ((AdvancedEdge *)edge->info);
            // AdvancedTriangle *otherTriangle = keyEdge->getBoundaryTriangle();

            boundaryTriangle->replaceEdge(edge, keyEdge);
            ((keyEdge->t1 == nullptr) ? (keyEdge->t1) : (keyEdge->t2)) = boundaryTriangle;

            edge->v1 = edge->v2 = nullptr;
            keyEdge->v1->e0 = keyEdge->v2->e0 = keyEdge;
        }
    }
    removeUnlinkedElements();

    return 1;
}

bool AdvancedMesh::fixConnectivity()
{
    // True if connectivity is fixed and false otherwise.
    bool returnValue = true;

    // Counting the number of issues
    int counter = 0;

    // Remove flying or unwanted vertices that are flying aroundand does NOT
    // contrinute to the geometry of the mesh
    if ((counter = removeVertices()))
    {
        _loadedNumberFloatingVertices = counter;

        returnValue = false;
        LOG_WARNING("Extra Vertices: [%d] isolated vertices have been removed.",
                    counter);
    }
    else
    {
        _loadedNumberFloatingVertices = 0;
    }

    // Try to turn to manifold
    counter = cutAndStitch();
    if (counter)
    {
        // if (!_initialValuesUpdated)
            _loadedNumberSingularEdges = counter;

        returnValue = false;
        LOG_WARNING("Manifoldness: Some cuts were necessary to cope with "
                    "input non-manifold configuration.");
    }
    else
    {
        _loadedNumberSingularEdges = 0;
    }

    // Duplicate the non-manifold vertices, if any
    counter = duplicateNonManifoldVertices();
    if (counter > 0)
    {
        LOG_WARNING("%d", _loadedNonManifoldVertices);
        _loadedNonManifoldVertices = counter;

        returnValue=false;
        LOG_WARNING("Manifoldness: [%d] non-manifold Vertices have been "
                    "Duplicated.",counter);
    }
    else
    {
        _loadedNonManifoldVertices = 0;
    }

    // Fix the orientation of the faces if were inverted
    counter = forceNormalConsistence();
    if (counter)
    {
        if (!_initialValuesUpdated)
            _loadedOrientationFlag = counter;

        returnValue = false;
        LOG_WARNING("Orientation: Some triangles have been reversed to achieve "
                    "correct orientation.");
    }
    else
    {
         _loadedOrientationFlag = 0;
    }

    counter = removeDuplicatedTriangles();
    if (counter > 0)
    {
        if (!_initialValuesUpdated)
            _loadedDuplicateTriangles = counter;
        returnValue=false;
        LOG_WARNING("Triangles: [%d] double-triangles have been removed.", counter);
    }

    _initialValuesUpdated = true;
    return returnValue;
}

bool AdvancedMesh::rebuildConnectivity(bool fixMeshConnectivity)
{
    // If no vertices, return
    if (_vertices.numberElements() == 0)
        return false;

    // Sort the vertices
    _vertices.sort(&xyzCompare);

    Node *iNode;
    AdvancedVertex *vertex, *pv=(AdvancedVertex *)_vertices.head()->data;

    FOR_EACH_VERTEX(vertex, iNode)
    {
        if ((*vertex)!=(*pv))
            pv = vertex;
        vertex->info = pv;
    }

    // At this point any vertex points (through 'info') to its unique representative (possibly itself)
    FOR_EACH_VERTEX(vertex, iNode)
    {
        vertex->e0 = nullptr;
    }

    AdvancedEdge *iEdge;
    FOR_EACH_EDGE(iEdge, iNode)
    {
        if (iEdge->v1->info != iEdge->v1)
            iEdge->v1 = (AdvancedVertex *)iEdge->v1->info;

        if (iEdge->v2->info != iEdge->v2)
            iEdge->v2 = (AdvancedVertex *)iEdge->v2->info;

        iEdge->v1->e0 = iEdge->v2->e0 = iEdge;
    }

    removeVertices();

    // At this point the mesh should no longer have duplicated vertices, but may have duplicated edges
    AdvancedTriangle *triangle;
    ExtendedVertex **extendedVertex = new ExtendedVertex *[_vertices.numberElements()];
    uint64_t i = 0;
    FOR_EACH_VERTEX(vertex, iNode)
    {
        vertex->e0 = nullptr;
        extendedVertex[i] = new ExtendedVertex(vertex);
        vertex->info = (void *)i;
        i++;
    }

    uint64_t numberTriangles = _triangles.numberElements();
    int *triangles = new int[numberTriangles * 3];
    i = 0;
    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        const int* x1 = (int*) triangle->v1()->info;
        const int* x2 = (int*) triangle->v2()->info;
        const int* x3 = (int*) triangle->v3()->info;

        triangles[i * 3 + 0] = *x1;
        triangles[i * 3 + 1] = *x2;
        triangles[i * 3 + 2] = *x3;
        i++;
    }
    _triangles.freeNodes();
    _edges.freeNodes();

    for (i = 0; i < numberTriangles; ++i)
    {
        const int v1 = triangles[i * 3 + 0];
        const int v2 = triangles[i * 3 + 1];
        const int v3 = triangles[i * 3 + 2];

        if (v1 != v2 && v2 != v3 && v1 != v3)
            createIndexedTriangle(extendedVertex, v1, v2, v3);
    }

    for (uint64_t j=0; j < _vertices.numberElements(); ++j)
    {
        delete(extendedVertex[j]);
    }

    delete [] extendedVertex;
    delete [] triangles;

    if(fixMeshConnectivity)
        return fixConnectivity();
    else
        return true;
}

int AdvancedMesh::removeDuplicatedTriangles()
{
    AdvancedEdge *edge;
    Node *node;
    AdvancedPoint point;

    uint64_t i = 0;
    FOR_EACH_EDGE(edge, node)
    {
        if (!edge->isOnBoundary() &&
             edge->t1->oppositeVertex(edge) == edge->t2->oppositeVertex(edge))
        {
            unlinkTriangle(edge->t2);
            i++;
        }
    }
    removeUnlinkedElements();

    if (i)
    {
        _dBoundaries = _dHandles = _dShells = 1;
    }

    return i;
}

int multiSplitEdge(AdvancedMesh *inputTriangulation, AdvancedEdge *inputEdge)
{
    List splitVertices;
    List triangles, toUnmark;

    MARK_BIT(inputEdge, 5);

    if (inputEdge->t1 != nullptr)
    {
        triangles.appendTail(inputEdge->t1);
        MARK_BIT(inputEdge->t1, 5);
    }

    if (inputEdge->t2 != nullptr)
    {
        triangles.appendTail(inputEdge->t2);
        MARK_BIT(inputEdge->t2, 5);
    }

    AdvancedTriangle *triangle;
    AdvancedTriangle *oppositrTriangle;
    AdvancedVertex *vertex;

    while ((triangle = (AdvancedTriangle *)triangles.popHead()) != nullptr)
    {
        toUnmark.appendHead(triangle);

        int numberFrontEdges = 0;

        if (IS_BIT(triangle->edge1, 5))
            numberFrontEdges++;

        if (IS_BIT(triangle->edge2, 5))
            numberFrontEdges++;

        if (IS_BIT(triangle->edge3, 5))
            numberFrontEdges++;

        if (numberFrontEdges == 3)
            continue;

        if (numberFrontEdges == 1)
        {
            AdvancedEdge *edge = (IS_BIT(triangle->edge1, 5)) ? (triangle->edge1) :
                                ((IS_BIT(triangle->edge2, 5)) ? (triangle->edge2) : (triangle->edge3));
            vertex = triangle->oppositeVertex(edge);

            if (!vertex->exactMisalignment(inputEdge->v1, inputEdge->v2))
            {
                if (!IS_BIT(vertex, 5) &&
                    AdvancedPoint::pointInInnerSegment(vertex, inputEdge->v1, inputEdge->v2))
                {
                    splitVertices.appendTail(vertex);
                    MARK_BIT(vertex, 5);
                }

                oppositrTriangle = triangle->nextEdge(edge)->oppositeTriangle(triangle);

                if (oppositrTriangle != nullptr && !IS_BIT(oppositrTriangle, 5))
                {
                    triangles.appendTail(oppositrTriangle);
                    MARK_BIT(oppositrTriangle, 5);
                }

                oppositrTriangle = triangle->prevEdge(edge)->oppositeTriangle(triangle);

                if (oppositrTriangle != nullptr && !IS_BIT(oppositrTriangle, 5))
                {
                    triangles.appendTail(oppositrTriangle);
                    MARK_BIT(oppositrTriangle, 5);
                }
            }
        }

        MARK_BIT(triangle->edge1, 5);
        MARK_BIT(triangle->edge2, 5);
        MARK_BIT(triangle->edge3, 5);
    }

    Node *iNode;
    FOR_EACH_VT_TRIANGLE((&toUnmark), triangle, iNode)
    {
        UNMARK_BIT(triangle, 5);
        UNMARK_BIT(triangle->edge1, 5);
        UNMARK_BIT(triangle->edge2, 5);
        UNMARK_BIT(triangle->edge3, 5);
    }

    FOR_EACH_VV_VERTEX((&splitVertices), vertex, iNode)
    {
        UNMARK_BIT(vertex, 5);
    }

    while (splitVertices.numberElements())
    {
        double distance, minDistance = DBL_MAX;
        AdvancedVertex *currentVertex;
        FOR_EACH_VV_VERTEX((&splitVertices), vertex, iNode)
        {
            if ((distance = vertex->squaredDistance(inputEdge->v2)) < minDistance)
            {
                currentVertex = vertex;
                minDistance = distance;
            }
        }

        splitVertices.removeNode(currentVertex);
        inputTriangulation->splitEdge(inputEdge, currentVertex);
    }

    return 1;
}

int AdvancedMesh::removeDegenerateTriangles()
{
    Node *iNode;
    AdvancedTriangle *triangle;
    AdvancedEdge *edge, *e1, *e2, *e3, *e4;
    AdvancedVertex *oppositeVertex1, *oppositeVertex2, *splitvs[2];
    int numberVertices;

    List edges(_edges);
    FOR_EACH_VE_EDGE((&edges), edge, iNode)
    {
        MARK_BIT(edge, 5);
    }

    // Split caps
    while ((edge = (AdvancedEdge *)edges.popHead()) != nullptr)
    {
        UNMARK_BIT(edge, 5);

        numberVertices = 0;

        oppositeVertex1 = (edge->t1 != nullptr) ? (edge->t1->oppositeVertex(edge)) : nullptr;
        oppositeVertex2 = (edge->t2 != nullptr) ? (edge->t2->oppositeVertex(edge)) : nullptr;

        if (oppositeVertex1 != nullptr &&
                AdvancedPoint::pointInInnerSegment(oppositeVertex1, edge->v1, edge->v2))
            splitvs[numberVertices++] = oppositeVertex1;

        if (oppositeVertex2 != nullptr &&
                AdvancedPoint::pointInInnerSegment(oppositeVertex2, edge->v1, edge->v2))
            splitvs[numberVertices++] = oppositeVertex2;

        if (numberVertices == 1)
            splitvs[1] = splitvs[0];

        if (numberVertices > 1 &&
                oppositeVertex1->squaredDistance(edge->v1) > oppositeVertex2->squaredDistance(edge->v1))
        {
            splitvs[0] = oppositeVertex2;
            splitvs[1] = oppositeVertex1;
        }

        if (numberVertices)
        {
            e1 = (edge->t1 != nullptr) ? (edge->t1->nextEdge(edge)) : nullptr;
            e2 = (edge->t1 != nullptr) ? (edge->t1->prevEdge(edge)) : nullptr;
            e3 = (edge->t2 != nullptr) ? (edge->t2->nextEdge(edge)) : nullptr;
            e4 = (edge->t2 != nullptr) ? (edge->t2->prevEdge(edge)) : nullptr;

            splitEdge(edge, splitvs[1]);

            if (numberVertices > 1)
                splitEdge(edge, splitvs[0]);

            edge = e1;
            if (edge != nullptr && !IS_BIT(edge, 5))
            {
                edges.appendTail(edge);
                MARK_BIT(edge, 5);
            }

            edge = e2;
            if (edge != nullptr && !IS_BIT(edge, 5))
            {
                edges.appendTail(edge);
                MARK_BIT(edge, 5);
            }

            edge = e3;
            if (edge != nullptr && !IS_BIT(edge, 5))
            {
                edges.appendTail(edge);
                MARK_BIT(edge, 5);
            }

            edge = e4;
            if (edge != nullptr && !IS_BIT(edge, 5))
            {
                edges.appendTail(edge);
                MARK_BIT(edge, 5);
            }
        }
    }

    // Num of collapses to remove needles
    int numberCollapses = 0;

    // Remove needles
    FOR_EACH_EDGE(edge, iNode)
    {
        if (edge->isLinked() && ((*edge->v1) == (*edge->v2)))
            if (edge->collapse())
                numberCollapses++;
    }

    FOR_EACH_EDGE(edge, iNode)
    {
        if (edge->isLinked() && ((*edge->v1) == (*edge->v2)))
        {
            if (edge->t1)
                unlinkTriangle(edge->t1);

            if (edge->t2)
                unlinkTriangle(edge->t2);
        }
    }
    removeUnlinkedElements();

    int degeneracies = 0;
    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        if (triangle->isExactlyDegenerate())
            degeneracies++;
    }

    if (degeneracies)
    {
        LOG_WARNING("removeDegenerateTriangles(): This should NOT happen!");

        FOR_EACH_TRIANGLE(triangle, iNode)
        {
            if (triangle->isExactlyDegenerate())
                MARK_VISIT(triangle);
            else
                UNMARK_VISIT(triangle);
        }
    }

    return (numberCollapses) * ((degeneracies) ? (-1) : (1));
}

bool AdvancedMesh::strongDegeneracyRemoval(int maxIterations)
{
    // Counte the number of terations
    int iterationCounter = 0;

    while ((++iterationCounter) <= maxIterations && removeDegenerateTriangles() < 0)
    {
        for (int n = 1; n < iterationCounter; ++n)
            growSelection();

        removeSelectedTriangles();

        removeSmallestComponents();

        fillHoles(_edges.numberElements(), false);
        coordBackApproximation();
    }

    if (iterationCounter > maxIterations)
        return false;
    return true;
}

int AdvancedMesh::removeSmallestComponents()
{
    /// If the mesh is composed of more than one connected partition, we will keep only the largest
    /// partition and remove all the other partitions.
    Node *iNode;
    List toProcess;
    List components;

    // Largest component (or partition) stays in the mesh
    List *component, *largestComponent = nullptr;
    AdvancedTriangle *triangle, *t1, *t2, *t3;

    // Unmark all
    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        UNMARK_BIT(triangle, 5);
    }

    triangle = ((AdvancedTriangle *)_triangles.head()->data);
    iNode = _triangles.head();

    do
    {
        component = new List;
        components.appendHead(component);
        toProcess.appendHead(triangle);

        while (toProcess.numberElements())
        {
            triangle = (AdvancedTriangle *)toProcess.head()->data;
            toProcess.removeCell(toProcess.head());

            if (!IS_BIT(triangle, 5))
            {
                t1 = triangle->t1();
                t2 = triangle->t2();
                t3 = triangle->t3();

                if (t1 != nullptr && !IS_BIT(t1, 5))
                    toProcess.appendHead(t1);

                if (t2 != nullptr && !IS_BIT(t2, 5))
                    toProcess.appendHead(t2);

                if (t3 != nullptr && !IS_BIT(t3, 5))
                    toProcess.appendHead(t3);

                MARK_BIT(triangle, 5);
                component->appendTail(triangle);
            }
        }
        toProcess.removeNodes();

        for (; iNode != nullptr; iNode = iNode->next())
        {
            triangle = ((AdvancedTriangle *) iNode->data);
            if (!IS_BIT(triangle, 5))
                break;
        }
    }
    while (iNode != nullptr);

    // Get the number of components in the mesh
    uint64_t numberComponents = components.numberElements();

    int numberElements = 0, largestNumber = 0;
    FOR_EACH_NODE(components, iNode)
    {
        if ((numberElements = ((List *)iNode->data)->numberElements()) > largestNumber)
        {
            largestNumber = numberElements;
            largestComponent = (List *)iNode->data;
        }
    }

    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        UNMARK_BIT(triangle, 5);
    }

    Node *vtNode;
    numberElements = 0;
    FOR_EACH_NODE(components, iNode)
    {
        if (((List *)iNode->data) != largestComponent)
        {
            FOR_EACH_VT_TRIANGLE(((List *)iNode->data), triangle, vtNode)
            {
                if (triangle->edge1->v1 != nullptr)
                    triangle->edge1->v1->e0 = nullptr;

                if (triangle->edge1->v2 != nullptr)
                    triangle->edge1->v2->e0 = nullptr;

                if (triangle->edge2->v1 != nullptr)
                    triangle->edge2->v1->e0 = nullptr;

                if (triangle->edge2->v2 != nullptr)
                    triangle->edge2->v2->e0 = nullptr;

                if (triangle->edge3->v1 != nullptr)
                    triangle->edge3->v1->e0 = nullptr;

                if (triangle->edge3->v2 != nullptr)
                    triangle->edge3->v2->e0 = nullptr;

                triangle->edge1->v1 = triangle->edge1->v2 = nullptr;
                triangle->edge2->v1 = triangle->edge2->v2 = nullptr;
                triangle->edge3->v1 = triangle->edge3->v2 = nullptr;
                triangle->edge1 = triangle->edge2 = triangle->edge3 = nullptr;
                numberElements++;
            }
        }
    }

    FOR_EACH_NODE(components, iNode)
    {
        delete((List *)iNode->data);
    }

    if (numberElements)
    {
        _dBoundaries = _dHandles = _dShells = 1;
        removeUnlinkedElements();
        return numberComponents - 1;
    }

    return 0;
}

int AdvancedMesh::removeSmallestComponents(double epsilonArea)
{
    // Generic
    Node *iNode;
    List toProcess, component;
    AdvancedTriangle *triangle, *adjacentTriangle;

    // Number of removed components
    int numberRemovedComponents=0;

    // If no triangles, return
    if (_triangles.numberElements() == 0)
        return 0;

    // Unmark all triangles
    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        UNMARK_BIT(triangle, 5);
    }

    // First triangle
    triangle = ((AdvancedTriangle *)_triangles.head()->data);

    // First node
    iNode = _triangles.head();

    double area;

    // Fill the toProcess list
    do
    {
        toProcess.appendTail(triangle);
        MARK_BIT(triangle, 5);

        // Accumulate
        area = 0.0;

        while ((triangle=(AdvancedTriangle *) toProcess.popHead())!=nullptr)
        {
            // Adjacent triangle (1)
            adjacentTriangle = triangle->t1();
            if (adjacentTriangle != nullptr && !IS_BIT(adjacentTriangle, 5))
            {
                toProcess.appendTail(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }

            // Adjacent triangle (2)
            adjacentTriangle = triangle->t2();
            if (adjacentTriangle != nullptr && !IS_BIT(adjacentTriangle, 5))
            {
                toProcess.appendTail(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }

            // Adjacent triangle (3)
            adjacentTriangle = triangle->t3();
            if (adjacentTriangle != nullptr && !IS_BIT(adjacentTriangle, 5))
            {
                toProcess.appendTail(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }
            component.appendTail(triangle);

            area += triangle->area();
        }

        // If the area is smaller than the defined one, remove
        if (area < epsilonArea)
        {
            numberRemovedComponents++;
            while ((triangle = (AdvancedTriangle *)component.popHead()) != nullptr)
                unlinkTriangle(triangle);
        }
        else
            component.removeNodes();

        for (; iNode != nullptr; iNode = iNode->next())
        {
            triangle = ((AdvancedTriangle *)iNode->data);

            if (!IS_BIT(triangle, 5))
                break;
        }
    } while (iNode != nullptr);

    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        UNMARK_BIT(triangle, 5);
    }

    if (numberRemovedComponents)
    {
        _dBoundaries = _dHandles = _dShells = 1;
        removeUnlinkedElements();
    }

    return numberRemovedComponents;
}

int AdvancedMesh::forceNormalConsistence()
{
    // Generic
    Node *node;
    AdvancedTriangle *triangle;

    int orientations = 0;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        if (!IS_BIT(triangle, 5))
        {
            orientations |= forceNormalConsistence(triangle);
        }
    }

    FOR_EACH_TRIANGLE(triangle, node)
    {
        UNMARK_BIT(triangle, 5);
    }

    return orientations;
}

int AdvancedMesh::forceNormalConsistence(AdvancedTriangle *inputTriangle)
{
    // Generic
    Node *node;
    AdvancedEdge *edge;
    List toProcess, edgeList;
    AdvancedTriangle *triangle;
    AdvancedTriangle *t1, *t2, *t3;

    int returnValue = 0;
    toProcess.appendHead(inputTriangle);
    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle *) toProcess.head()->data;
        toProcess.removeCell(toProcess.head());
        if (!IS_BIT(triangle, 5))
        {
            t1 = triangle->t1();
            t2 = triangle->t2();
            t3 = triangle->t3();

            if (!IS_BIT(triangle->edge1, 5))
            {
                MARK_BIT(triangle->edge1, 5);
                edgeList.appendHead(triangle->edge1);
            }

            if (!IS_BIT(triangle->edge2, 5))
            {
                MARK_BIT(triangle->edge2, 5);
                edgeList.appendHead(triangle->edge2);
            }

            if (!IS_BIT(triangle->edge3, 5))
            {
                MARK_BIT(triangle->edge3, 5);
                edgeList.appendHead(triangle->edge3);
            }

            if (t1 != nullptr && !IS_BIT(t1, 5))
            {
                toProcess.appendHead(t1);
                if (!triangle->checkAdjNor(t1))
                {
                    t1->invert();
                    returnValue = 1;
                }
            }

            if (t2 != nullptr && !IS_BIT(t2, 5))
            {
                toProcess.appendHead(t2);
                if (!triangle->checkAdjNor(t2))
                {
                    t2->invert();
                    returnValue = 1;
                }
            }

            if (t3 != nullptr && !IS_BIT(t3, 5))
            {
                toProcess.appendHead(t3);
                if (!triangle->checkAdjNor(t3))
                {
                    t3->invert();
                    returnValue = 1;
                }
            }

            MARK_BIT(triangle, 5);
        }
    }


    int warning = 0;
    int isClosed = 1;
    FOR_EACH_VE_EDGE((&(edgeList)), edge, node)
    {
        UNMARK_BIT(edge, 5);

        if (isClosed && edge->isOnBoundary())
            isClosed = 0;

        const int factor1 = (edge->t1 != nullptr) ?
              ((edge->commonVertex(edge->t1->nextEdge(edge)) == edge->v1) ? (-1) : (1)) : (0);

        const int factor2 = (edge->t2 != nullptr) ?
                ((edge->commonVertex(edge->t2->nextEdge(edge)) == edge->v2) ? (-1) : (1)) : (0);

        if (factor1 * factor2 < 0)
        {
            warning++;

            if (factor1 == -1)
            {
                swapPointers((void **)(&(edge->v1)), (void **)(&(edge->v2)));
            }

            AdvancedEdge *nEdge = newEdge(edge->v2, edge->v1);
            _edges.appendHead(nEdge);
            edge->t2->replaceEdge(edge, nEdge);
            nEdge->t2 = edge->t2;
            edge->t2 = nullptr;
        }
        else if (factor1 == -1 || factor2 == -1)
        {
            swapPointers((void **)(&(edge->v1)), (void **)(&(edge->v2)));
        }
    }

    if (warning)
    {
        _dBoundaries = _dHandles = _dShells = 1;
        LOG_WARNING("forceNormalConsistence: The mesh was NOT orientable. Cut performed.");
    }

    // Though useful in some easy cases, the flip below destroys the orientation
    // when it is set on purpose (e.g. for internal cavities)

    // if (isclosed)
    // {
    //      t = topTriangle(t0);
    //      if (t->getNormal().z < 0) {flipNormals(t0); r=1;}
    // }

    if (warning)
        returnValue |= 2;

    return returnValue;
}

int AdvancedMesh::removeOverlappingTriangles()
{
    // Generic
    Node *node;
    AdvancedEdge *edge;
    List overlappingEdges;

    FOR_EACH_EDGE(edge, node)
    {
        if (edge->overlaps())
            overlappingEdges.appendHead(edge);
    }
    overlappingEdges.sort(edgeCompare);

    // Check the overlapping edges
    for (node = overlappingEdges.tail(); node != nullptr; node=node->prev())
    {
        edge = (AdvancedEdge *) node->data;
        if (edge->overlaps() && edge->swap())
        {
            if (edge->t1->isExactlyDegenerate() ||
                edge->t2->isExactlyDegenerate())
            {
                edge->swap(1);
                continue;
            }

            if (edge->t1->nextEdge(edge)->overlaps())
            {
                edge->swap(1);
                continue;
            }

            if (edge->t1->prevEdge(edge)->overlaps())
            {
                edge->swap(1);
                continue;
            }

            if (edge->t2->nextEdge(edge)->overlaps())
            {
                edge->swap(1);
                continue;
            }

            if (edge->t2->prevEdge(edge)->overlaps())
            {
                edge->swap(1);
                continue;
            }
        }
    }

    // unlink the triangles that contained the removed edges
    int numberRemovedEdges = 0;
    for (node = overlappingEdges.tail(); node != nullptr; node=node->prev())
    {
        edge = (AdvancedEdge *) node->data;
        if (edge->overlaps())
        {
            unlinkTriangle(edge->t1);
            unlinkTriangle(edge->t2);
            numberRemovedEdges++;
        }
    }
    if (numberRemovedEdges)
    {
        // Remove the unlinked elements
        removeUnlinkedElements();

        // Update the status
        _dBoundaries = _dHandles = _dShells = 1;
    }

    // Return the number of removed edges
    return numberRemovedEdges * 2;
}

int AdvancedMesh::removeBoundaryTriangles()
{
    LOG_STATUS("Cleaning Boundary Triangles");
    TIMER_SET;

    // Deselect all the triangles
    deselectTriangles();

    // Select the boundary triangles only
    int numberBoundaryTriangles = selectBoundaryTriangles();

    // Remove the boundary triangles
    removeSelectedTriangles();

    LOG_STATUS_IMPORTANT("Cleaning Boundary Triangles Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    if (numberBoundaryTriangles)
        LOG_WARNING("The mesh has [%d] boundary triangles. DIRTY MESH!",
                    numberBoundaryTriangles);
    else
        LOG_SUCCESS("No Boundary Edges Detected");

    return numberBoundaryTriangles;
}

bool AdvancedMesh::ensureWatertightness()
{
    LOG_TITLE("Ensuring Watertightness");
    TIMER_SET;

    // Filling holes with boundary edges
    fillHoles();

    // Run geometry correction
    const int numberBoundaryEdges = boundaries();
    if (numberBoundaryEdges > 0)
    {
        LOG_WARNING("The mesh has [ %d ] boundary edge(s)", numberBoundaryEdges);
    }

    const bool isMeshClean = cleanMesh();
    if (!isMeshClean)
    {
        LOG_WARNING("Repairing the mesh was not possible at this stage!");
    }

    LOG_STATUS_IMPORTANT("Ensuring Watertightness Stats.");
    LOG_STATS(GET_TIME_SECONDS);

    return isMeshClean || (numberBoundaryEdges > 0);
}

bool AdvancedMesh::cleanMesh(uint64_t maxIterations, int innerLoops)
{
    LOG_STATUS("Cleaning Mesh");

    // Generic
    AdvancedTriangle *triangle;
    Node *node;

    // Deselect all the triangles
    deselectTriangles();

    // Invert the selection and select all the triangles
    invertSelection();

    bool intersectionRemovalFlag, degeneracyRemovalFlag;
    for (uint64_t n = 0; n < maxIterations; ++n)
    {
        LOG_STATUS("Removing Self Intersections [%d]", n);

        // Remove the degeneracies
        degeneracyRemovalFlag = strongDegeneracyRemoval(innerLoops);

        // Deselect the triangles
        deselectTriangles();

        // Invert the selection and select all the triangles again
        invertSelection();

        // Removing self-intersections
        intersectionRemovalFlag = strongIntersectionRemoval(innerLoops);

        // Removing boundary triangles
        int boundaryTriangles = removeBoundaryTriangles();

        // The mesh must have self intersections and degeneracies to process
        if (intersectionRemovalFlag && degeneracyRemovalFlag)
        {
            LOOP_COUNTER_SET;
            TIMER_SET;
            LOOP_STARTS("Final Processing")
            FOR_EACH_TRIANGLE(triangle, node)
            {
                LOOP_PROGRESS_FRACTION(++COUNTER, _triangles.numberElements());

                if (triangle->isExactlyDegenerate())
                {
                    intersectionRemovalFlag = false;
                }
            }
            LOOP_DONE;
            LOG_STATS(GET_TIME_SECONDS);

            // If all the intersections are removed, then return. Process DONE
            if (intersectionRemovalFlag && boundaryTriangles == 0)
            {
                return true;
            }
        }
    }

    // The cleaning process was unable to clean the mesh and remove the self intersections
    return false;
}

}
