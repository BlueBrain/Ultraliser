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
#include <common/Common.h>
#include <data/meshes/simple/MeshStatistics.h>
#include <data/meshes/advanced/Defines.hh>
#include <data/meshes/simple/primitives/Primitives.h>
#include <utilities/Timer.h>
#include <utilities/Data.h>

namespace Ultraliser
{

AdvancedVertex*	AdvancedMesh::newVertex()
{
    return new AdvancedVertex();
}

AdvancedVertex*	AdvancedMesh::newVertex(const double &x, const double &y, const double &z)
{
    return new AdvancedVertex(x, y, z);
}

AdvancedVertex*	AdvancedMesh::newVertex(AdvancedPoint *p)
{
    return new AdvancedVertex(p);
}

AdvancedVertex*	AdvancedMesh::newVertex(AdvancedPoint &p)
{
    return new AdvancedVertex(p);
}

AdvancedVertex*	AdvancedMesh::newVertex(AdvancedVertex* v)
{
    return new AdvancedVertex(v);
}

AdvancedEdge* AdvancedMesh::newEdge(AdvancedVertex*s, AdvancedVertex*d)
{
    return new AdvancedEdge(s, d);
}

AdvancedEdge* AdvancedMesh::newEdge(AdvancedEdge* e)
{
    return new AdvancedEdge(e->v1,e->v2);
}

AdvancedTriangle* AdvancedMesh::newTriangle()
{
    return new AdvancedTriangle();
}

AdvancedTriangle* AdvancedMesh::newTriangle(AdvancedEdge* a, AdvancedEdge* b, AdvancedEdge* c)
{
    return new AdvancedTriangle(a, b, c);
}

AdvancedMesh::AdvancedMesh()
{
    info = nullptr;

    _nBoundaries = _nHandles = _nShells = 0;
    _dBoundaries = _dHandles = _dShells = 0;
}

AdvancedMesh::AdvancedMesh(const std::string filePath)
{
    importMesh(filePath.c_str());
}

AdvancedMesh::AdvancedMesh(Vertices vertices, Triangles triangles)
{
    LOG_TITLE("Building an Advanced Mesh");
    LOG_STATUS("Cloning");

    Ultraliser::Utilities::Timer statsTimer;
    statsTimer.start();

    // Generic
    AdvancedVertex* vertex;
    Node* node;

    // Fill the _vertices list
    for(uint64_t i = 0; i < vertices.size(); ++i)
        _vertices.appendTail(
                newVertex(vertices[i].x(), vertices[i].y(), vertices[i].z()));

    // TODO: Use OpenMP to fill and delete
    int numberVertices = _vertices.numberElements();
    ExtendedVertex** vertexList = nullptr;
    vertexList = (ExtendedVertex**) malloc(
                sizeof(ExtendedVertex* )*numberVertices);

    TIMER_SET;
    LOOP_COUNTER_SET;
    LOOP_STARTS("Creating Vertex List");
    FOR_EACH_VERTEX(vertex, node)
    {
        LONG_LOOP_PROGRESS(COUNTER, numberVertices);

        // Add the vertex to the vertex list
        vertexList[COUNTER++] = new ExtendedVertex(vertex);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Set the faces
    TIMER_RESET;
    LOOP_STARTS("Creating Face List");
    for(uint64_t i = 0; i < triangles.size(); ++i)
    {
        LONG_LOOP_PROGRESS(i, triangles.size());

        Ultraliser::Triangle triangle = triangles[i];
        if (createIndexedTriangle(vertexList, triangle[0], triangle[1], triangle[2])) { /* NOTHING */ }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Remove the extended vertices
    TIMER_RESET;
    LOOP_STARTS("Cleaning Extended Data");
    if (vertexList != nullptr)
    {
        for (int i = 0; i < numberVertices; ++i)
            delete(vertexList[i]);
        free(vertexList);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Fix the connectivity, to reconstruct the topology
    fixConnectivity();

    LOG_STATUS_IMPORTANT("Building an Advanced Mesh Stats.");
    LOG_STATS(statsTimer.elapsedTimeInSeconds());
}

AdvancedMesh::AdvancedMesh(const Vertex *vertices,
                           const uint64_t &numberVertices,
                           const Ultraliser::Triangle *triangles,
                           const uint64_t &numberTriangles)
{
    LOG_TITLE("Building an Advanced Mesh");
    LOG_STATUS("Cloning");

    Ultraliser::Utilities::Timer statsTimer;
    statsTimer.start();

    // Generic
    AdvancedVertex* vertex;
    Node* node;

    // Fill the _vertices list
    for(uint64_t i = 0; i < numberVertices; ++i)
        _vertices.appendTail(newVertex(vertices[i].x(), vertices[i].y(), vertices[i].z()));

    // TODO: Use OpenMP to fill and delete
    ExtendedVertex** vertexList = nullptr;
    vertexList = (ExtendedVertex**) malloc(sizeof(ExtendedVertex*) * numberVertices);

    TIMER_SET;
    LOOP_COUNTER_SET;
    LOOP_STARTS("Creating Vertex List");
    FOR_EACH_VERTEX(vertex, node)
    {
        LONG_LOOP_PROGRESS(COUNTER, numberVertices);

        // Add the vertex to the vertex list
        vertexList[COUNTER++] = new ExtendedVertex(vertex);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Set the faces
    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Creating Face List");
    for(uint64_t i = 0; i < numberTriangles; ++i)
    {
        Ultraliser::Triangle triangle = triangles[i];
        if (createIndexedTriangle(vertexList, triangle[0], triangle[1], triangle[2])) { /* NOTHING */ }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Remove the extended vertices
    if (vertexList != nullptr)
    {
        TIMER_RESET;
        LOOP_STARTS("Cleaning Extended Data");
        for (uint64_t i = 0; i < numberVertices; ++i)
            delete(vertexList[i]);
        free(vertexList);

        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);

        // Fix the connectivity of the mesh
        fixConnectivity();

        // Update the data
        eulerUpdate();

        // _dBoundaries = _dHandles = _dShells = 1;
    }

    LOG_STATUS_IMPORTANT("Building an Advanced Mesh Stats.");
    LOG_STATS(statsTimer.elapsedTimeInSeconds());
}

AdvancedMesh::AdvancedMesh(const char *inputMeshDefinition)
{
    init(inputMeshDefinition);
}

void AdvancedMesh::init(const char *inputMeshDefinition)
{
    info = nullptr;

    if (!strcmp(inputMeshDefinition, "triangle"))
    {
        // New vertices
        AdvancedVertex* v1 = newVertex(0, 0, 0);
        AdvancedVertex* v2 = newVertex(2, 0, 0);
        AdvancedVertex* v3 = newVertex(1, 1, 0);

        // New edges
        AdvancedEdge* e1 = newEdge(v1, v2); v1->e0 = e1;
        AdvancedEdge* e2 = newEdge(v2, v3); v2->e0 = e2;
        AdvancedEdge* e3 = newEdge(v3, v1); v3->e0 = e3;

        // New triangle
        AdvancedTriangle* t1 = newTriangle(e1, e2, e3);

        // Update edges
        e1->t1 = t1; e1->t2 = nullptr;
        e2->t1 = t1; e2->t2 = nullptr;
        e3->t1 = t1; e3->t2 = nullptr;

        // Update the vertices
        _vertices.appendHead(v1);
        _vertices.appendHead(v2);
        _vertices.appendHead(v3);

        // Update the triangle
        _triangles.appendHead(t1);

        // Update the edges
        _edges.appendHead(e1);
        _edges.appendHead(e2);
        _edges.appendHead(e3);

        // State
        _nBoundaries = 1;
        _nHandles = 0;
        _nShells = 1;

        _dBoundaries = _dHandles = _dShells = 0;
    }
    else if (!strcmp(inputMeshDefinition, "tetrahedron"))
    {
        // Vertices
        AdvancedVertex* v1 = newVertex(-1, -1.4142136, 0);
        AdvancedVertex* v2 = newVertex(-1, 1.4142136, 0);
        AdvancedVertex* v3 = newVertex( 1, 0, -1.4142136);
        AdvancedVertex* v4 = newVertex( 1, 0, 1.4142136);

        // Edges
        AdvancedEdge* e1 = newEdge(v1, v2); v1->e0 = e1;
        AdvancedEdge* e2 = newEdge(v2, v3); v2->e0 = e2;
        AdvancedEdge* e3 = newEdge(v3, v1); v3->e0 = e3;
        AdvancedEdge* e4 = newEdge(v1, v4); v4->e0 = e4;
        AdvancedEdge* e5 = newEdge(v2, v4);
        AdvancedEdge* e6 = newEdge(v3, v4);

        // Triangles
        AdvancedTriangle* t1 = newTriangle(e1, e2, e3);
        AdvancedTriangle* t2 = newTriangle(e1, e4, e5);
        AdvancedTriangle* t3 = newTriangle(e2, e5, e6);
        AdvancedTriangle* t4 = newTriangle(e3, e6, e4);

        // Edges
        e1->t1 = t1; e1->t2 = t2;
        e2->t1 = t1; e2->t2 = t3;
        e3->t1 = t1; e3->t2 = t4;
        e4->t1 = t2; e4->t2 = t4;
        e5->t1 = t3; e5->t2 = t2;
        e6->t1 = t4; e6->t2 = t3;

        // Vertices
        _vertices.appendHead(v1);
        _vertices.appendHead(v2);
        _vertices.appendHead(v3);
        _vertices.appendHead(v4);

        // Triangles
        _triangles.appendHead(t1);
        _triangles.appendHead(t2);
        _triangles.appendHead(t3);
        _triangles.appendHead(t4);

        // Edges
        _edges.appendHead(e1);
        _edges.appendHead(e2);
        _edges.appendHead(e3);
        _edges.appendHead(e4);
        _edges.appendHead(e5);
        _edges.appendHead(e6);

        // State
        _nBoundaries = 0;
        _nHandles = 0;
        _nShells = 1;
        _dBoundaries = _dHandles = _dShells = 0;
    }
    else if (!strcmp(inputMeshDefinition, "cube"))
    {
        const double coordinates[8][3] = {{0, 0, 0},
                                          {1, 0, 0},
                                          {1, 1, 0},
                                          {0, 1, 0},
                                          {0, 0, 1},
                                          {1, 0, 1},
                                          {1, 1, 1},
                                          {0, 1, 1}};

        const int triangles[12][3] = {{3, 2, 1},
                                      {3, 1, 0},
                                      {4, 5, 6},
                                      {4, 6, 7},
                                      {7, 6, 2},
                                      {7, 2, 3},
                                      {0, 1, 5},
                                      {0, 5, 4},
                                      {1, 2, 6},
                                      {1, 6, 5},
                                      {3, 0, 4},
                                      {3, 4, 7}};

        ExtendedVertex* extendedVertexv[8];
        for (int i = 0; i < 8; ++i)
        {
            AdvancedVertex* v = newVertex(coordinates[i][0], coordinates[i][1], coordinates[i][2]);
            extendedVertexv[i] = new ExtendedVertex(v);
            _vertices.appendTail(v);
        }

        for (int i = 0; i < 12; ++i)
            createIndexedTriangle(extendedVertexv,
                                  triangles[i][0], triangles[i][1], triangles[i][2]);

        for (int i = 0; i < 8; ++i)
            delete extendedVertexv[i];

        // State
        _nBoundaries = 1;
        _nHandles = 0;
        _nShells = 1;
        _dBoundaries = _dHandles = _dShells = 0;
    }
    else if (!strcmp(inputMeshDefinition, "cylinder"))
    {
        const double coordinates[8][2] =
        {{1,0}, {0.7,0.7}, {0,1}, {-0.7,0.7}, {-1,0}, {-0.7,-0.7}, {0,-1}, {0.7,-0.7}};

        ExtendedVertex*extendedVertexv[8];
        for (int i = 0; i < 16; ++i)
        {
            AdvancedVertex* v = newVertex(
                coordinates[i%8][0], coordinates[i%8][1], (i < 8) ? (-1) : (1));
            extendedVertexv[i] = new ExtendedVertex(v);
            _vertices.appendTail(v);
        }

        for (int i = 0; i < 8; ++i)
        {
            createIndexedTriangle(extendedVertexv, i, (i + 1) % 8, i + 8);
            createIndexedTriangle(extendedVertexv, i + 8, (i + 1) % 8, 8 + (i + 1) % 8);
        }
        for (int i = 0; i < 6; ++i)
        {
            createIndexedTriangle(extendedVertexv, 0, i + 2, i + 1);
            createIndexedTriangle(extendedVertexv, 8, i + 9, i + 10);
        }

        for (int i = 0; i < 8; ++i)
            delete extendedVertexv[i];

        _nBoundaries = 1;
        _nHandles = 0;
        _nShells = 1;
        _dBoundaries = _dHandles = _dShells = 0;
    }
    else
    {
        LOG_ERROR("Unknown triangulation type [ %s ]", inputMeshDefinition);
    }
}


AdvancedMesh::AdvancedMesh(const AdvancedMesh *input, const bool cloneInfo)
{
    init(input, cloneInfo);
}

void AdvancedMesh::init(const AdvancedMesh *input, const bool cloneInfo)
{
    info = nullptr;

    Node *iNode;
    AdvancedVertex *iVertex, *vertex;
    AdvancedEdge *iEdge, *edge;
    AdvancedTriangle *iTriangle, *triangle;

    int i;

    // Triangles
    i = 0;
    void **tInfo = new void *[input->_triangles.numberElements()];
    FOR_EACH_VT_TRIANGLE((&(input->_triangles)), iTriangle, iNode)
    {
        tInfo[i++] = iTriangle->info;
    }

    // Edges
    i = 0;
    void **eInfo = new void *[input->_edges.numberElements()];
    FOR_EACH_VE_EDGE((&(input->_edges)), iEdge, iNode)
    {
        eInfo[i++]=iEdge->info;
    }

    // Vertices
    i = 0;
    void **vInfo = new void *[input->_vertices.numberElements()];
    FOR_EACH_VV_VERTEX((&(input->_vertices)), iVertex, iNode)
    {
        vInfo[i++] = iVertex->info;
    }

    FOR_EACH_VV_VERTEX((&(input->_vertices)), iVertex, iNode)
    {
        vertex = newVertex(iVertex);
        _vertices.appendTail(vertex);
        iVertex->info = vertex;
    }

    FOR_EACH_VE_EDGE((&(input->_edges)), iEdge, iNode)
    {
        edge=newEdge((AdvancedVertex*)iEdge->v1->info,
                     (AdvancedVertex*)iEdge->v2->info);

        _edges.appendTail(edge);
        iEdge->info = edge;
    }

    FOR_EACH_VT_TRIANGLE((&(input->_triangles)), iTriangle, iNode)
    {
        triangle=newTriangle((AdvancedEdge*)iTriangle->edge1->info,
                       (AdvancedEdge*)iTriangle->edge2->info,
                       (AdvancedEdge*)iTriangle->edge3->info);

        _triangles.appendTail(triangle);
        iTriangle->info = triangle;
    }

    FOR_EACH_VV_VERTEX((&(input->_vertices)), iVertex, iNode)
    {
        ((AdvancedVertex*)iVertex->info)->e0 = (AdvancedEdge*)iVertex->e0->info;
        iVertex->info = nullptr;
    }

    FOR_EACH_VE_EDGE((&(input->_edges)), iEdge, iNode)
    {
        ((AdvancedEdge*)iEdge->info)->t1 = (iEdge->t1) ?
        ((AdvancedTriangle*)iEdge->t1->info) : (nullptr);

        ((AdvancedEdge*)iEdge->info)->t2 = (iEdge->t2) ?
        ((AdvancedTriangle*)iEdge->t2->info) : (nullptr);

        iEdge->info = nullptr;
    }

    i = 0;
    FOR_EACH_VT_TRIANGLE((&(input->_triangles)), iTriangle, iNode)
    {
        iTriangle->info=tInfo[i++];
    }

    i = 0;
    FOR_EACH_VE_EDGE((&(input->_edges)), iEdge, iNode)
    {
        iEdge->info=eInfo[i++];
    }

    i = 0;
    FOR_EACH_VV_VERTEX((&(input->_vertices)), iVertex, iNode)
    {
        iVertex->info=vInfo[i++];
    }

    if (cloneInfo)
    {
        i = 0;

        FOR_EACH_TRIANGLE(iTriangle, iNode)
        {
            iTriangle->info=tInfo[i++];
        }

        i = 0;
        FOR_EACH_EDGE(iEdge, iNode)
        {
            iEdge->info=eInfo[i++];
        }

        i = 0;
        FOR_EACH_VERTEX(iVertex, iNode)
        {
            iVertex->info=vInfo[i++];
        }
    }

    delete [] tInfo;
    delete [] eInfo;
    delete [] vInfo;

    _dBoundaries = _dHandles = _dShells = 1;
}

AdvancedMesh::AdvancedMesh(const AdvancedTriangle* t0, const bool keepReferences)
{
    init(t0, keepReferences);
}

void AdvancedMesh::init(const AdvancedTriangle* t0, const bool keepReferences)
{
    info = nullptr;

    // Collecting list
    List collecting(t0);

    // Lists for the vertices, triangles and edges
    List listTriangles, listVertices, listEdges;

    // Iterator node
    Node* iNode;
    AdvancedTriangle* iTriangle, *triangle;
    AdvancedVertex* iVertex, *vertex;
    AdvancedEdge* iEdge, *edge;

    // Initially
    iTriangle=(AdvancedTriangle*) t0;
    MARK_VISIT2(iTriangle);

    while (collecting.numberElements())
    {
        iTriangle = (AdvancedTriangle*)collecting.popHead();
        listTriangles.appendHead(iTriangle);

        triangle=iTriangle->t1();
        if (triangle != nullptr && !IS_VISITED2(triangle))
        {
            MARK_VISIT2(triangle);
            collecting.appendHead(triangle);
        }

        triangle=iTriangle->t2();
        if (triangle != nullptr && !IS_VISITED2(triangle))
        {
            MARK_VISIT2(triangle);
            collecting.appendHead(triangle);
        }

        triangle=iTriangle->t3();
        if (triangle != nullptr && !IS_VISITED2(triangle))
        {
            MARK_VISIT2(triangle);
            collecting.appendHead(triangle);
        }
    }

    FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
    {
        UNMARK_VISIT2(iTriangle);

        iEdge = iTriangle->edge1;
        if (!IS_VISITED2(iEdge))
        {
            MARK_VISIT2(iEdge);
            listEdges.appendHead(iEdge);
        }

        iEdge = iTriangle->edge2;
        if (!IS_VISITED2(iEdge))
        {
            MARK_VISIT2(iEdge);
            listEdges.appendHead(iEdge);
        }

        iEdge = iTriangle->edge3;
        if (!IS_VISITED2(iEdge))
        {
            MARK_VISIT2(iEdge);
            listEdges.appendHead(iEdge);
        }

        iVertex = iTriangle->v1();
        if (!IS_VISITED2(iVertex))
        {
            MARK_VISIT2(iVertex);
            listVertices.appendHead(iVertex);
        }

        iVertex = iTriangle->v2();
        if (!IS_VISITED2(iVertex))
        {
            MARK_VISIT2(iVertex);
            listVertices.appendHead(iVertex);
        }

        iVertex = iTriangle->v3();
        if (!IS_VISITED2(iVertex))
        {
            MARK_VISIT2(iVertex);
            listVertices.appendHead(iVertex);
        }
    }

    FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
    {
        UNMARK_VISIT2(iVertex);
        vertex = newVertex(iVertex);
        _vertices.appendTail(vertex);
        iVertex->info = vertex;
    }

    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        UNMARK_VISIT2(iEdge);
        edge = newEdge((AdvancedVertex*)iEdge->v1->info,
                       (AdvancedVertex*)iEdge->v2->info);

        _edges.appendTail(edge); iEdge->info = edge;
    }

    FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
    {
        triangle = newTriangle((AdvancedEdge*)iTriangle->edge1->info,
                               (AdvancedEdge*)iTriangle->edge2->info,
                               (AdvancedEdge*)iTriangle->edge3->info);

        _triangles.appendTail(triangle);
        iTriangle->info = triangle;
    }

    FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
    {
        ((AdvancedVertex*)iVertex->info)->e0 = (AdvancedEdge*)iVertex->e0->info;
    }

    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        ((AdvancedEdge*)iEdge->info)->t1 = (iEdge->t1) ?
        ((AdvancedTriangle*)iEdge->t1->info) : (nullptr);

        ((AdvancedEdge*)iEdge->info)->t2 = (iEdge->t2) ?
        ((AdvancedTriangle*)iEdge->t2->info) : (nullptr);
    }

    if (!keepReferences)
    {
        FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
        {
            iVertex->info = nullptr;
        }

        FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
        {
            iEdge->info = nullptr;
        }

        FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
        {
            iTriangle->info = nullptr;
        }
    }

    eulerUpdate();
}

AdvancedMesh *AdvancedMesh::split()
{
    std::cout << "1 \n";
    if (_triangles.numberElements() == 0)
        return nullptr;

    std::cout << "2 \n";
    // Deselect all the triangles
    deselectTriangles();

    std::cout << "3 \n";
    // Get a list of all the triangles
    AdvancedTriangle* triangles = (AdvancedTriangle*) _triangles.head()->data;

    std::cout << "4 \n";
    // Select the connected component
    selectConnectedComponent(triangles);

    std::cout << "5 \n";
    // Create a new sub mesh from the original mesh
    AdvancedMesh* splitMesh = createSubMeshFromSelection(triangles);

    std::cout << "6 \n";
    // Remove the selected mesh from the original mesh
    removeSelectedTriangles();

    std::cout << "7 \n";
    // Return the split mesh
    return splitMesh;
}

std::vector< AdvancedMesh* > AdvancedMesh::splitPartitions()
{
    LOG_TITLE("Splitting Partitions");

    // A list of all the partitions in the mesh
    std::vector< AdvancedMesh* > partitions;

    // Get the number of partitions
    eulerUpdate();

    // Ensure that the mesh has no partitions
    while (_nShells > 1)
    {
        // Split and proceed
        printf("before \n");
        partitions.push_back(split());
        printf("after \n");

        // Update the mesh
        eulerUpdate();
    }

    // Return a list of all the partitions
    return partitions;
}

void AdvancedMesh::appendMeshes(std::vector< AdvancedMesh* > listMeshes)
{
    for (auto mesh : listMeshes)
    {
        append(mesh);
    }
}

AdvancedMesh::~AdvancedMesh()
{
    // Free the data
    _triangles.freeNodes();
    _vertices.freeNodes();
    _edges.freeNodes();
}

AdvancedEdge* AdvancedMesh::createEdge(AdvancedVertex* v1, AdvancedVertex* v2)
{
    // A new edge
    AdvancedEdge* edge;

    // Ensure a valid edge
    if ((edge = v1->getEdge(v2)) != nullptr)
        return edge;

    // Create a new edge
    edge = newEdge(v1, v2);

    // Propagate the information to the vertices
    v1->e0 = edge;
    v2->e0 = edge;

    // Append the edge to the list
    _edges.appendHead(edge);

    // Return a pointer to the edge
    return edge;
}

AdvancedEdge* AdvancedMesh::createEdge(ExtendedVertex* v1, ExtendedVertex* v2,
                        const bool check)
{
    AdvancedEdge* edge;
    Node* node;

    if (check)
    {
        FOR_EACH_VE_EDGE((&(v1->VE)), edge, node)
        {
            if (edge->oppositeVertex(v1->v) == v2->v)
                return edge;
        }
    }

    // Create the new edge
    edge = newEdge(v1->v, v2->v);

    // Propagate the information to the vertices
    if (v1->v->e0 == nullptr)
        v1->v->e0 = edge;
    if (v2->v->e0 == nullptr)
        v2->v->e0 = edge;

    // Append the data
    v1->VE.appendHead(edge);
    v2->VE.appendHead(edge);
    _edges.appendHead(edge);

    // Return a pointer to the edge
    return edge;
}

AdvancedTriangle* AdvancedMesh::createTriangle(AdvancedEdge* e1, AdvancedEdge* e2, AdvancedEdge* e3)
{
    AdvancedTriangle* triangle, **t1, **t2, **t3;

    if (e1->commonVertex(e2) == e1->v2 && e1->t1 == nullptr)
        t1 = &(e1->t1);
    else if (e1->commonVertex(e2) == e1->v1 && e1->t2 == nullptr)
        t1 = &(e1->t2);
    else
        return nullptr;

    if (e2->commonVertex(e3) == e2->v2 && e2->t1 == nullptr)
        t2 = &(e2->t1);
    else if (e2->commonVertex(e3) == e2->v1 && e2->t2 == nullptr)
        t2 = &(e2->t2);
    else
        return nullptr;

    if (e3->commonVertex(e1) == e3->v2 && e3->t1 == nullptr)
        t3 = &(e3->t1);
    else if (e3->commonVertex(e1) == e3->v1 && e3->t2 == nullptr)
        t3 = &(e3->t2);
    else
        return nullptr;

    // Create the triangle
    triangle = newTriangle(e1, e2, e3);

    // Append
    *t1 = *t2 = *t3 = triangle;
    _triangles.appendHead(triangle);

    // Mark thos triangle as visited
    MARK_VISIT(triangle);

    // Update the data
    _dBoundaries = _dHandles = _dShells = 1;

    // Return a pointer to the created triangle
    return triangle;
}

AdvancedTriangle*AdvancedMesh::createUnorientedTriangle(AdvancedEdge* e1,
                                                        AdvancedEdge* e2,
                                                        AdvancedEdge* e3)
{
    AdvancedTriangle* triangle, **t1, **t2, **t3;

    if (e1->t1 == nullptr)
        t1 = &(e1->t1);
    else if (e1->t2 == nullptr)
        t1 = &(e1->t2);
    else
        return nullptr;

    if (e2->t1 == nullptr)
        t2 = &(e2->t1);
    else if (e2->t2 == nullptr)
        t2 = &(e2->t2);
    else
        return nullptr;

    if (e3->t1 == nullptr)
        t3 = &(e3->t1);
    else if (e3->t2 == nullptr)
        t3 = &(e3->t2);
    else
        return nullptr;

    // Create the triangle
    triangle = newTriangle(e1, e2, e3);

    // Append
    *t1 = *t2 = *t3 = triangle;
    _triangles.appendHead(triangle);

    // Return a pointer to the created triangle
    return triangle;
}

AdvancedTriangle*AdvancedMesh::eulerEdgeTriangle(AdvancedEdge* e2, AdvancedEdge* e3)
{
    // Get the common vertex
    AdvancedVertex* commonVertex = e2->commonVertex(e3);

    // Get the adjacent triangle
    AdvancedTriangle* adjacentTriangle = (e2->t1 == nullptr) ? (e2->t2) : (e2->t1);

    // It must be valid
    if (commonVertex == nullptr || !e2->isOnBoundary() || !e3->isOnBoundary())
        return nullptr;

    // Create the edge
    AdvancedEdge* e1 = createEdge(e2->oppositeVertex(commonVertex),
                                  e3->oppositeVertex(commonVertex));

    // Create the triangle
    if (adjacentTriangle->nextEdge(e2)->hasVertex(commonVertex))
        return createTriangle(e1, e3, e2);

    return createTriangle(e1, e2, e3);
}

AdvancedEdge* AdvancedMesh::bridgeBoundaries(AdvancedEdge* edge1, AdvancedEdge* edge2)
{
    // Make sure that its is valid
    if (edge1 == edge2 || !edge1->isOnBoundary() || !edge2->isOnBoundary())
        return nullptr;

    AdvancedTriangle* triangle;
    AdvancedVertex* vertex = edge1->commonVertex(edge2);
    if (vertex != nullptr)
    {
        triangle = eulerEdgeTriangle(edge1, edge2);
        return edge1;
    }

    AdvancedVertex* bVertex1 = (edge1->t1) ? (edge1->v1) : (edge1->v2);
    AdvancedVertex* bVertex2 = (edge2->t1) ? (edge2->v2) : (edge2->v1);

    AdvancedVertex* obVertex1 = edge1->oppositeVertex(bVertex1);
    AdvancedVertex* obVertex2 = edge2->oppositeVertex(bVertex2);

    AdvancedEdge* je = createEdge(bVertex1, bVertex2);
    AdvancedEdge* je2 = createEdge(obVertex2, obVertex1);
    AdvancedEdge* je1 = createEdge(bVertex1, obVertex2);

    triangle = createTriangle(je, edge2, je1);
    triangle = createTriangle(je1, je2, edge1);

    return je1;
}

void AdvancedMesh::unlinkTriangle(AdvancedTriangle* triangle)
{
    // Get vertices
    AdvancedVertex* v1 = triangle->v1(), *v2 = triangle->v2(), *v3 = triangle->v3();

    // Get edges
    AdvancedEdge* e1 = triangle->edge1, *e2 = triangle->edge2, *e3 = triangle->edge3;

    int v1nm = (v1->isOnBoundary() && !e1->isOnBoundary() && !e2->isOnBoundary());
    int v2nm = (v2->isOnBoundary() && !e2->isOnBoundary() && !e3->isOnBoundary());
    int v3nm = (v3->isOnBoundary() && !e3->isOnBoundary() && ! e1->isOnBoundary());

    v1->e0 = ((e2->isOnBoundary()) ? (e1) : (e2));
    v2->e0 = ((e3->isOnBoundary()) ? (e2) : (e3));
    v3->e0 = ((e1->isOnBoundary()) ? (e3) : (e1));

    // Replace edges and vertices by nullptr
    e1->replaceTriangle(triangle, nullptr);
    e2->replaceTriangle(triangle, nullptr);
    e3->replaceTriangle(triangle, nullptr);

    if (e1->isIsolated() && e2->isIsolated())
        v1->e0 = nullptr;
    if (e2->isIsolated() && e3->isIsolated())
        v2->e0 = nullptr;
    if (e3->isIsolated() && e1->isIsolated())
        v3->e0 = nullptr;
    if (e1->isIsolated())
        e1->v1 = e1->v2 = nullptr;
    if (e2->isIsolated())
        e2->v1 = e2->v2 = nullptr;
    if (e3->isIsolated())
        e3->v1 = e3->v2 = nullptr;

    triangle->edge1 = triangle->edge2 = triangle->edge3 = nullptr;

    AdvancedVertex* vertex;
    AdvancedEdge* iEdge;
    List* incidentEdges;
    Node* iNode;

    if (v1nm)
    {
        vertex = newVertex(v1->x, v1->y, v1->z);
        vertex->e0 = v1->e0;
        incidentEdges = v1->getIncidentEdges();

        FOR_EACH_VE_EDGE(incidentEdges, iEdge, iNode)
        {
            iEdge->replaceVertex(v1, vertex);
        }

        delete(incidentEdges);
        v1->e0 = e1;
        _vertices.appendHead(vertex);
    }

    if (v2nm)
    {
        vertex = newVertex(v2->x, v2->y, v2->z);
        vertex->e0 = v2->e0;
        incidentEdges = v2->getIncidentEdges();

        FOR_EACH_VE_EDGE(incidentEdges, iEdge, iNode)
        {
            iEdge->replaceVertex(v2, vertex);
        }

        delete(incidentEdges);
        v2->e0 = e2;
        _vertices.appendHead(vertex);
    }

    if (v3nm)
    {
        vertex = newVertex(v3->x, v3->y, v3->z);
        vertex->e0 = v3->e0;
        incidentEdges = v3->getIncidentEdges();

        FOR_EACH_VE_EDGE(incidentEdges, iEdge, iNode)
        {
            iEdge->replaceVertex(v3, vertex);
        }

        delete(incidentEdges);
        v3->e0 = e3;
        _vertices.appendHead(vertex);
    }
}

void AdvancedMesh::unlinkTriangleNoManifold(AdvancedTriangle* triangle)
{
    // Get the edges of the triangles
    AdvancedEdge* edge1 = triangle->edge1;
    AdvancedEdge* edge2 = triangle->edge2;
    AdvancedEdge* edge3 = triangle->edge3;

    // Replace by nullptr
    edge1->replaceTriangle(triangle, nullptr);
    edge2->replaceTriangle(triangle, nullptr);
    edge3->replaceTriangle(triangle, nullptr);

    // If already isolated, set to nullptr
    if (edge1->isIsolated())
        edge1->v1 = edge1->v2 = nullptr;
    if (edge2->isIsolated())
        edge2->v1 = edge2->v2 = nullptr;
    if (edge3->isIsolated())
        edge3->v1 = edge3->v2 = nullptr;

    // Set all the triangles to nullptr to confirm the operation
    triangle->edge1 = triangle->edge2 = triangle->edge3 = nullptr;
}

int AdvancedMesh::removeTriangles()
{
    // Counter for the removed triangles
    int numberRemovedTriangles = 0;

    // Get the root node
    Node* node = _triangles.head();

    // Ensure that the root is not nullptr to proceed
    while (node != nullptr)
    {
        // Get the triangle
        AdvancedTriangle* triangle = (AdvancedTriangle*) node->data;

        // Get the next node
        node = node->next();

        if (triangle->edge1 == nullptr || triangle->edge2 == nullptr || triangle->edge3 == nullptr)
        {
            numberRemovedTriangles++;

            // Remove the node
            _triangles.removeCell((node != nullptr) ? (node->prev()) : _triangles.tail());

            // Delete the triangle
            delete triangle;
        }
    }

    // Update the data
    _dBoundaries = _dHandles = _dShells = 1;

    // Return the number of removed triangles
    return numberRemovedTriangles;
}

int AdvancedMesh::removeEdges()
{
    // Counter for the removed edges
    int numberRemovedEdges = 0;

    // Get the root node
    Node* node = _edges.head();

    // Ensure that the root is not nullptr to proceed
    while (node != nullptr)
    {
        // Get the edge
        AdvancedEdge* edge = (AdvancedEdge*) node->data;

        // Get the next node
        node = node->next();

        if (edge->v1 == nullptr || edge->v2 == nullptr)
        {
            numberRemovedEdges++;

            // Remove the node
            _edges.removeCell((node!=nullptr) ? (node->prev()):_edges.tail());

            // Delete the edge
            delete edge;
        }
    }

    // Update the data
    _dBoundaries = _dHandles = _dShells = 1;

    // Return the number of removed edges
    return numberRemovedEdges;
}

int AdvancedMesh::removeVertices()
{
    // Counter for the removed edges
    int numberRemovedVertices = 0;

    // Get the root node
    Node* node = _vertices.head();

    // Ensure that the root is not nullptr to proceed
    while (node != nullptr)
    {
        // Get the vertex
        AdvancedVertex* vertex = static_cast<AdvancedVertex*>(node->data);

        // Get the next node
        node = node->next();

        // If the vertex does NOT have or belong to an edge, then it should
        // be removed
        if (vertex->e0 == nullptr)
        {
            numberRemovedVertices++;

            // Remove the node
            _vertices.removeCell((node!=nullptr) ? (node->prev()):_vertices.tail());

            // Delete the vertex
            delete vertex;
        }
    }

    if (numberRemovedVertices > 0)
        LOG_WARNING("[%d] floating vertices removed from the mesh.", numberRemovedVertices);

    // Update the data
    _dBoundaries = _dHandles = _dShells = 1;

    // Return the number of removed vertices
    return numberRemovedVertices;
}

int AdvancedMesh::removeRedundantVertices()
{
    int numberRedundantVertices = 0;

    Node* node;
    AdvancedVertex* vertex;

    // Do it vertex by vertex
    FOR_EACH_VERTEX(vertex, node)
    {
        // If the veretex is redundant, remove it
        if (vertex->removeIfRedundant())
            numberRedundantVertices++;
    }

    // Delete the unlinked elements that were deleted before
    removeUnlinkedElements();

    // Return the counted vertices
    return numberRedundantVertices;
}

void AdvancedMesh::deselectTriangles()
{
    AdvancedTriangle* triangle;
    Node* node;

    // Unmark all the triangles
    TIMER_SET;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        // Unmark the triangle
        UNMARK_VISIT(triangle);
    }
}

void AdvancedMesh::removeSelectedTriangles()
{
    Node* node;
    AdvancedTriangle* triangle;

    FOR_EACH_TRIANGLE(triangle, node)
    {
        // Make sure that the triangle was selected to be removed
        if (IS_VISITED(triangle))
            unlinkTriangle(triangle);
    }

    // Remove the unlinked elements
    removeUnlinkedElements();
}

uint64_t AdvancedMesh::getNumberBoundaryEdges()
{
    uint64_t numberBoundaryEdges = 0;
    Node* node;
    AdvancedEdge* edge;

    FOR_EACH_EDGE(edge, node)
    {
        if (edge->isOnBoundary())
        {
            numberBoundaryEdges++;
        }
    }

    return numberBoundaryEdges;
}

int AdvancedMesh::selectBoundaryTriangles()
{
    Node* node;
    AdvancedEdge* edge;
    AdvancedVertex* vertex;
    AdvancedVertex *vertex1, *vertex2, *vertex3;
    AdvancedTriangle* triangle;

    int numberBoundaryTriangles = 0;

    // Mark the boundary triangles only
    TIMER_SET;
    uint64_t counter = 0;
    LOOP_STARTS("Selecting Boundary Edges");
    FOR_EACH_EDGE(edge, node)
    {
        LONG_LOOP_PROGRESS(++counter, _edges.numberElements());

        if (edge->isOnBoundary())
        {
            MARK_VISIT(edge->v1);
            MARK_VISIT(edge->v2);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    counter = 0;
    TIMER_RESET;
    LOOP_STARTS("Selecting Boundary Triangles");
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LONG_LOOP_PROGRESS(++counter, _triangles.numberElements());

        if (!IS_VISITED(triangle))
        {
            vertex1 = triangle->v1();
            vertex2 = triangle->v2();
            vertex3 = triangle->v3();

            if (IS_VISITED(vertex1) || IS_VISITED(vertex2) || IS_VISITED(vertex3))
            {
                // Mark
                MARK_VISIT(triangle);

                // Increment
                numberBoundaryTriangles++;
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Unmark
    TIMER_RESET;
    counter = 0;
    LOOP_STARTS("Unmarking Vertices");
    FOR_EACH_VERTEX(vertex, node)
    {
        LONG_LOOP_PROGRESS(++counter, _vertices.numberElements());
        UNMARK_VISIT(vertex);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Return the number of boundary triangles
    return numberBoundaryTriangles;
}

int AdvancedMesh::growSelection()
{
    LOG_STATUS("Growing Selection");

    Node* node;
    AdvancedVertex* vertex;
    AdvancedVertex* vertex1, *vertex2, *vertex3;
    AdvancedTriangle* triangle;

    int numberSelectedTriangles = 0;

    TIMER_SET;
    LOOP_STARTS("Marking Triangles");
    LOOP_COUNTER_SET;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _triangles.numberElements());

        if (IS_VISITED(triangle))
        {
            vertex1 = triangle->v1();
            vertex2 = triangle->v2();
            vertex3 = triangle->v3();

            MARK_VISIT(vertex1);
            MARK_VISIT(vertex2);
            MARK_VISIT(vertex3);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_STARTS("Selecting Triangles");
    LOOP_COUNTER_RESET;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _triangles.numberElements());

        if (!IS_VISITED(triangle))
        {
            vertex1 = triangle->v1();
            vertex2 = triangle->v2();
            vertex3 = triangle->v3();

            if (IS_VISITED(vertex1) || IS_VISITED(vertex2) || IS_VISITED(vertex3))
            {
                MARK_VISIT(triangle);
                numberSelectedTriangles++;
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_STARTS("Unmarking Vertices");
    LOOP_COUNTER_RESET;
    FOR_EACH_VERTEX(vertex, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _vertices.numberElements());

        UNMARK_VISIT(vertex);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return numberSelectedTriangles;
}

void AdvancedMesh::shrinkSelection()
{
    Node* n;
    AdvancedVertex* v, *v1, *v2, *v3;
    AdvancedTriangle* t;

    FOR_EACH_TRIANGLE(t, n)
    {
        if (!IS_VISITED(t))
        {
            v1 = t->v1();
            v2 = t->v2();
            v3 = t->v3();

            MARK_VISIT(v1);
            MARK_VISIT(v2);
            MARK_VISIT(v3);
        }
    }

    FOR_EACH_TRIANGLE(t, n)
    {
        if (IS_VISITED(t))
        {
            v1 = t->v1();
            v2 = t->v2();
            v3 = t->v3();

            if (IS_VISITED(v1) || IS_VISITED(v2) || IS_VISITED(v3))
                UNMARK_VISIT(t);
        }
    }

    FOR_EACH_VERTEX(v, n)
    {
        UNMARK_VISIT(v);
    }
}

void AdvancedMesh::invertSelection(AdvancedTriangle* inputTriangles)
{
    Node* node;
    AdvancedTriangle* triangle;

    if (inputTriangles != nullptr)
    {
        List toToggle(inputTriangles);

        // Selected triangle
        AdvancedTriangle* selected;

        bool unmark = IS_VISITED(inputTriangles);

        if (unmark)
            UNMARK_VISIT(inputTriangles);
        else
            MARK_VISIT(inputTriangles);

        while ((triangle = (AdvancedTriangle*)toToggle.popHead()) != nullptr)
        {
            if ((selected = triangle->t1()) != nullptr &&
                ((IS_VISITED(selected) && unmark) || (!IS_VISITED(selected) && !unmark)))
            {
                if (unmark)
                    UNMARK_VISIT(selected);
                else
                    MARK_VISIT(selected);
                toToggle.appendTail(selected);
            }

            if ((selected = triangle->t2()) != nullptr &&
               ((IS_VISITED(selected) && unmark) || (!IS_VISITED(selected) && !unmark)))
            {
                if (unmark)
                    UNMARK_VISIT(selected);
                else
                    MARK_VISIT(selected);
                toToggle.appendTail(selected);
            }

            if ((selected = triangle->t3()) != nullptr &&
                ((IS_VISITED(selected) && unmark) || (!IS_VISITED(selected) && !unmark)))
            {
                if (unmark)
                    UNMARK_VISIT(selected);
                else
                    MARK_VISIT(selected);
                toToggle.appendTail(selected);
            }
        }
    }
    else
    {
        TIMER_SET;
        uint64_t counter = 0;
        LOOP_STARTS("Inverting Selection");
        FOR_EACH_TRIANGLE(triangle, node)
        {
            LONG_LOOP_PROGRESS(++counter, _triangles.numberElements());

            if (IS_VISITED(triangle))
                UNMARK_VISIT(triangle);
            else
                MARK_VISIT(triangle);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void AdvancedMesh::reselectSelection(AdvancedTriangle* t0)
{
    if (!IS_VISITED(t0))
        return;

    Node* node;
    AdvancedTriangle* t, *s;
    List triList(t0);
    MARK_VISIT2(t0);

    while(triList.numberElements())
    {
        t = (AdvancedTriangle*)triList.popHead();
        if ((s = t->t1()) != nullptr && !IS_VISITED2(s) && IS_VISITED(s))
        {
            triList.appendHead(s);
            MARK_VISIT2(s);
        }

        if ((s = t->t2()) != nullptr && !IS_VISITED2(s) && IS_VISITED(s))
        {
            triList.appendHead(s);
            MARK_VISIT2(s);
        }

        if ((s = t->t3()) != nullptr && !IS_VISITED2(s) && IS_VISITED(s))
        {
            triList.appendHead(s);
            MARK_VISIT2(s);
        }
    }

    FOR_EACH_TRIANGLE(t, node)
    {
        if (!IS_VISITED2(t))
            UNMARK_VISIT(t);
        else
            UNMARK_VISIT2(t);
    }
}

AdvancedMesh *AdvancedMesh::createSubMeshFromSelection(AdvancedTriangle* selection,
                                                       bool keepReference)
{
    Node* iNode;

    // If no selection, return a nullptr
    if (selection != nullptr && !IS_VISITED(selection))
        return nullptr;

    // Create a new mesh object
    AdvancedMesh *subMeshObject = newObject();
    AdvancedVertex* iVertex,*vertex;
    AdvancedEdge* iEdge, *edge;
    AdvancedTriangle* iTriangle,*adjacentTriangle, *triangle;
    List listTriangles, listEdges, listVertices;
    List currentList;

    // Selection is not null
    if (selection != nullptr)
    {
        // Add the selection to the current list to start processing
        currentList.appendHead(selection);
        MARK_BIT(selection, 5);

        // Still elements
        while (currentList.numberElements())
        {
            // New triangle
            iTriangle = (AdvancedTriangle*) currentList.popHead();

            // Add it to the list
            listTriangles.appendHead(iTriangle);

            // Adjacent triangle (1)
            if ((adjacentTriangle = iTriangle->t1()) != nullptr &&
                !IS_BIT(adjacentTriangle, 5) && IS_VISITED(adjacentTriangle))
            {
                currentList.appendHead(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }

            // Adjacent triangle (2)
            if ((adjacentTriangle = iTriangle->t2()) != nullptr &&
                !IS_BIT(adjacentTriangle, 5) && IS_VISITED(adjacentTriangle))
            {
                currentList.appendHead(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }

            // Adjacent triangle (3)
            if ((adjacentTriangle = iTriangle->t3()) != nullptr &&
                !IS_BIT(adjacentTriangle, 5) && IS_VISITED(adjacentTriangle))
            {
                currentList.appendHead(adjacentTriangle);
                MARK_BIT(adjacentTriangle, 5);
            }
        }

        // Unmakrs edges
        FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
        {
            UNMARK_BIT(iTriangle->edge1, 5);
            UNMARK_BIT(iTriangle->edge2, 5);
            UNMARK_BIT(iTriangle->edge3, 5);

            UNMARK_BIT(iTriangle->v1(), 5);
            UNMARK_BIT(iTriangle->v2(), 5);
            UNMARK_BIT(iTriangle->v3(), 5);
        }

        // Update the edge list
        FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
        {
            if (!IS_BIT(iTriangle->edge1, 5))
            {
                listEdges.appendHead(iTriangle->edge1);
                MARK_BIT(iTriangle->edge1, 5);
            }

            if (!IS_BIT(iTriangle->edge2, 5))
            {
                listEdges.appendHead(iTriangle->edge2);
                MARK_BIT(iTriangle->edge2, 5);
            }

            if (!IS_BIT(iTriangle->edge3, 5))
            {
                listEdges.appendHead(iTriangle->edge3);
                MARK_BIT(iTriangle->edge3, 5);
            }

            if ((iVertex = iTriangle->v1()) && !IS_BIT(iVertex, 5))
            {
                listVertices.appendHead(iVertex);
                MARK_BIT(iVertex, 5);
            }

            if ((iVertex = iTriangle->v2()) && !IS_BIT(iVertex, 5))
            {
                listVertices.appendHead(iVertex);
                MARK_BIT(iVertex, 5);
            }

            if ((iVertex = iTriangle->v3()) && !IS_BIT(iVertex, 5))
            {
                listVertices.appendHead(iVertex);
                MARK_BIT(iVertex, 5);
            }
        }
    }

    // Null selection
    else
    {
        FOR_EACH_EDGE(iEdge, iNode)
        {
            UNMARK_BIT(iEdge, 5);
        }

        FOR_EACH_VERTEX(iVertex, iNode)
        {
            UNMARK_BIT(iVertex, 5);
        }

        FOR_EACH_TRIANGLE(iTriangle, iNode)
        {
            if (IS_VISITED(iTriangle))
            {
                listTriangles.appendHead(iTriangle);

                MARK_BIT(iTriangle->edge1, 5);
                MARK_BIT(iTriangle->edge2, 5);
                MARK_BIT(iTriangle->edge3, 5);
            }
        }

        FOR_EACH_EDGE(iEdge, iNode)
        {
            if (IS_BIT(iEdge,5))
            {
                listEdges.appendHead(iEdge);

                MARK_BIT(iEdge->v1, 5);
                MARK_BIT(iEdge->v2, 5);
            }
        }

        FOR_EACH_VERTEX(iVertex, iNode)
        {
            if (IS_BIT(iVertex,5))
                listVertices.appendHead(iVertex);
        }
    }

    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        iEdge->v1->e0 = iEdge->v2->e0 = iEdge;
    }

    int i;
    void **vertexInfo = nullptr, **edgeInfo = nullptr, **triangleInfo = nullptr;
    if (!keepReference)
    {
        vertexInfo = new void *[listVertices.numberElements()];

        // Vertex
        i = 0;
        FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
        {
            vertexInfo[i++] = iVertex->info;
        }

        edgeInfo = new void *[listEdges.numberElements()];

        // Edge
        i = 0;
        FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
        {
            edgeInfo[i++] = iEdge->info;
        }

        // Triangle
        triangleInfo = new void *[listTriangles.numberElements()];

        i = 0;
        FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
        {
            triangleInfo[i++] = iTriangle->info;
        }
    }

    // New vertex in the sub-mesh
    FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
    {
        vertex = newVertex(iVertex);
        subMeshObject->_vertices.appendTail(vertex);
        iVertex->info = vertex;
    }

    // New edge in the sub-mesh
    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        edge = newEdge((AdvancedVertex*)iEdge->v1->info,
                       (AdvancedVertex*)iEdge->v2->info);

        subMeshObject->_edges.appendTail(edge);
        iEdge->info = edge;
    }

    // New triangle in the sub-mesh
    FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
    {
        triangle = newTriangle((AdvancedEdge*)iTriangle->edge1->info,
                               (AdvancedEdge*)iTriangle->edge2->info,
                               (AdvancedEdge*)iTriangle->edge3->info);

        subMeshObject->_triangles.appendTail(triangle);
        iTriangle->info = triangle;
        triangle->info = iTriangle;
    }

    FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
    {
        ((AdvancedVertex*)iVertex->info)->e0 = (AdvancedEdge*)iVertex->e0->info;
    }

    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        ((AdvancedEdge*)iEdge->info)->t1 = (iEdge->t1 && IS_VISITED(iEdge->t1)) ?
        ((AdvancedTriangle*)iEdge->t1->info) : (nullptr);

        ((AdvancedEdge*)iEdge->info)->t2 = (iEdge->t2 && IS_VISITED(iEdge->t2)) ?
        ((AdvancedTriangle*)iEdge->t2->info) : (nullptr);
    }

    i = 0;
    if (!keepReference)
    {
        FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
        {
            iVertex->info = vertexInfo[i++];
        }
    }

    i = 0;
    if (!keepReference)
    {
        FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
        {
            iEdge->info = edgeInfo[i++];
        }
    }

    i = 0;
    if (!keepReference)
    {
        FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
        {
            iTriangle->info = triangleInfo[i++];
        }
    }

    /// Done and unmarking selections
    // Triangles
    FOR_EACH_VT_TRIANGLE((&listTriangles), iTriangle, iNode)
    {
        UNMARK_BIT(iTriangle, 5);
    }

    // Edges
    FOR_EACH_VE_EDGE((&listEdges), iEdge, iNode)
    {
        UNMARK_BIT(iEdge, 5);
    }

    // Vertices
    FOR_EACH_VV_VERTEX((&listVertices), iVertex, iNode)
    {
        UNMARK_BIT(iVertex, 5);
    }

    // If no elements in the listTriangles
    if (!listTriangles.numberElements())
    {
        // Release the sub-mesh obejct and return nullptr
        delete(subMeshObject);
        return nullptr;
    }

    // Update the sub-mesh object and return it
    subMeshObject->duplicateNonManifoldVertices();
    subMeshObject->eulerUpdate();

    return subMeshObject;
}

AdvancedMesh *AdvancedMesh::createSubMeshFromTriangle(AdvancedTriangle* t0)
{
    // A new sub-mesh object
    AdvancedMesh *subMesh = newObject("triangle");

    // Vertices
    AdvancedVertex* v1 = (AdvancedVertex*) subMesh->_vertices.head()->data;
    AdvancedVertex* v2 = (AdvancedVertex*) subMesh->_vertices.head()->next()->data;
    AdvancedVertex* v3 = (AdvancedVertex*) subMesh->_vertices.head()->next()->next()->data;

    v1->setValue(t0->v1());
    v2->setValue(t0->v3());
    v3->setValue(t0->v2());

    ((AdvancedTriangle*) subMesh->_triangles.head()->data)->info = t0->info;

    return subMesh;
}

int AdvancedMesh::selectSphericalRegion(AdvancedTriangle* inputMesh,
                                        const double &distance,
                                        const AdvancedPoint *center)
{
    // Get a region (list of triangles) of the input mesh for a radius of distance from the center
    List *region = getRegion(inputMesh, distance, center);

    Node* iNodes;
    AdvancedTriangle* triangle;

    // Count the number of selected triangles
    int numberSelectedTriangles = 0;
    FOR_EACH_VT_TRIANGLE(region, triangle, iNodes)
    {
        MARK_VISIT(triangle);
        numberSelectedTriangles++;
    }
    delete(region);

    // Return the number of selected triangles
    return numberSelectedTriangles;
}

int AdvancedMesh::deselectSphericalRegion(AdvancedTriangle* inputMesh,
                                          const double &distance,
                                          const AdvancedPoint *center)
{
    List *region = getRegion(inputMesh, distance, center);

    Node* iNode;
    AdvancedTriangle* triangle;

    // Count the number of deselected triangles
    int numberDeselectedTriangles=0;
    FOR_EACH_VT_TRIANGLE(region, triangle, iNode)
    {
        UNMARK_VISIT(triangle);
        numberDeselectedTriangles++;
    }
    delete(region);

    return numberDeselectedTriangles;
}

void AdvancedMesh::reselectSphericalRegion(AdvancedTriangle* inputMesh,
                                           const double &distance,
                                           const AdvancedPoint *center)
{
    List *region = getRegion(inputMesh, distance, center);

    Node* iNode;
    AdvancedTriangle* triangle;

    FOR_EACH_VT_TRIANGLE(region, triangle, iNode)
    {
        MARK_VISIT2(triangle);
    }

    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        if (IS_VISITED(triangle) && !IS_VISITED2(triangle))
        {
            UNMARK_VISIT(triangle);
        }
    }

    FOR_EACH_VT_TRIANGLE(region, triangle, iNode)
    {
        UNMARK_VISIT2(triangle);
    }
    delete(region);
}

bool AdvancedMesh::retriangulateSelectedRegion()
{
    List currentList;

    Node* iNode;
    AdvancedTriangle*iTriangle;
    AdvancedPoint normal;

    FOR_EACH_TRIANGLE(iTriangle, iNode)
    {
        if (IS_VISITED(iTriangle))
        {
            currentList.appendHead(iTriangle);
            normal = normal + (iTriangle->getNormal() * iTriangle->area());
        }
    }

    if (currentList.numberElements() < 2)
    {
        LOG_WARNING("retriangulateSelectedRegion: Nothing to retriangulate!");
        return 0;
    }

    FOR_EACH_VT_TRIANGLE((&(currentList)), iTriangle, iNode)
    {
        if (iTriangle->getNormal()*normal <= 0.0)
        {
            LOG_WARNING("retriangulateSelectedRegion: Complex geometry, CANNOT retriangulate!");
            return 0;
        }
    }

    if (!isSelectionSimple(&currentList))
    {
        LOG_WARNING("retriangulateSelectedRegion: Non-simple region. CANNOT retriangulate!");
        return 0;
    }

    List *internalVertices = getRegionInternalVertices(&currentList);
    FOR_EACH_VT_TRIANGLE((&(currentList)), iTriangle, iNode)
    {
        unlinkTriangle(iTriangle);
    }

    AdvancedEdge* edge = ((AdvancedEdge*)internalVertices->head()->data);
    List *vertices = ((List *)internalVertices->head()->next()->data);
    TriangulateHole(edge, vertices);
    delete(vertices);
    delete(internalVertices);

    // Remove unlinked elements from the mesh
    removeUnlinkedElements();

    // Process has been done successfully
    return 1;
}

bool AdvancedMesh::isSelectionSimple(List *selectionList)
{
    // Empty region is not simple
    if (!selectionList->numberElements())
        return 0;

    Node* iNode;
    AdvancedTriangle* adjacentTriangle, *triangle = (AdvancedTriangle*)selectionList->head()->data;
    List bdr, top(triangle);
    MARK_VISIT2(triangle);

    uint64_t numberElements = 0;

    while (top.numberElements())
    {
        triangle = (AdvancedTriangle*) top.popHead();
        numberElements++;

        // Adjacent triangle (1)
        adjacentTriangle = triangle->t1();
        if (adjacentTriangle && IS_VISITED(adjacentTriangle) && !IS_VISITED2(adjacentTriangle))
        {
            MARK_VISIT2(adjacentTriangle);
            top.appendHead(adjacentTriangle);
        }
        else if (adjacentTriangle == nullptr)
            break;
        else if (!IS_VISITED(adjacentTriangle))
            bdr.appendHead(triangle->edge1);

        // Adjacent triangle (2)
        adjacentTriangle = triangle->t2();
        if (adjacentTriangle && IS_VISITED(adjacentTriangle) && !IS_VISITED2(adjacentTriangle))
        {
            MARK_VISIT2(adjacentTriangle);
            top.appendHead(adjacentTriangle);
        }
        else if (adjacentTriangle == nullptr)
            break;
        else if (!IS_VISITED(adjacentTriangle))
            bdr.appendHead(triangle->edge2);

        // Adjacent triangle (3)
        adjacentTriangle = triangle->t3();
        if (adjacentTriangle && IS_VISITED(adjacentTriangle) && !IS_VISITED2(adjacentTriangle))
        {
            MARK_VISIT2(adjacentTriangle);
            top.appendHead(adjacentTriangle);
        }
        else if (adjacentTriangle == nullptr)
            break;
        else if (!IS_VISITED(adjacentTriangle))
            bdr.appendHead(triangle->edge3);
    }

    FOR_EACH_VT_TRIANGLE(selectionList, triangle, iNode)
    {
        UNMARK_VISIT2(triangle);
    }

    // Mesh-boundary in selection
    if (top.numberElements())
        return 0;

    // Disconnected selection
    if (numberElements != selectionList->numberElements())
        return 0;

    AdvancedEdge* iEdge, *edge, *ge = nullptr, *e0;
    List *incidentEdges;

    FOR_EACH_VE_EDGE((&(bdr)), iEdge, iNode)
    {
        MARK_VISIT(iEdge);
    }

    uint64_t numberEdges;
    numberElements = 0;
    iEdge = e0 = (AdvancedEdge*)bdr.head()->data;
    AdvancedVertex* v = iEdge->v1;

    do
    {
        numberElements++;

        // Opposite vertex
        v = iEdge->oppositeVertex(v);

        // Incident edges
        incidentEdges = v->getIncidentEdges();

        numberEdges = 0;
        FOR_EACH_VE_EDGE(incidentEdges, edge, iNode)
        {
            if (edge!=iEdge && IS_VISITED(edge))
            {
                ge = edge;
                numberEdges++;
            }
        }

        delete(incidentEdges);

        if (numberEdges > 1)
            break;

        iEdge = ge;
    } while (iEdge != e0);

    FOR_EACH_VE_EDGE((&(bdr)), iEdge, iNode)
    {
        UNMARK_VISIT(iEdge);
    }

    // Non-simple selection
    if (numberElements != bdr.numberElements())
        return 0;

    // Indeed, a simple selection
    return 1;
}

void AdvancedMesh::unmarkEverythingButSelections()
{
    AdvancedVertex* vertex;
    AdvancedEdge* edge;
    AdvancedTriangle* triangle;
    Node* node;

    // Vertices
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->mask = 0;
    }

    // Edges
    FOR_EACH_EDGE(edge, node)
    {
        edge->mask = 0;
    }

    // Triangles
    FOR_EACH_TRIANGLE(triangle, node)
    {
        triangle->mask &= (unsigned char) 1;
    }
}

int AdvancedMesh::selectConnectedComponent(AdvancedTriangle* t0, bool sos)
{
    List toProcess;
    AdvancedTriangle* triangle;
    AdvancedTriangle *t1, *t2, *t3;

    // Count the number of selected triangles
    int numberSelectedTriangles = 0;

    toProcess.appendHead(t0);
    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle*) toProcess.popHead();
        if (!IS_VISITED(triangle))
        {
            t1 = triangle->t1();
            t2 = triangle->t2();
            t3 = triangle->t3();

            if (t1 != nullptr && !IS_VISITED(t1) && (!(sos && IS_SHARPEDGE(triangle->edge1))))
                toProcess.appendHead(t1);

            if (t2 != nullptr && !IS_VISITED(t2) && (!(sos && IS_SHARPEDGE(triangle->edge2))))
                toProcess.appendHead(t2);

            if (t3 != nullptr && !IS_VISITED(t3) && (!(sos && IS_SHARPEDGE(triangle->edge3))))
                toProcess.appendHead(t3);

            MARK_VISIT(triangle);
            numberSelectedTriangles++;
        }
    }

    return numberSelectedTriangles;
}

int AdvancedMesh::deselectConnectedComponent(AdvancedTriangle* inputMesh, bool stopOnSharp)
{
    List toProcess;
    AdvancedTriangle* triangle;
    AdvancedTriangle*t1, *t2, *t3;
    int numberComponents = 0;

    toProcess.appendHead(inputMesh);
    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle*) toProcess.popHead();
        if (IS_VISITED(triangle))
        {
            t1 = triangle->t1();
            t2 = triangle->t2();
            t3 = triangle->t3();

            if (t1 != nullptr && IS_VISITED(t1) &&
                    (!(stopOnSharp && IS_SHARPEDGE(triangle->edge1))))
                toProcess.appendHead(t1);

            if (t2 != nullptr && IS_VISITED(t2) &&
                    (!(stopOnSharp && IS_SHARPEDGE(triangle->edge2))))
                toProcess.appendHead(t2);

            if (t3 != nullptr && IS_VISITED(t3) &&
                    (!(stopOnSharp && IS_SHARPEDGE(triangle->edge3))))
                toProcess.appendHead(t3);

            UNMARK_VISIT(triangle);
            numberComponents++;
        }
    }

    return numberComponents;
}

void AdvancedMesh::append(AdvancedMesh *input)
{
    // Deselect all
    deselectTriangles();

    // A copy of the input
    AdvancedMesh copy(input);

    // Select all
    copy.invertSelection();

    // Append the data
    _vertices.joinTailList(&(copy._vertices));
    _edges.joinTailList(&(copy._edges));
    _triangles.joinTailList(&(copy._triangles));

    // Update the state
    _dBoundaries = _dHandles = _dShells = 1;
}

void AdvancedMesh::moveMeshElements(AdvancedMesh *inputMesh, bool deleteInput)
{
    _vertices.joinTailList(&(inputMesh->_vertices));
    _edges.joinTailList(&(inputMesh->_edges));
    _triangles.joinTailList(&(inputMesh->_triangles));

    _dBoundaries = _dHandles = _dShells = 1;

    if(deleteInput)
        delete inputMesh;
}

List *AdvancedMesh::getRegion(AdvancedTriangle* inputMesh,
                              const double &distance,
                              const AdvancedPoint *center)
{
    List trianglesList, *toRemove = new List;

    if (inputMesh->v1()->distance(center) > distance)
        return toRemove;

    if (inputMesh->v2()->distance(center) > distance)
        return toRemove;

    if (inputMesh->v3()->distance(center) > distance)
        return toRemove;

    AdvancedTriangle* triangle;
    Node* iNode;

    // Initially
    trianglesList.appendHead(inputMesh);
    MARK_BIT(inputMesh, 3);

    while(trianglesList.numberElements() > 0)
    {
        inputMesh = (AdvancedTriangle*)trianglesList.head()->data;
        trianglesList.removeCell(trianglesList.head());
        toRemove->appendHead(inputMesh);

        if ((triangle = inputMesh->t1()) != nullptr && !IS_BIT(triangle,3) &&
             triangle->oppositeVertex(inputMesh->edge1)->distance(center) <= distance)
        {
            trianglesList.appendHead(triangle);
            MARK_BIT(triangle, 3);
        }

        if ((triangle = inputMesh->t2()) != nullptr && !IS_BIT(triangle,3) &&
             triangle->oppositeVertex(inputMesh->edge2)->distance(center) <= distance)
        {
            trianglesList.appendHead(triangle);
            MARK_BIT(triangle,3);
        }

        if ((triangle = inputMesh->t3()) != nullptr && !IS_BIT(triangle,3) &&
             triangle->oppositeVertex(inputMesh->edge3)->distance(center) <= distance)
        {
            trianglesList.appendHead(triangle);
            MARK_BIT(triangle,3);
        }
    }

    FOR_EACH_VT_TRIANGLE(toRemove, triangle, iNode)
    {
        UNMARK_BIT(triangle, 3);
    }

    return toRemove;
}

void AdvancedMesh::removeRegion(AdvancedTriangle* inputTriangulaton,
                                const double &distance,
                                const AdvancedPoint *center)
{
    List triangleList, toRemove;
    Node* iNode;
    AdvancedTriangle* triangle;

    // Initially
    triangleList.appendHead(inputTriangulaton);
    MARK_VISIT(inputTriangulaton);

    while(triangleList.numberElements() > 0)
    {
        inputTriangulaton = (AdvancedTriangle*)triangleList.head()->data;
        triangleList.removeCell(triangleList.head());
        toRemove.appendHead(inputTriangulaton);

        if ((triangle = inputTriangulaton->t1()) != nullptr && !IS_VISITED(triangle) &&
             triangle->oppositeVertex(inputTriangulaton->edge1)->distance(center) <= distance)
        {
            triangleList.appendHead(triangle);
            MARK_VISIT(triangle);
        }

        if ((triangle = inputTriangulaton->t2()) != nullptr && !IS_VISITED(triangle) &&
             triangle->oppositeVertex(inputTriangulaton->edge2)->distance(center) <= distance)
        {
            triangleList.appendHead(triangle);
            MARK_VISIT(triangle);
        }

        if ((triangle = inputTriangulaton->t3()) != nullptr && !IS_VISITED(triangle) &&
             triangle->oppositeVertex(inputTriangulaton->edge3)->distance(center) <= distance)
        {
            triangleList.appendHead(triangle);
            MARK_VISIT(triangle);
        }
    }

    for (iNode = toRemove.tail(); iNode != nullptr; iNode=iNode->prev())
    {
        triangle = ((AdvancedTriangle*)iNode->data);
        unlinkTriangle(triangle);
    }

    removeUnlinkedElements();
}

AdvancedVertex* AdvancedMesh::nextVertexOnRegionBoundary(AdvancedVertex*vertex) const
{
    AdvancedTriangle*leftTriangle, *rightTriangle;
    Node* iNode;
    AdvancedEdge* iEdge;
    List *incidentEdges = vertex->getIncidentEdges();

    FOR_EACH_VE_EDGE(incidentEdges, iEdge, iNode)
    {
        leftTriangle = iEdge->leftTriangle(vertex);
        rightTriangle = iEdge->rightTriangle(vertex);

        if (leftTriangle != nullptr && IS_VISITED(leftTriangle) &&
           (rightTriangle == nullptr || !IS_VISITED(rightTriangle)))
        {
            delete(incidentEdges);
            return iEdge->oppositeVertex(vertex);
        }
    }
    delete(incidentEdges);

    // Nothing is obtained
    return nullptr;
}

List *AdvancedMesh::getRegionInternalVertices(List *region)
{
    List *vertexList = new List;
    List *selectionList = new List;
    AdvancedEdge* edge = nullptr;

    AdvancedTriangle* triangle, *iTriangle;
    Node* iNode;
    AdvancedVertex* v1, *v2, *v3;

    FOR_EACH_VT_TRIANGLE(region, iTriangle, iNode)
    {
        MARK_VISIT(iTriangle);
        MARK_BIT(iTriangle, 3);
    }

    FOR_EACH_VT_TRIANGLE(region, iTriangle, iNode)
    {
        if (IS_BIT(iTriangle,3))
        {
            UNMARK_BIT(iTriangle,3);

            if ((triangle = iTriangle->t1()) != nullptr && !IS_VISITED(triangle))
            {
                edge = iTriangle->edge1;
                MARK_BIT(iTriangle->edge1->v1, 3);
                MARK_BIT(iTriangle->edge1->v2, 3);
            }

            if ((triangle = iTriangle->t2()) != nullptr && !IS_VISITED(triangle))
            {
                edge = iTriangle->edge2;
                MARK_BIT(iTriangle->edge2->v1, 3);
                MARK_BIT(iTriangle->edge2->v2, 3);
            }

            if ((triangle = iTriangle->t3()) != nullptr && !IS_VISITED(triangle))
            {
                edge = iTriangle->edge3;
                MARK_BIT(iTriangle->edge3->v1, 3);
                MARK_BIT(iTriangle->edge3->v2, 3);
            }
        }
    }

    FOR_EACH_VT_TRIANGLE(region, triangle, iNode)
    {
        v1 = triangle->v1();
        v2 = triangle->v2();
        v3 = triangle->v3();

        if (!IS_BIT(v1, 3))
        {
            vertexList->appendHead(v1);
            MARK_BIT(v1, 3);
        }

        if (!IS_BIT(v2, 3))
        {
            vertexList->appendHead(v2);
            MARK_BIT(v2, 3);
        }

        if (!IS_BIT(v3, 3))
        {
            vertexList->appendHead(v3);
            MARK_BIT(v3, 3);
        }
    }

    FOR_EACH_VT_TRIANGLE(region, triangle, iNode)
    {
        v1 = triangle->v1();
        v2 = triangle->v2();
        v3 = triangle->v3();

        UNMARK_BIT(v1, 3);
        UNMARK_BIT(v2, 3);
        UNMARK_BIT(v3, 3);
    }

    selectionList->appendHead(vertexList);
    selectionList->appendHead(edge);

    return selectionList;
}

void AdvancedMesh::transformShell(AdvancedTriangle* inputTriangulation, const Matrix4x4& m)
{
    List toProcess(inputTriangulation), listTriangles, listVertices;
    AdvancedTriangle* triangle, *adjacentTriangle;
    AdvancedVertex* vertex;

    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle*)toProcess.popHead();
        listTriangles.appendHead(triangle);

        adjacentTriangle = triangle->t1();
        if (adjacentTriangle != nullptr && !IS_VISITED(adjacentTriangle))
        {
            MARK_VISIT(adjacentTriangle);
            toProcess.appendHead(adjacentTriangle);
        }

        adjacentTriangle=triangle->t2();
        if (adjacentTriangle != nullptr && !IS_VISITED(adjacentTriangle))
        {
            MARK_VISIT(adjacentTriangle);
            toProcess.appendHead(adjacentTriangle);
        }

        adjacentTriangle=triangle->t3();
        if (adjacentTriangle != nullptr && !IS_VISITED(adjacentTriangle))
        {
            MARK_VISIT(adjacentTriangle);
            toProcess.appendHead(adjacentTriangle);
        }
    }

    while (listTriangles.numberElements())
    {
        triangle = (AdvancedTriangle*)listTriangles.popHead();
        UNMARK_VISIT(triangle);
        vertex = triangle->v1();
        if (!IS_VISITED(vertex))
        {
            MARK_VISIT(vertex);
            listVertices.appendHead(vertex);
        }

        vertex = triangle->v2();
        if (!IS_VISITED(vertex))
        {
            MARK_VISIT(vertex);
            listVertices.appendHead(vertex);
        }

        vertex = triangle->v3();
        if (!IS_VISITED(vertex))
        {
            MARK_VISIT(vertex);
            listVertices.appendHead(vertex);
        }
    }

    while (listVertices.numberElements())
    {
        vertex = (AdvancedVertex*)listVertices.popHead();

        UNMARK_VISIT(vertex);
        const double x = ((*vertex)*AdvancedPoint(
            m.elements[0][0], m.elements[1][0], m.elements[2][0])) + m.elements[3][0];

        const double y = ((*vertex)*AdvancedPoint(
            m.elements[0][1], m.elements[1][1], m.elements[2][1])) + m.elements[3][1];

        const double z = ((*vertex)*AdvancedPoint(
            m.elements[0][2], m.elements[1][2], m.elements[2][2])) + m.elements[3][2];

        const double w = ((*vertex)*AdvancedPoint(
            m.elements[0][3], m.elements[1][3], m.elements[2][3])) + m.elements[3][3];

        // Normalize
        vertex->x = x / w;
        vertex->y = y / w;
        vertex->z = z / w;
    }
}

void AdvancedMesh::removeShell(AdvancedTriangle* inputTriangles)
{
    List toProcess(inputTriangles);
    AdvancedTriangle* triangle, *t1, *t2, *t3;

    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle*) toProcess.popHead();

        t1 = triangle->t1();
        t2 = triangle->t2();
        t3 = triangle->t3();

        if (t1 != nullptr && !IS_VISITED2(t1))
        {
            MARK_VISIT2(t1);
            toProcess.appendHead(t1);
        }

        if (t2 != nullptr && !IS_VISITED2(t2))
        {
            MARK_VISIT2(t2);
            toProcess.appendHead(t2);
        }
        if (t3 != nullptr && !IS_VISITED2(t3))
        {
            MARK_VISIT2(t3);
            toProcess.appendHead(t3);
        }

        // unlink the triangle
        unlinkTriangle(triangle);
    }

    // Delete all the unlinked elements
    removeUnlinkedElements();
}

void AdvancedMesh::sharpEdgeTagging(const double taggingAngle)
{
    // Generic
    Node* node;
    AdvancedEdge* edge;

    FOR_EACH_EDGE(edge, node)
    {
        if (edge->curvature() > taggingAngle)
            TAG_SHARPEDGE(edge);
        else
            UNTAG_SHARPEDGE(edge);
    }
}

void AdvancedMesh::unmarkEverything()
{
    // Generic
    Node* node;
    AdvancedVertex* vertex;
    AdvancedEdge* edge;
    AdvancedTriangle* triangle;

    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->mask = 0;
    }

    FOR_EACH_EDGE(edge, node)
    {
        edge->mask = 0;
    }

    FOR_EACH_TRIANGLE(triangle, node)
    {
        triangle->mask = 0;
    }
}

double AdvancedMesh::getBoundingBox(AdvancedPoint& pMin, AdvancedPoint& pMax) const
{
    AdvancedVertex* v; Node* n;

    pMax.x = -DBL_MAX;
    pMin.x = DBL_MAX;

    pMax.y = -DBL_MAX;
    pMin.y = DBL_MAX;

    pMax.z = -DBL_MAX;
    pMin.z = DBL_MAX;

    FOR_EACH_VERTEX(v, n)
    {
        if (v->x < pMin.x)
            pMin.x = v->x;

        if (v->x > pMax.x)
            pMax.x = v->x;

        if (v->y < pMin.y)
            pMin.y = v->y;

        if (v->y > pMax.y)
            pMax.y = v->y;

        if (v->z < pMin.z)
            pMin.z = v->z;

        if (v->z > pMax.z)
            pMax.z = v->z;
    }

    return MAX(pMax.x - pMin.x, MAX(pMax.y - pMin.y, pMax.z - pMin.z));
}

double AdvancedMesh::getBoundingBallRadius() const
{
    // Generic
    AdvancedVertex* vertex;
    Node* node;

    AdvancedPoint tc, pMin, pMax;

    double tb, bsr = NODE_TO_DOUBLE(getBoundingBox(pMin, pMax)) / 2.0f;
    AdvancedPoint centerPoint = (pMax + pMin) / 2.0f;

    FOR_EACH_VERTEX(vertex, node)
    {
        if ((tb = ((*vertex) - centerPoint).length()) > bsr)
        {
            tc = ((*vertex) - centerPoint);
            tc.normalize();
            tb = ((tb - bsr) / 2.0f);
            centerPoint = centerPoint + (tc * tb);
            bsr += tb;
        }
    }

    return bsr;
}

double AdvancedMesh::area() const
{
    // Generic
    AdvancedTriangle* triangle;
    Node* node;

    // Total mesh surface area
    double totalArea = 0.0;

    // Compute the area of each triangle
    FOR_EACH_TRIANGLE(triangle, node)
    {
        totalArea  += triangle->area();
    }

    // Return the total area
    return totalArea ;
}

double AdvancedMesh::volume() const
{
    // Generic
    AdvancedTriangle* triangle;
    Node* node;

    // Total mesh volume
    double meshVolume = 0.0;

    FOR_EACH_TRIANGLE(triangle, node)
    {
        meshVolume += NODE_TO_DOUBLE((triangle->getCenter() * triangle->getNormal())) *
                triangle->area();
    }

    return meshVolume / 3.0f;
}

void AdvancedMesh::normalize(double scaleFactor)
{
    AdvancedVertex* vertex;
    Node* node;
    AdvancedPoint pMin, pMax;

    const double mel = getBoundingBox(pMin, pMax) / scaleFactor;
    FOR_EACH_VERTEX(vertex, node)
    {
        // Shift and normalize
        vertex->setValue(((*vertex) - pMin) / mel);
    }
}

void AdvancedMesh::quantize(int nc)
{
    AdvancedVertex* vertex;
    Node* node;
    normalize(nc);

    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = double(NODE_TO_INT(vertex->x));
        vertex->y = double(NODE_TO_INT(vertex->y));
        vertex->z = double(NODE_TO_INT(vertex->z));
    }
}

void AdvancedMesh::transform(const Matrix4x4& m)
{
    Node* node;
    AdvancedVertex* vertex;
    double x, y, z, w;

    FOR_EACH_VERTEX(vertex, node)
    {
        x = ((*vertex) * AdvancedPoint(m.elements[0][0], m.elements[1][0], m.elements[2][0]))
                + m.elements[3][0];
        y = ((*vertex) * AdvancedPoint(m.elements[0][1], m.elements[1][1], m.elements[2][1]))
                + m.elements[3][1];
        z = ((*vertex) * AdvancedPoint(m.elements[0][2], m.elements[1][2], m.elements[2][2]))
                + m.elements[3][2];
        w = ((*vertex) * AdvancedPoint(m.elements[0][3], m.elements[1][3], m.elements[2][3]))
                + m.elements[3][3];

        vertex->x = x / w;
        vertex->y = y / w;
        vertex->z = z / w;
    }
}

void AdvancedMesh::translate(const AdvancedPoint& translationVector)
{
    AdvancedVertex* vertex;
    Node *iNode;

    FOR_EACH_VERTEX(vertex, iNode)
    {
        vertex->setValue((*vertex)+translationVector);
    }
}

AdvancedPoint AdvancedMesh::getCenter() const
{
    AdvancedPoint centerMass;
    AdvancedTriangle* triangle;
    Node *iNode;
    double av, tvol=0.0;

    FOR_EACH_TRIANGLE(triangle, iNode)
    {
        av = triangle->area();
        tvol += av;
        centerMass += (triangle->getCenter() * av);
    }

    return centerMass / tvol;
}

void AdvancedMesh::addNormalNoise(double ns)
{
    AdvancedVertex* v;
    Node* n;
    AdvancedPoint np;
    int i;
    double noise;
    double *xyz = (double *)malloc(sizeof(double)*_vertices.numberElements()*3);
    ns *= (getBoundingBallRadius()/100.0);

    i = 0; FOR_EACH_VERTEX(v, n)
    {
        noise = ns*(((((double)rand()))-(((double) RAND_MAX) / 2.0f)) / ((double) RAND_MAX));
        np = (*v)+((v->getNormal())*noise);

        xyz[i++] = np.x;
        xyz[i++] = np.y;
        xyz[i++] = np.z;
    }

    i = 0;
    FOR_EACH_VERTEX(v, n)
    {
        v->x = xyz[i++];
        v->y = xyz[i++];
        v->z = xyz[i++];
    }

    free(xyz);
}

bool AdvancedMesh::iterativeEdgeSwaps()
{
    Node* node;
    AdvancedEdge* edge1, *edge2;
    double l;
    int swaps = 1, totits=1;
    AdvancedPoint n1, n2, nor;
    List toswap;

    bool selection=0;
    AdvancedTriangle* t;
    FOR_EACH_TRIANGLE(t, node)
    {
        if (IS_VISITED(t))
        {
            selection=1;
            break;
        }
    }

    FOR_EACH_EDGE(edge1, node)
    {
        if (!IS_SHARPEDGE(edge1) && !edge1->isOnBoundary())
        {
            MARK_VISIT(edge1);

            if ((!selection || (IS_VISITED(edge1->t1) && IS_VISITED(edge1->t2))))
                toswap.appendTail(edge1);
        }
    }


    while (swaps && totits++ < 10)
    {
        swaps = 0; for (node=toswap.head(); node!=nullptr; )
        {
            edge1 = (AdvancedEdge*)node->data;
            if (node==toswap.tail()) {toswap.removeCell(toswap.tail()); node=nullptr;}
            else {node=node->next(); toswap.removeCell(node->prev());}
            UNMARK_VISIT(edge1);

            n1 = edge1->t1->getNormal();
            n2 = edge1->t2->getNormal();
            nor = n1 + n2;
            l = edge1->delaunayMinAngle();
            if (edge1->swap())
            {
                if (edge1->delaunayMinAngle() <= l * 1.000001 ||
                    nor * edge1->t1->getNormal() <= 0 ||
                    nor * edge1->t2->getNormal() <= 0)
                    edge1->swap(1);
                else
                {
                    swaps++;
                    edge2 = edge1->t1->nextEdge(edge1);

                    if (!IS_VISITED(edge2) && !IS_SHARPEDGE(edge2) && !edge2->isOnBoundary())
                    {
                        MARK_VISIT(edge2);
                        toswap.appendHead(edge2);
                    }

                    edge2 = edge1->t1->prevEdge(edge1);
                    if (!IS_VISITED(edge2) && !IS_SHARPEDGE(edge2) && !edge2->isOnBoundary())
                    {
                        MARK_VISIT(edge2);
                        toswap.appendHead(edge2);
                    }

                    edge2 = edge1->t2->nextEdge(edge1);
                    if (!IS_VISITED(edge2) && !IS_SHARPEDGE(edge2) && !edge2->isOnBoundary())
                    {
                        MARK_VISIT(edge2); toswap.appendHead(edge2);
                    }

                    edge2 = edge1->t2->prevEdge(edge1);
                    if (!IS_VISITED(edge2) && !IS_SHARPEDGE(edge2) && !edge2->isOnBoundary())
                    {
                        MARK_VISIT(edge2);
                        toswap.appendHead(edge2);
                    }
                }
            }

        }
    }


    FOR_EACH_EDGE(edge1, node)
    {
        UNMARK_VISIT(edge1);
    }

    if (totits >= 10)
    {
        LOG_WARNING("Optimization did not converge after 10 iterations! Stopping.");
        LOG_WARNING("You may try to run the method again.");
        return 0;
    }

    return 1;
}

void AdvancedMesh::flipNormals()
{
    // Generic data
    Node* node;
    AdvancedEdge* edge;
    AdvancedTriangle* triangle;

    FOR_EACH_TRIANGLE(triangle, node)
    {
        triangle->invert();
    }

    FOR_EACH_EDGE(edge, node)
    {
        swapPointers((void **)(&(edge->v1)), (void **)(&(edge->v2)));
    }
}

void AdvancedMesh::flipNormals(AdvancedTriangle* t0)
{
    List triangleToCheck;
    AdvancedTriangle* triangle, *t1, *t2, *t3;

    triangleToCheck.appendHead(t0);
    while (triangleToCheck.numberElements())
    {
        triangle = (AdvancedTriangle*) triangleToCheck.popHead();
        if (!IS_BIT(triangle, 6))
        {
            t1 = triangle->t1();
            t2 = triangle->t2();
            t3 = triangle->t3();

            if (t1 != nullptr && !IS_BIT(t1, 6))
                triangleToCheck.appendHead(t1);

            if (t2 != nullptr && !IS_BIT(t2, 6))
                triangleToCheck.appendHead(t2);

            if (t3 != nullptr && !IS_BIT(t3, 6))
                triangleToCheck.appendHead(t3);

            // Flip the triangle
            triangle->invert();

            if (!IS_BIT(triangle->edge1, 6))
                swapPointers((void **)(&(triangle->edge1->v1)), (void **)(&(triangle->edge1->v2)));

            if (!IS_BIT(triangle->edge2, 6))
                swapPointers((void **)(&(triangle->edge2->v1)), (void **)(&(triangle->edge2->v2)));

            if (!IS_BIT(triangle->edge3, 6))
                swapPointers((void **)(&(triangle->edge3->v1)), (void **)(&(triangle->edge3->v2)));

            MARK_BIT(triangle->edge1, 6);
            MARK_BIT(triangle->edge2, 6);
            MARK_BIT(triangle->edge3, 6);
            MARK_BIT(triangle, 6);
        }
    }

    triangleToCheck.appendHead(t0);
    while (triangleToCheck.numberElements())
    {
        triangle = (AdvancedTriangle*) triangleToCheck.popHead();

        if (IS_BIT(triangle, 6))
        {
            t1 = triangle->t1();
            t2 = triangle->t2();
            t3 = triangle->t3();

            if (t1 != nullptr && IS_BIT(t1, 6))
                triangleToCheck.appendHead(t1);

            if (t2 != nullptr && IS_BIT(t2, 6))
                triangleToCheck.appendHead(t2);

            if (t3 != nullptr && IS_BIT(t3, 6))
                triangleToCheck.appendHead(t3);

            UNMARK_BIT(triangle->edge1, 6);
            UNMARK_BIT(triangle->edge2, 6);
            UNMARK_BIT(triangle->edge3, 6);
            UNMARK_BIT(triangle, 6);
        }
    }
}

AdvancedTriangle*AdvancedMesh::topTriangle(AdvancedTriangle* t0)
{
    // Generic
    AdvancedVertex* vertex;
    AdvancedVertex* v1, *v2, *v3;
    AdvancedEdge* edge;
    Node* node;
    AdvancedTriangle* triangle;
    AdvancedTriangle *t1, *t2, *t3;

    double az, Mz = -DBL_MAX;

    List toProcess;
    List triangleList, edgeList, vertexList;

    toProcess.appendHead(t0);
    MARK_BIT(t0, 2);

    while (toProcess.numberElements())
    {
        triangle = (AdvancedTriangle*) toProcess.popHead();
        triangleList.appendHead(triangle);

        t1 = triangle->t1();
        t2 = triangle->t2();
        t3 = triangle->t3();

        v1 = triangle->v1();
        v2 = triangle->v2();
        v3 = triangle->v3();

        if (!IS_VISITED(v1))
        {
            MARK_VISIT(v1);
            vertexList.appendHead(v1);
        }

        if (!IS_VISITED(v2))
        {
            MARK_VISIT(v2);
            vertexList.appendHead(v2);
        }

        if (!IS_VISITED(v3))
        {
            MARK_VISIT(v3);
            vertexList.appendHead(v3);
        }

        if (!IS_VISITED(triangle->edge1))
        {
            MARK_VISIT(triangle->edge1);
            edgeList.appendHead(triangle->edge1);
        }

        if (!IS_VISITED(triangle->edge2))
        {
            MARK_VISIT(triangle->edge2);
            edgeList.appendHead(triangle->edge2);
        }

        if (!IS_VISITED(triangle->edge3))
        {
            MARK_VISIT(triangle->edge3);
            edgeList.appendHead(triangle->edge3);
        }

        if (t1 != nullptr && !IS_BIT(t1, 2))
        {
            MARK_BIT(t1, 2);
            toProcess.appendHead(t1);
        }

        if (t2 != nullptr && !IS_BIT(t2, 2))
        {
            MARK_BIT(t2, 2);
            toProcess.appendHead(t2);
        }

        if (t3 != nullptr && !IS_BIT(t3, 2))
        {
            MARK_BIT(t3, 2);
            toProcess.appendHead(t3);
        }
    }

    List* ve = new List;

    AdvancedVertex *hv = nullptr;
    FOR_EACH_VV_VERTEX((&(vertexList)), vertex, node)
    {
        UNMARK_VISIT(vertex);
        if ((az = vertex->z) > Mz)
        {
            Mz=az;
            hv = vertex;
        }
    }

    Mz = DBL_MAX;
    FOR_EACH_VE_EDGE((&(edgeList)), edge, node)
    {
        UNMARK_VISIT(edge);
        if (edge->hasVertex(hv) && edge->length() != 0)
            ve->appendHead(edge);
    }

    FOR_EACH_VT_TRIANGLE((&(triangleList)), triangle, node)
    {
        UNMARK_BIT(triangle, 2);
    }

    AdvancedEdge *fe = nullptr;
    FOR_EACH_VE_EDGE(ve, edge, node)
    {
        if ((az = (hv->z - edge->oppositeVertex(hv)->z)/edge->length()) < Mz)
        {
            Mz=az;
            fe = edge;
        }
    }

    delete(ve);

    if (fe == nullptr)
        fe = hv->e0;

    if (fe->t1 == nullptr || fe->t2 == nullptr)
        return nullptr;

    return (FABS(fe->t1->getNormal().z) > FABS(fe->t2->getNormal().z)) ? (fe->t1) : (fe->t2);
}

void AdvancedMesh::eulerUpdate()
{
    // Generic data
    Node* node;
    AdvancedVertex* v1, *v2;
    AdvancedEdge* edge;
    AdvancedTriangle* t1, *t2;
    List triangleList;

    // Initially
    _nBoundaries = _nShells = _nHandles = 0;

    // Starting the timer
    TIMER_SET;

    // Unmark all the triangles
    LOOP_STARTS("Unmarking Triangles");
    LOOP_COUNTER_SET;
    FOR_EACH_TRIANGLE(t1, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _triangles.numberElements());

        UNMARK_BIT(t1, 5);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Unmark all the vertices
    TIMER_RESET;
    LOOP_STARTS("Unmarking VERTICES");
    LOOP_COUNTER_RESET;
    FOR_EACH_VERTEX(v1, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _vertices.numberElements());

        UNMARK_BIT(v1, 5);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_STARTS("Collecting Triangles");
    LOOP_COUNTER_RESET;
    FOR_EACH_TRIANGLE(t1, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _triangles.numberElements());

        if (!IS_BIT(t1, 5))
        {
            _nShells++;
            triangleList.appendHead(t1);
            MARK_BIT(t1, 5);

            while (triangleList.numberElements())
            {
                t1 = (AdvancedTriangle*) triangleList.popHead();
                if ((t2 = t1->t1()) != nullptr && !IS_BIT(t2, 5))
                {
                    triangleList.appendHead(t2);
                    MARK_BIT(t2, 5);
                }

                if ((t2 = t1->t2()) != nullptr && !IS_BIT(t2, 5))
                {
                    triangleList.appendHead(t2);
                    MARK_BIT(t2, 5);
                }

                if ((t2 = t1->t3()) != nullptr && !IS_BIT(t2, 5))
                {
                    triangleList.appendHead(t2);
                    MARK_BIT(t2, 5);
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Unmarking Triangles");
    FOR_EACH_TRIANGLE(t1, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _triangles.numberElements());

        UNMARK_BIT(t1, 5);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    bool hasBoundary = false;
    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Checking Boundary Edges");
    FOR_EACH_EDGE(edge, node)
    {
        LONG_LOOP_PROGRESS(++COUNTER, _edges.numberElements());
        if (edge->isOnBoundary())
        {
            hasBoundary = true;

            MARK_BIT(edge->v1, 5);
            MARK_BIT(edge->v2, 5);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Compute the boundaries
    if (hasBoundary)
    {
        TIMER_RESET;
        LOOP_COUNTER_RESET;
        LOOP_STARTS("Checking Boundary Vertices");
        FOR_EACH_VERTEX(v1, node)
        {
            LONG_LOOP_PROGRESS(++COUNTER, _vertices.numberElements());

            if (IS_BIT(v1, 5))
            {
                _nBoundaries++;
                for (v2 = v1; IS_BIT(v2, 5); v2 = v2->nextOnBoundary())
                {
                    UNMARK_BIT(v2, 5);
                }
            }
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    if (_nShells > 1)
        LOG_WARNING("The mesh has [%d] partition(s) [INVALID]", _nShells);
    else
        LOG_SUCCESS("Mesh with 1 Partition", _nShells);

    // Update the data
    _nHandles = (_edges.numberElements() - _vertices.numberElements() -
                 _triangles.numberElements() + 2 * _nShells - _nBoundaries) / 2;
    _dBoundaries = _dHandles = _dShells = 0;
}

void AdvancedMesh::openToDisk()
{
    AdvancedTriangle* t = (AdvancedTriangle*)_triangles.head()->data;
    AdvancedTriangle*s;
    List triList, *ve;
    AdvancedVertex* v, *w;
    AdvancedEdge* e, *ne;
    Node* n;
    triList.appendHead(t);
    MARK_BIT(t,3);

    while(triList.numberElements())
    {
        t = (AdvancedTriangle*)triList.popHead();
        if ((s = t->t1()) != nullptr && !IS_BIT(s,3))
        {
            triList.appendTail(s);
            MARK_BIT(s,3);
            MARK_BIT(t->edge1,3);
        }

        if ((s = t->t2()) != nullptr && !IS_BIT(s,3))
        {
            triList.appendTail(s);
            MARK_BIT(s,3);
            MARK_BIT(t->edge2,3);
        }

        if ((s = t->t3()) != nullptr && !IS_BIT(s,3))
        {
            triList.appendTail(s);
            MARK_BIT(s,3);
            MARK_BIT(t->edge3,3);
        }
    }

    FOR_EACH_TRIANGLE (t, n)
    {
        UNMARK_BIT(t, 3);
    }

    FOR_EACH_VERTEX(v, n)
    {
        v->info = new List;
    }

    FOR_EACH_EDGE(e, n) if (!IS_BIT(e, 3))
    {
        ((List *) e->v1->info)->appendHead(e);
        ((List *) e->v2->info)->appendHead(e);
    }

    FOR_EACH_VERTEX(v, n)
    {
        if (((List *)v->info)->numberElements()==1)
            triList.appendHead(v);
    }

    if (!triList.numberElements())
    {
        LOG_ERROR("openToDisk: Can NOT find a root!");
    }

    while(triList.numberElements())
    {
        v = (AdvancedVertex*)triList.popHead();
        ve = ((List *)v->info);
        if (ve->numberElements())
        {
            e = (AdvancedEdge*)(ve->head()->data);
            MARK_BIT(e, 3);
            ve->popHead();
            w = e->oppositeVertex(v);
            ve = ((List *)w->info);
            ve->removeNode(e);

            if (ve->numberElements() == 1)
                triList.appendHead(w);
        }
        else
        {
            ve = v->getIncidentEdges();
            e = (AdvancedEdge*) ve->head()->data;
            UNMARK_BIT(e, 3);
            ((List *)v->info)->appendHead(e);

            e = (AdvancedEdge*) ve->head()->next()->data;
            UNMARK_BIT(e, 3);
            ((List *)v->info)->appendHead(e);

            delete(ve);
        }
    }

    FOR_EACH_EDGE(e, n) if (!IS_BIT(e, 3) && !e->isOnBoundary())
    {
        ne = newEdge(e->v1, e->v2);
        ne->t1 = e->t1; e->t1 = nullptr; _edges.appendHead(ne);
        ne->t1->replaceEdge(e, ne);
    }

    FOR_EACH_EDGE(e, n)
    {
        UNMARK_BIT(e, 3);
    }

    FOR_EACH_VERTEX(v, n)
    {
        if (v->info)
        {
            delete(((List *)v->info));
            v->info = nullptr;
        }
    }

    duplicateNonManifoldVertices();

    _dBoundaries = _dHandles = _dShells = 1;
}

AdvancedVertex* AdvancedMesh::splitEdge(AdvancedEdge* e, AdvancedPoint *p, bool copyMask)
{
    if ((*p) == (*(e->v1)))
        return e->v1;
    if ((*p) == (*(e->v2)))
        return e->v2;

    AdvancedVertex* v3 = (e->t1 != nullptr) ? (e->t1->oppositeVertex(e)) : (nullptr);
    AdvancedVertex* v4 = (e->t2 != nullptr) ? (e->t2->oppositeVertex(e)) : (nullptr);

    AdvancedEdge* be1 = (e->t1 != nullptr) ? (e->t1->nextEdge(e)) : (nullptr);
    AdvancedEdge* be4 = (e->t2 != nullptr) ? (e->t2->prevEdge(e)) : (nullptr);

    AdvancedVertex* vertex = newVertex(p->x, p->y, p->z);
    AdvancedEdge* ne = newEdge(vertex, e->v2);
    AdvancedEdge* ne1 = (e->t1 != nullptr) ? (newEdge(vertex, v3)) : (nullptr);
    AdvancedEdge* ne2 = (e->t2 != nullptr) ? (newEdge(vertex, v4)) : (nullptr);

    AdvancedTriangle* nt1 = (e->t1 != nullptr) ? (newTriangle(ne1, ne,be1)) : (nullptr);
    AdvancedTriangle* nt2 = (e->t2 != nullptr) ? (newTriangle(ne, ne2,be4)) : (nullptr);

    ne->t1 = nt1;
    ne->t2 = nt2;

    if (ne1 != nullptr)
    {
        ne1->t1 = e->t1;
        ne1->t2 = nt1;
    }

    if (ne2 != nullptr)
    {
        ne2->t1 = nt2;
        ne2->t2 = e->t2;
    }

    if (be1 != nullptr)
        be1->replaceTriangle(e->t1, nt1);

    if (be4 != nullptr)
        be4->replaceTriangle(e->t2, nt2);

    e->v2->e0 = (be1 != nullptr) ? (be1) : (be4);
    e->v2 = vertex;
    vertex->e0 = e;

    if (e->t1 != nullptr)
        e->t1->replaceEdge(be1, ne1);
    if (e->t2 != nullptr)
        e->t2->replaceEdge(be4, ne2);

    if (copyMask)
    {
        ne->mask = e->mask;
        if (nt1 != nullptr)
            nt1->mask = e->t1->mask;

        if (nt2 != nullptr)
            nt2->mask = e->t2->mask;
    }

    // Update the data
    _vertices.appendHead(vertex);
    _edges.appendHead(ne);

    if (ne1 != nullptr)
        _edges.appendHead(ne1);

    if (ne2 != nullptr)
        _edges.appendHead(ne2);

    if (nt1 != nullptr)
        _triangles.appendHead(nt1);

    if (nt2 != nullptr)
        _triangles.appendHead(nt2);

    return vertex;
}

AdvancedVertex* AdvancedMesh::splitTriangle(AdvancedTriangle* triangle,
                                            AdvancedPoint *point,
                                            bool copyMask)
{
    AdvancedVertex* vertex1 = triangle->v1();
    AdvancedVertex* vertex2 = triangle->v2();
    AdvancedVertex* vertex3 = triangle->v3();

    // New vertex
    AdvancedVertex* vertex = newVertex(point->x, point->y, point->z);

    // New edges
    AdvancedEdge* newEdge1 = newEdge(vertex, vertex1);
    AdvancedEdge* newEdge2 = newEdge(vertex, vertex2);
    AdvancedEdge* newEdge3 = newEdge(vertex, vertex3);

    // New triangles
    AdvancedTriangle* newTriangle1 = newTriangle(newEdge2, triangle->edge3, newEdge3);
    AdvancedTriangle* newTriangle2 = newTriangle(newEdge3, triangle->edge1, newEdge1);

    // Replace the triangles
    triangle->edge3->replaceTriangle(triangle, newTriangle1);
    triangle->edge1->replaceTriangle(triangle, newTriangle2);

    // Replace the edges
    triangle->replaceEdge(triangle->edge3, newEdge2);
    triangle->replaceEdge(triangle->edge1, newEdge1);

    // Edges
    newEdge1->t1 = triangle;
    newEdge1->t2 = newTriangle2;
    newEdge2->t1 = newTriangle1;
    newEdge2->t2 = triangle;
    newEdge3->t1 = newTriangle2;
    newEdge3->t2 = newTriangle1;

    // Propagate to the vertices
    vertex->e0 = newEdge1;

    // Update the data
    _vertices.appendHead(vertex);
    _edges.appendHead(newEdge1);
    _edges.appendHead(newEdge2);
    _edges.appendHead(newEdge3);
    _triangles.appendHead(newTriangle1);
    _triangles.appendHead(newTriangle2);

    // Copy masks
    if (copyMask)
    {
        newTriangle1->mask = triangle->mask;
        newTriangle2->mask = triangle->mask;
    }

    // Return the newly created vertex
    return vertex;
}

void AdvancedMesh::printReport()
{
    eulerUpdate();

    LOG_STATUS("Mesh Stats.");

    LOG_INFO("\t* Number Vertices       | %s",
             FORMAT(_vertices.numberElements()));
    LOG_INFO("\t* Number Triangles      | %s",
             FORMAT(_triangles.numberElements()));
    LOG_INFO("\t* Number Edges      | %s",
             FORMAT(_triangles.numberElements()));
    LOG_INFO("\t* Number Partitions      | %s",
             FORMAT(shells()));
}


// Unless there are bugs, the following should be exact ...

bool AdvancedMesh::isInnerPoint(AdvancedPoint& p) const
{
    // An ampty mesh does not enclose anything
    if (_triangles.numberElements() == 0)
        return false;

    // We assume that the mesh is correctly oriented.
    Node* iNode;
    AdvancedVertex* v1, *v2, *v3;
    AdvancedTriangle* iTriangle;

    // Ray casting along positive X direction
    AdvancedPoint p2(1, 0, 0);
    p2 += p;
    double ad, minDistance = DBL_MAX;
    AdvancedPoint ip;
    AdvancedTriangle* closestTriangle = nullptr;
    AdvancedEdge* e, *closestEdge = nullptr;
    AdvancedVertex* closestVertex = nullptr;

    double o1, o2, o3;
    FOR_EACH_TRIANGLE(iTriangle, iNode)
    {
        v1 = iTriangle->v1();
        v2 = iTriangle->v2();
        v3 = iTriangle->v3();

        if (((v1->y > p.y && v2->y > p.y) || (v1->y < p.y && v2->y < p.y)) &&
            ((v1->y > p.y && v3->y > p.y) || (v1->y < p.y && v3->y < p.y)))
            continue;

        if (((v1->z > p.z && v2->z > p.z) || (v1->z < p.z && v2->z < p.z)) &&
            ((v1->z > p.z && v3->z > p.z) || (v1->z < p.z && v3->z < p.z)))
            continue;

        o1 = orient2D(p.y, p.z, v1->y, v1->z, v2->y, v2->z);
        o2 = orient2D(p.y, p.z, v2->y, v2->z, v3->y, v3->z);
        o3 = orient2D(p.y, p.z, v3->y, v3->z, v1->y, v1->z);

        // Degenerate triangle. Skip.
        if (o1 == 0 && o2 == 0 && o3 == 0)
            continue;
        else if (o1 == 0 && o2 == 0)
        {
            // Point is on surface
            if ((ad = v2->x - p.x) == 0)
                return false;
            else if (ad > 0 && ad<minDistance)
            {
                closestVertex = v2;
                closestEdge = nullptr;
                 closestTriangle = nullptr;
                minDistance = ad;
            }
        }
        else if (o1 == 0 && o3 == 0)
        {
            // Point is on surface
            if ((ad = v1->x - p.x) == 0)
                return false;
            else if (ad > 0 && ad<minDistance)
            {
                closestVertex = v1;
                closestEdge = nullptr;
                closestTriangle = nullptr;
                minDistance = ad;
            }
        }
        else if (o2 == 0 && o3 == 0)
        {
            // Point is on surface
            if ((ad = v3->x - p.x) == 0)
                return false;
            else if (ad > 0 && ad<minDistance)
            {
                closestVertex = v3;
                closestEdge = nullptr;
                closestTriangle = nullptr;
                minDistance = ad;
            }
        }
        else if (o1 == 0 || o2 == 0 || o3 == 0)
        {
            e = (o1 == 0) ? (iTriangle->edge2) : ((o2 == 0) ? (iTriangle->edge3) : (iTriangle->edge1));

            if ((p.y < e->v1->y && p.y < e->v2->y) ||
                (p.y > e->v1->y && p.y > e->v2->y) ||
                (p.z < e->v1->z && p.z < e->v2->z) ||
                (p.z > e->v1->z && p.z > e->v2->z))
                continue;

            ip = AdvancedPoint::lineLineIntersection(p, p2, *e->v1, *e->v2);

            ad = ip.x - p.x;

            // Point is on surface
            if (ad == 0)
                return false;
            else if (ad > 0 && ad<minDistance)
            {
                closestVertex = nullptr;
                closestEdge = e;
                closestTriangle = nullptr;
                minDistance = ad;
            }
        }
        else if ((o1 > 0 && o2 > 0 && o3 > 0) || (o1 < 0 && o2 < 0 && o3 < 0))
        {
            ip = AdvancedPoint::linePlaneIntersection(p, p2, *v1, *v2, *v3);
            ad = ip.x - p.x;

            // Point is on surface
            if (ad == 0)
                return false;
            else if (ad > 0 && ad<minDistance)
            {
                closestVertex = nullptr;
                closestEdge = nullptr;
                closestTriangle = iTriangle;
                minDistance = ad;
            }
        }
    }

    if (closestVertex != nullptr)
    {
        List *incidentEdges = closestVertex->getIncidentEdges();
        double ad, mind = DBL_MAX;
        AdvancedPoint a = p - (*closestVertex);
        FOR_EACH_VE_EDGE(incidentEdges, e, iNode)
        {
            AdvancedVertex*ov1 = e->oppositeVertex(closestVertex);
            AdvancedPoint b = (*ov1) - (*closestVertex);

            ad = (((a&b)*(a&b)) / ((b*b)*(a*a))) - 1;

            if (a*b < 0)
                ad = -ad;

            if (ad < mind)
            {
                mind = ad;
                closestEdge = e;
            }
        }
        delete incidentEdges;
    }

    // If cp is in the interior of an edge, select one of its incident triangles
    if (closestEdge != nullptr)
    {
        if (closestEdge->isOnBoundary())
            return false;

        AdvancedVertex* ov1 = closestEdge->t1->oppositeVertex(closestEdge);
        AdvancedVertex* ov2 = closestEdge->t2->oppositeVertex(closestEdge);

        double o1 = p.exactOrientation(closestEdge->v1, closestEdge->v2, ov1);
        double o2 = ov2->exactOrientation(closestEdge->v1, closestEdge->v2, ov1);

        closestTriangle = ((o1 >= 0 && o2 >= 0) || (o1 <= 0 && o2 <= 0)) ?
                    (closestEdge->t2) : (closestEdge->t1);
    }

    // If closest point is in the interior of a triangle, just check orientation
    if ( closestTriangle != nullptr)
    {
        return ( closestTriangle->getVector().x > 0);
    }

    return false;
}

double AdvancedMesh::_cotangentAngle(const Vector3f& pivot, const Vector3f& a, const Vector3f& b)
{
  const Vector3f pa = (a - pivot).normalized();
  const Vector3f pb = (b - pivot).normalized();

  const double sinA = Vector3f::cross(pa, pb).abs();
  const double cosA = Vector3f::dot(pa, pb);

  return (cosA / sinA);
}

AdvancedIndexedMesh AdvancedMesh::_buildIndexedMesh(AdvancedMesh& mesh)
{
    AdvancedIndexedMesh tempMesh;
    tempMesh._triangles.reserve(mesh._triangles.numberElements());
    tempMesh._vertices.resize(mesh._vertices.numberElements());

    AdvancedVertex** vertices = (AdvancedVertex **)mesh._vertices.toArray();
    const uint64_t numVertices = mesh._vertices.numberElements();

    std::unordered_map<const AdvancedVertex*, uint64_t> vertexMap;

    for(uint64_t i = 0; i < numVertices; i++)
    {
        vertexMap[vertices[i]] = i;
        tempMesh._vertices[i] = vertices[i];
    }

    const AdvancedTriangle** triangles = (const AdvancedTriangle **)mesh._triangles.toArray();

    const uint64_t numFaces = mesh._triangles.numberElements();

    for (uint64_t i = 0; i < numFaces; i++)
    {
        const AdvancedTriangle* tri = triangles[i];
        const AdvancedVertex* v1 = tri->v1();
        const AdvancedVertex* v2 = tri->v2();
        const AdvancedVertex* v3 = tri->v3();

        auto v1it = vertexMap.find(v1);
        auto v2it = vertexMap.find(v2);
        auto v3it = vertexMap.find(v3);

        if (v1it != vertexMap.end() && v2it != vertexMap.end() && v3it != vertexMap.end())
        {
            tempMesh._triangles.emplace_back(tri, v1it->second, v2it->second, v3it->second);
        }
    }

    return tempMesh;
}

void AdvancedMesh::_computeNeighborhoods(AdvancedIndexedMesh& mesh,
                                         Neighborhood& vertexNeighbors,
                                         Neighborhood& faceNeighbors)
{
    const uint64_t numV = mesh._vertices.size();

    vertexNeighbors.resize(numV);
    faceNeighbors.resize(numV);

    for(uint64_t i = 0; i < mesh._triangles.size(); i++)
    {
        const AdvancedIndexedTriangle & face = mesh._triangles[i];

        vertexNeighbors[face._v1].insert(face._v2);
        vertexNeighbors[face._v1].insert(face._v3);

        vertexNeighbors[face._v2].insert(face._v1);
        vertexNeighbors[face._v2].insert(face._v3);

        vertexNeighbors[face._v3].insert(face._v1);
        vertexNeighbors[face._v3].insert(face._v2);

        faceNeighbors[face._v1].insert(i);
        faceNeighbors[face._v2].insert(i);
        faceNeighbors[face._v3].insert(i);
    }
}

Vector3f AdvancedMesh::_convertVertexToVector3f(const AdvancedVertex* v)
{
    return Vector3f(D2F(v->x), D2F(v->y), D2F(v->z));
}

Vector3f AdvancedMesh::_computeKernel(AdvancedIndexedMesh& mesh,
                                      const uint32_t &vertexIndex,
                                      Neighbors& verticesN,
                                      Neighbors& facesN)
{
    double sumWeights = 0.0;
    Vector3f sumPos (0.0, 0.0, 0.0);

    const Vector3f currentVertex =
            _convertVertexToVector3f(mesh._vertices[vertexIndex]);

    for(const auto & nv : verticesN)
    {
        if(nv == vertexIndex)
            continue;

        const double weight =
                _computeCotangentWeight(mesh, vertexIndex, nv, facesN);
        sumWeights += weight;
        sumPos += _convertVertexToVector3f(mesh._vertices[nv]) * weight;
    }

    const Vector3f kernel = (sumPos / sumWeights) - currentVertex;

    return kernel;
}

float AdvancedMesh::_computeCotangentWeight(AdvancedIndexedMesh& mesh,
                                            const uint32_t &vertexIndex,
                                            const uint32_t &neighborIndex,
                                            Neighbors& faceN)
{
    uint32_t e1 = 0, e2 = 0;
    bool firstFound = false, secondFound = false;
    for(const auto & nvFace : faceN)
    {
        const AdvancedIndexedTriangle& f = mesh._triangles[nvFace];
        if(f.containsVertex(vertexIndex) && f.containsVertex(neighborIndex))
        {
            if(!firstFound)
            {
                firstFound = true;
                e1 = f._v1 != vertexIndex &&
                        f._v1 != neighborIndex? f._v1 : (f._v2 != vertexIndex &&
                        f._v2 != neighborIndex? f._v2 : f._v3);
            }
            else if(!secondFound)
            {
                secondFound = true;
                e2 = f._v1 != vertexIndex &&
                        f._v1 != neighborIndex? f._v1 : (f._v2 != vertexIndex &&
                        f._v2 != neighborIndex? f._v2 : f._v3);
            }
        }

        if(firstFound && secondFound)
            break;
    }

    if(!firstFound || !secondFound)
    {
        throw std::runtime_error("LaplacianSmoothOperator: Mesh has boundaries");
    }

    const Vector3f pivotA = _convertVertexToVector3f(mesh._vertices[e1]);
    const Vector3f pivotB = _convertVertexToVector3f(mesh._vertices[e2]);

    const Vector3f a = _convertVertexToVector3f(mesh._vertices[vertexIndex]);
    const Vector3f b = _convertVertexToVector3f(mesh._vertices[neighborIndex]);

    const double ctgA = _cotangentAngle(pivotA, a, b);
    const double ctgB = _cotangentAngle(pivotB, a, b);

    return 0.5f * D2F(ctgA + ctgB);
}

AdvancedVertex AdvancedMesh::_smoothVertex(AdvancedIndexedMesh& mesh,
                     const uint32_t vertexIndex,
                     const Vector3f kernel,
                     const float &smoothParam,
                     const float &inflateParam)
{
    const AdvancedVertex* v = mesh._vertices[vertexIndex];
    Vector3f newPosition = _convertVertexToVector3f(v);
    newPosition = (1.0f - smoothParam) * newPosition + (smoothParam * kernel);
    newPosition = (1.0f - inflateParam) * newPosition + (inflateParam * kernel);

    AdvancedVertex newVertex (v);
    newVertex.setValue(F2D(newPosition.x()),
                       F2D(newPosition.y()),
                       F2D(newPosition.z()));
    return newVertex;
}

void AdvancedMesh::applyLaplacianSmooth(const uint64_t &numIterations,
                                        const float &smoothLambda,
                                        const float &inflateMu)
{
    if (numIterations < 1)
    {
        return;
    }

    AdvancedIndexedMesh mesh = _buildIndexedMesh(*this);

    Neighborhood vertexN, faceN;
    _computeNeighborhoods(mesh, vertexN, faceN);

    std::vector<AdvancedVertex> smoothedVertices(mesh._vertices.size());
    for(uint64_t i = 0; i < numIterations; ++i)
    {
        // Store vertices updates in separate vector to allow for parallel processing
        #pragma omp parallel for
        for(std::size_t v = 0; v < mesh._vertices.size(); ++v)
        {
            const Vector3f kernel = _computeKernel(mesh, v, vertexN[v], faceN[v]);
            smoothedVertices[v] = _smoothVertex(mesh, v, kernel, smoothLambda, inflateMu);
        }

        // Update vertices
        #pragma omp parallel for
        for(std::size_t v = 0; v < mesh._vertices.size(); ++v)
        {
            AdvancedVertex* vertex = mesh._vertices[v];
            vertex->setValue(smoothedVertices[v]);
        }
    }
}

void AdvancedMesh::getVerticesAndTrianglesArray(Vertex *& vertexArray,
                                                Triangle *& triangleArray,
                                                uint64_t& numberVertices,
                                                uint64_t& numberTriangles)
{
    // Generic
    Node *node;
    AdvancedVertex *vertex;

    numberVertices = I2UI64(_vertices.numberElements());
    numberTriangles = I2UI64(_triangles.numberElements());

    vertexArray = new Vertex[numberVertices];
    triangleArray = new Triangle[numberTriangles];

    // Write vertices
    uint64_t i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertexArray[i++] = Vertex(vertex->x, vertex->y, vertex->z);
    }

    // Copy the vertices
    float *auxVertices = new float[I2UI64(_vertices.numberElements())];

    // Construct the faces from the vertices
    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        auxVertices[i++] = D2F(vertex->x);
    }

    // Reset the counter
    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = i++;
    }

    i = 0;
    FOR_EACH_NODE(_triangles, node)
    {
        triangleArray[i][0] = VERTEX_1(node);
        triangleArray[i][1] = VERTEX_2(node);
        triangleArray[i][2] = VERTEX_3(node);
        i++;
    }

    // Reset the counter
    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = F2D(auxVertices[i++]);
    }

    // Free the auxiliary array
    delete[] auxVertices;
}

void AdvancedMesh::printMeshStats(const std::string &reference,
                                  const std::string *prefix)
{
    LOG_TITLE("Mesh Statistics");

    LOG_STATUS("Collecting Stats.");
    double area = this->area();
    double volume = this->volume();

    // Get the data
    Vertex *vertexArray;
    Triangle *triangleArray;
    uint64_t numberVertices;
    uint64_t numberTriangles;
    getVerticesAndTrianglesArray(vertexArray, triangleArray, numberVertices, numberTriangles);
    Ultraliser::MeshStatistics stats(vertexArray, triangleArray, numberVertices, numberTriangles);
    stats.writeStatsDistributions(*prefix + "-" + reference);

    AdvancedPoint pMin, pMax;
    getBoundingBox(pMin, pMax);
    AdvancedPoint bounds = pMax - pMin;

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + MESH_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                bounds.x, bounds.y, bounds.z);
        fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                pMin.x, pMin.y, pMin.z);
        fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                pMax.x, pMax.y, pMax.z);
        fprintf(info, "\t* Number Vertices       | %s \n",
                FORMAT(_vertices.numberElements()));
        fprintf(info, "\t* Number Triangles      | %s \n",
                FORMAT(_triangles.numberElements()));
        fprintf(info, "\t* Surface Area          | %f² \n",
                F2D(area));
        fprintf(info, "\t* Volume                | %f³ \n",
                F2D(volume));

        // Close the file
        fclose(info);
    }

    LOG_STATUS_IMPORTANT("Mesh Stats. [ %s ]", reference.c_str());
    LOG_INFO("\t* Bounding Box:         | [%f, %f, %f]",
             bounds.x, bounds.y, bounds.z);
    LOG_INFO("\t* pMin:                 | [%f, %f, %f]",
             pMin.x, pMin.y, pMin.z);
    LOG_INFO("\t* pMax:                 | [%f, %f, %f]",
            pMax.x, pMax.y, pMax.z);
    LOG_INFO("\t* Number Vertices       | %s",
             FORMAT(_vertices.numberElements()));
    LOG_INFO("\t* Number Triangles      | %s",
             FORMAT(_triangles.numberElements()));
    LOG_INFO("\t* Surface Area          | %f²",
             F2D(area));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volume));
}

}

