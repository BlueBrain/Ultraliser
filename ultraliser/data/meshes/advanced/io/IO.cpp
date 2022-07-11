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

#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/meshes/advanced/Defines.hh>
#include <utilities/Utilities.h>

namespace Ultraliser
{

void AdvancedMesh::importMesh(const std::string &filePath, const bool update)
{
    // Start the timer
    TIMER_SET;

    LOG_TITLE("Importing Mesh");

    std::ifstream fileStream;
    fileStream.open(filePath.c_str());
    if (!fileStream.good())
    {
        LOG_ERROR("Cannot load mesh file [ %s ]", filePath.c_str());
    }
    fileStream.close();

    // Switch to the corresponding loader based on the extension
    switch(filePath[strlen(filePath.c_str()) - 1])
    {
    case 'y':
    case 'Y':
    {
        importPLY(filePath);
        break;
    }
    case 'j':
    case 'J':
    {
        importOBJ(filePath);
        break;
    }
    case 'l':
    case 'L':
    {
        importSTL(filePath);
        break;
    }
    case 'f':
    case 'F':
    {
        importOFF(filePath);
        break;
    }
    default:
        LOG_ERROR("The file format of this mesh file [ %s ] is NOT supported.", filePath.c_str());
    }

    // If no error is detected and the update is requested, update the mesh
    if (update)
    {
        eulerUpdate();
    }

    // Statistics
    LOG_STATUS_IMPORTANT("Importing Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void AdvancedMesh::append(const std::string &filePath, const bool update)
{
    // If no meshes exist, then just load the input mesh.
    if (!_triangles.numberElements())
        return importMesh(filePath, update);

    // Create a new mesh object and import the given mesh into it
    AdvancedMesh newMesh;
    newMesh.importMesh(filePath, 0);

    // Append the vertices
    _vertices.joinTailList(&(newMesh._vertices));

    // Append the edges
    _edges.joinTailList(&(newMesh._edges));

    // Append the triangles
    _triangles.joinTailList(&(newMesh._triangles));
    if (update)
        eulerUpdate();
    else
        _dBoundaries = _dHandles = _dShells = 1;
}

void AdvancedMesh::exportMesh(const std::string &prefix,
                              const bool& formatOBJ,
                              const bool& formatPLY,
                              const bool& formatOFF,
                              const bool& formatSTL)
{
    // Start the timer
    TIMER_SET;

    LOG_TITLE("Exporting Mesh");


    if (formatOBJ)
    {
        std::stringstream fileName;
        fileName << prefix << OBJ_EXTENSION;
        exportOBJ(fileName.str());
    }

    if (formatOFF)
    {
        std::stringstream fileName;
        fileName << prefix << OFF_EXTENSION;
        exportOFF(fileName.str());
    }

    if (formatPLY)
    {
        std::stringstream fileName;
        fileName << prefix << PLY_EXTENSION;
        exportPLY(fileName.str());
    }

    if (formatSTL)
    {
        std::stringstream fileName;
        fileName << prefix << STL_EXTENSION;
        exportSTL(fileName.str());
    }

    // Statistics
    LOG_STATUS_IMPORTANT("Exporting Mesh Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

AdvancedTriangle* AdvancedMesh::createTriangleFromVertices(ExtendedVertex *vertex1,
                                                           ExtendedVertex *vertex2,
                                                           ExtendedVertex *vertex3)
{
    // Generic
    AdvancedEdge *edge1, *edge2, *edge3;
    AdvancedTriangle *triangle = nullptr;

    // Create the first edge
    edge1 = createEdge(vertex1, vertex2);
    if (edge1->t1 != nullptr && edge1->t2 != nullptr)
    {
        MARK_BIT(edge1, 5);
    }

    // Create the second edge
    edge2 = createEdge(vertex2, vertex3);
    if (edge2->t1 != nullptr && edge2->t2 != nullptr)
    {
        MARK_BIT(edge2, 5);
    }

    // Create the third edge
    edge3 = createEdge(vertex3, vertex1);
    if (edge3->t1 != nullptr && edge3->t2 != nullptr)
    {
        MARK_BIT(edge3, 5);
    }

    // Verify edge 1
    if (IS_BIT(edge1, 5))
    {
        edge1 = createEdge(vertex1, vertex2, 0);
        MARK_BIT(edge1, 5);
    }

    // Verify edge 2
    if (IS_BIT(edge2, 5))
    {
        edge2 = createEdge(vertex2, vertex3, 0);
        MARK_BIT(edge2, 5);
    }

    // Verify edge 3
    if (IS_BIT(edge3, 5))
    {
        edge3 = createEdge(vertex3, vertex1, 0);
        MARK_BIT(edge3, 5);
    }

    // Create an unoriented triangle from the three edges
    triangle = createUnorientedTriangle(edge1, edge2, edge3);

    // If the triangle is not created (not valid triangle), then remove the edges
    if (triangle == nullptr)
    {
        LOG_WARNING("EMPTY TRIANGLE");
        if (edge3->t1 == nullptr && edge3->t2 == nullptr)
        {
            _edges.freeNode(edge3);
            vertex3->VE.removeNode(edge3);
            vertex1->VE.removeNode(edge3);

            if (vertex3->v->e0 == edge3)
                vertex3->v->e0 = nullptr;
            if (vertex1->v->e0 == edge3)
                vertex1->v->e0 = nullptr;
        }
        if (edge2->t1 == nullptr && edge2->t2 == nullptr)
        {
            _edges.freeNode(edge2);
            vertex2->VE.removeNode(edge2);
            vertex3->VE.removeNode(edge2);

            if (vertex2->v->e0 == edge2)
                vertex2->v->e0 = nullptr;
            if (vertex3->v->e0 == edge2)
                vertex3->v->e0 = nullptr;
        }
        if (edge1->t1 == nullptr && edge1->t2 == nullptr)
        {
            _edges.freeNode(edge1);
            vertex1->VE.removeNode(edge1);
            vertex2->VE.removeNode(edge1);

            if (vertex1->v->e0 == edge1)
                vertex1->v->e0 = nullptr;
            if (vertex2->v->e0 == edge1)
                vertex2->v->e0 = nullptr;
        }
    }

    return triangle;
}

AdvancedTriangle* AdvancedMesh::createIndexedTriangle(ExtendedVertex **vertexList,
                                                      int64_t i1, int64_t i2, int64_t i3)
{
    return createTriangleFromVertices(vertexList[i1], vertexList[i2], vertexList[i3]);
}

void AdvancedMesh::_closeLoadingSession(FILE *file, int loadedFaces,
                                        ExtendedVertex **var, bool triangulate)
{
    // Close the file
    fclose(file);

    // Delete the extended vertices
    if (var != nullptr)
    {
        for (uint64_t i = 0; i < _vertices.numberElements(); ++i)
            delete(var[i]);
        free(var);
    }

    if (loadedFaces)
    {
        if (triangulate)
        {
            LOG_WARNING("Some polygonal faces needed to be triangulated.");
        }

        // Fix the connectivity of the mesh
        fixConnectivity();

        // Update the data
        eulerUpdate();

        // Update the data
        _dBoundaries = _dHandles = _dShells = 1;
    }
}

void AdvancedMesh::coordBackApproximation()
{
    Node *node;
    AdvancedVertex *vertex;
    char floatVertex[32];
    float value;

    FOR_EACH_VERTEX(vertex, node)
    {
        sprintf(floatVertex, "%f", NODE_TO_DOUBLE(vertex->x));
        sscanf(floatVertex, "%f", &value); vertex->x = F2D(value);
        sprintf(floatVertex, "%f", NODE_TO_DOUBLE(vertex->y));
        sscanf(floatVertex, "%f", &value); vertex->y = F2D(value);
        sprintf(floatVertex, "%f", NODE_TO_DOUBLE(vertex->z));
        sscanf(floatVertex, "%f", &value); vertex->z = F2D(value);
    }
}

/**
 * @brief Mesh::_pinch
 * Implements the cutting and stitching procedure to convert to manifold mesh.
 * Assumes that singular edges to be cut and stitched are marked as BIT5.
 *
 * @param inputEdge
 * @param withCommonVertex
 * @return
 */
bool AdvancedMesh::_pinch(AdvancedEdge *inputEdge, bool withCommonVertex)
{
    List *boundaryEdges = (List *) inputEdge->info;
    if (boundaryEdges == nullptr)
        return false;

    Node *node = nullptr;
    AdvancedEdge *edge = nullptr;
    List *vertexEdges;

    if (withCommonVertex)
    {
        inputEdge->v1->e0 = inputEdge;
        vertexEdges = inputEdge->v1->getIncidentEdges();

        FOR_EACH_VE_EDGE(vertexEdges, edge, node)
        {
            if (edge != inputEdge && edge->isOnBoundary() &&
             (*(edge->oppositeVertex(inputEdge->v1))) == (*(inputEdge->v2)) &&
                inputEdge->merge(edge))
                break;
        }

        // Delete the vertex
        delete vertexEdges;

        // Deleted node
        if (node == nullptr)
        {
            inputEdge->v2->e0 = inputEdge;
            vertexEdges = inputEdge->v2->getIncidentEdges();

            FOR_EACH_VE_EDGE(vertexEdges, edge, node)
            {
                if (edge != inputEdge && edge->isOnBoundary() &&
                 (*(edge->oppositeVertex(inputEdge->v2))) == (*(inputEdge->v1)) &&
                    inputEdge->merge(edge))
                    break;
            }

            // Delete the vertex
            delete vertexEdges;
        }
    }
    else
    {
        if (inputEdge->t1 != nullptr)
        {
            FOR_EACH_VE_EDGE(boundaryEdges, edge, node)
            {
                if (edge != inputEdge &&
               (((*(edge->v1)) == (*(inputEdge->v1)) && edge->t2 != nullptr) ||
                ((*(edge->v1)) == (*(inputEdge->v2)) && edge->t1 != nullptr)) &&
                    inputEdge->merge(edge))
                    break;
            }
        }
        else
        {
            FOR_EACH_VE_EDGE(boundaryEdges, edge, node)
            {
                if (edge != inputEdge &&
               (((*(edge->v1)) == (*(inputEdge->v1)) && edge->t1 != nullptr) ||
                ((*(edge->v1)) == (*(inputEdge->v2)) && edge->t2 != nullptr)) &&
                    inputEdge->merge(edge))
                    break;
            }
        }
    }

    if (node == nullptr)
        return false;

    boundaryEdges->removeNode(inputEdge);
    boundaryEdges->removeNode(edge);
    inputEdge->info = edge->info = nullptr;

    if (boundaryEdges->numberElements() == 0)
        delete boundaryEdges;

    AdvancedEdge *gEdge, *edge1 = nullptr, *edge2 = nullptr;

    vertexEdges = inputEdge->v1->getIncidentEdges();
    for (node = vertexEdges->head(); node != nullptr; node = node->next())
    {
        if ((gEdge = (AdvancedEdge *)node->data)->info != nullptr)
        {
            edge1 = gEdge;
            break;
        }
    }

    for (node = vertexEdges->tail(); node != nullptr; node = node->prev())
    {
        if ((gEdge = (AdvancedEdge *)node->data)->info != nullptr)
        {
            if ((*(gEdge->oppositeVertex(inputEdge->v1))) !=
                (*(edge1->oppositeVertex(inputEdge->v1))))
                edge1 = nullptr;
            break;
        }
    }
    delete vertexEdges;

    vertexEdges = inputEdge->v2->getIncidentEdges();
    for (node = vertexEdges->head(); node != nullptr; node = node->next())
    {
        if ((gEdge = (AdvancedEdge *)node->data)->info != nullptr)
        {
            edge2 = gEdge;
            break;
        }
    }

    for (node = vertexEdges->tail(); node != nullptr; node = node->prev())
    {
        if ((gEdge = (AdvancedEdge *)node->data)->info != nullptr)
        {
            if ((*(gEdge->oppositeVertex(inputEdge->v2))) !=
                (*(edge2->oppositeVertex(inputEdge->v2))))
                edge2 = nullptr;
            break;
        }
    }
    delete vertexEdges;

    if (edge1 != nullptr)
        _pinch(edge1, true);
    if (edge2 != nullptr)
        _pinch(edge2, true);

    // Done
    return true;
}

AdvancedEdge *AdvancedMesh::duplicateEdge(AdvancedEdge *e1)
{
    // Make sure that the edge exists
    if (e1->t1 == nullptr || e1->t2 == nullptr)
        return nullptr;

    // Create a new edge
    AdvancedEdge *e2 = newEdge(e1);

    // Append it
    _edges.appendHead(e2);

    e1->t2->replaceEdge(e1, e2);
    e2->t2 = e1->t2;
    e1->t2 = nullptr;

    // Return the created edge
    return e2;
}

int AdvancedMesh::cutAndStitch()
{
    // Generic
    AdvancedEdge *edge, *auxEdge;
    Node *node;

    // A list that will collect the edges that will be removed
    List singularEdges;

    // Marking a duplicate (or an auxiliary) edge
    FOR_EACH_EDGE(edge, node)
    {
        if (IS_BIT(edge, 5) && ((auxEdge = duplicateEdge(edge)) != nullptr))
            MARK_BIT(auxEdge, 5);
    }

    // Adding the singular edge to the list of edges to be removed
    int numberSingularEdges = 0;
    FOR_EACH_EDGE(edge, node)
    {
        if (IS_BIT(edge, 5))
        {
            numberSingularEdges++;
            singularEdges.appendHead(edge);
            UNMARK_BIT(edge, 5);
        }
    }

    if (numberSingularEdges > 0)
        LOG_WARNING("The mesh has [%d] singular edges", numberSingularEdges);

    // Check the orientation of the faces
    forceNormalConsistence();

    // Duplicate the non manifold vertices
    _loadedNonManifoldVertices= duplicateNonManifoldVertices();

    // Sort the singular edges
    singularEdges.sort(&lexEdgeCompare);

    // Free the edge
    FOR_EACH_EDGE(edge, node)
    {
        edge->info = nullptr;
    }

    // Delete the singular edges
    auxEdge = nullptr;
    FOR_EACH_VE_EDGE((&singularEdges), edge, node)
    {
        if (auxEdge == nullptr || lexEdgeCompare(edge, auxEdge) != 0)
        {
            edge->info = new List();
            auxEdge = edge;
        }

        ((List *)auxEdge->info)->appendTail(edge);
        edge->info = auxEdge->info;
    }

    // Now each edge is either 'regular' or has the info field pointing to a
    // list of coincident boundary edges
    // First, pinch bounded chains of singular edges starting from one endpoint
    FOR_EACH_VE_EDGE((&singularEdges), edge, node)
    {
        if (edge->isLinked())
            _pinch(edge, true);
    }

    // Then, pinch the remaining unbounded chains starting from any of the edges
    FOR_EACH_VE_EDGE((&singularEdges), edge, node)
    {
        if (edge->isLinked())
            _pinch(edge, false);
    }

    // Delete the edges
    removeUnlinkedElements();

    // Update the state to clean
    _dBoundaries = _dHandles = _dShells = 1;

    return singularEdges.numberElements();
}

}
