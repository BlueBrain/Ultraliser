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
#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/meshes/advanced/Macros.hh>
#include <data/meshes/advanced/Defines.hh>

namespace Ultraliser
{

double subSurfBetaLoop(int k)
{
    const double beta = (cos((2.0 * M_PI) / ((double) k)) / 4.0) + (3.0 / 8.0);
    return ((5.0 / 8.0) - (beta * beta)) / ((double) k);
}

void loopRelaxOriginal(AdvancedVertex *vertex)
{
    // Generic
    Node *node;
    AdvancedEdge *edge;
    List *incidentEdges;
    AdvancedPoint newPoint;

    // If the vertex is on a boundary edge
    if (vertex->isOnBoundary())
    {
        // List of incident edges
        incidentEdges = vertex->getIncidentEdges();

        // Create a new point
        newPoint = (*(((AdvancedEdge *) incidentEdges->head()->data)->oppositeVertex(vertex)));
        newPoint += (*(((AdvancedEdge *) incidentEdges->tail()->data)->oppositeVertex(vertex)));
        newPoint = (((*vertex) * 6.0) + newPoint) / 8.0;

        delete(incidentEdges);
    }
    else
    {
         // List of incident edges
        incidentEdges = vertex->getIncidentEdges();

        const size_t k = incidentEdges->numberElements();
        const double betaFactor = subSurfBetaLoop(k);

        FOR_EACH_VE_EDGE(incidentEdges, edge, node)
        {
            newPoint += (*(edge->oppositeVertex(vertex)));
        }

        newPoint = ((*vertex) * (1.0 - k * betaFactor)) + (newPoint * betaFactor);
        delete(incidentEdges);
    }

    vertex->info = new AdvancedPoint(&newPoint);
}

bool remintsAppendCubeToList(AdvancedTriangle *inputTriangulation, List& inputList)
{
    if (!IS_VISITED(inputTriangulation) || IS_BIT(inputTriangulation, 6))
        return false;

    AdvancedTriangle *triangle, *currentTriangle;
    AdvancedVertex *vertex;
    List triangleList(inputTriangulation);

    MARK_BIT(inputTriangulation,6);

    double minX = DBL_MAX, maxX = -DBL_MAX;
    double minY = DBL_MAX, maxY = -DBL_MAX;
    double minZ = DBL_MAX, maxZ = -DBL_MAX;

    while(triangleList.numberElements())
    {
        triangle = (AdvancedTriangle *)triangleList.popHead();

        vertex = triangle->v1();
        minX = MIN(minX, vertex->x);
        minY = MIN(minY, vertex->y);
        minZ = MIN(minZ, vertex->z);
        maxX = MAX(maxX, vertex->x);
        maxY = MAX(maxY, vertex->y);
        maxZ = MAX(maxZ, vertex->z);

        vertex = triangle->v2();
        minX = MIN(minX, vertex->x);
        minY = MIN(minY, vertex->y);
        minZ = MIN(minZ, vertex->z);
        maxX = MAX(maxX, vertex->x);
        maxY = MAX(maxY, vertex->y);
        maxZ = MAX(maxZ, vertex->z);

        vertex = triangle->v3();
        minX = MIN(minX, vertex->x);
        minY = MIN(minY, vertex->y);
        minZ = MIN(minZ, vertex->z);
        maxX = MAX(maxX, vertex->x);
        maxY = MAX(maxY, vertex->y);
        maxZ = MAX(maxZ, vertex->z);

        if ((currentTriangle = triangle->t1()) != nullptr &&
                !IS_BIT(currentTriangle, 6) && IS_VISITED(currentTriangle))
        {
            triangleList.appendHead(currentTriangle);
            MARK_BIT(currentTriangle, 6);
        }

        if ((currentTriangle = triangle->t2()) != nullptr &&
                !IS_BIT(currentTriangle, 6) && IS_VISITED(currentTriangle))
        {
            triangleList.appendHead(currentTriangle);
            MARK_BIT(currentTriangle, 6);
        }

        if ((currentTriangle = triangle->t3()) != nullptr &&
                !IS_BIT(currentTriangle, 6) && IS_VISITED(currentTriangle))
        {
            triangleList.appendHead(currentTriangle);
            MARK_BIT(currentTriangle, 6);
        }
    }

    inputList.appendTail(new AdvancedPoint(minX, minY, minZ));
    inputList.appendTail(new AdvancedPoint(maxX, maxY, maxZ));

    return true;
}

bool remintsIsVertexInCube(AdvancedVertex *vertex, List& inputList)
{
    Node *iNode;
    AdvancedPoint *p1, *p2;
    FOR_EACH_NODE(inputList, iNode)
    {
        p1 = (AdvancedPoint *)iNode->data;
        iNode = iNode->next();

        p2 = (AdvancedPoint *)iNode->data;
        if (!(vertex->x < p1->x || vertex->y < p1->y || vertex->z < p1->z ||
              vertex->x > p2->x || vertex->y > p2->y || vertex->z > p2->z))
            return true;
    }

    return false;
}

void remintsSelectTrianglesInCubes(AdvancedMesh *inputMesh)
{
    AdvancedTriangle *triangle;
    AdvancedVertex *vertex;
    Node *node;
    List loc;

    FOR_EACH_VT_TRIANGLE((&(inputMesh->_triangles)), triangle, node)
    {
        remintsAppendCubeToList(triangle, loc);
    }

    FOR_EACH_VV_VERTEX((&(inputMesh->_vertices)), vertex, node)
    {
        if (remintsIsVertexInCube(vertex, loc))
        {
            MARK_BIT(vertex, 5);
        }
    }

    FOR_EACH_VT_TRIANGLE((&(inputMesh->_triangles)), triangle, node)
    {
        UNMARK_BIT(triangle, 6);
        if (IS_BIT(triangle->v1(), 5) || IS_BIT(triangle->v2(), 5) || IS_BIT(triangle->v3(), 5))
        {
            MARK_VISIT(triangle);
        }
    }

    FOR_EACH_VV_VERTEX((&(inputMesh->_vertices)), vertex, node)
    {
        UNMARK_BIT(vertex, 5);
    }
    loc.freeNodes();
}

void jitterIncrease(char *charArray)
{
    bool isNegative = (charArray[0] == '-');
    int stringLength = strlen(charArray);

    if (isNegative)
    {
        for (int i = stringLength - 1; i >= 1; i--)
            if (charArray[i] == '0')
                charArray[i] = '9';
            else if (charArray[i] == '.')
                continue;
            else
            {
                charArray[i]--;
                break;
            }
    }
    else
    {
        for (int i = stringLength - 1; i >= 0; i--)
        {
            if (charArray[i] == '9')
                charArray[i] = '0';
            else if (charArray[i] == '.')
                continue;
            else
            {
                charArray[i]++;
                break;
            }
        }
    }
}

void jitterDecrease(char *charArray)
{
    bool isNegative = (charArray[0] == '-');
    int stringLength = strlen(charArray);

    if (isNegative)
    {
        for (int i = stringLength - 1; i >= 1; i--)
        {
            if (charArray[i] == '9')
                charArray[i] = '0';
            else if (charArray[i] == '.')
                continue;
            else
            {
                charArray[i]++;
                break;
            }
        }
    }
    else
    {
        for (int i = stringLength - 1; i >= 0; i--)
        {
            if (charArray[i] == '0')
                charArray[i] = '9';
            else if (charArray[i] == '.')
                continue;
            else
            {
                charArray[i]--;
                break;
            }
        }
    }
}

void jitterCoordinate(double& coordinate, int j)
{
    char floatVersion[32];
    float floatValue;

    sprintf(floatVersion, "%f", NODE_TO_FLOAT(coordinate));

    if (j > 0)
        jitterIncrease(floatVersion);
    else if (j < 0)
        jitterDecrease(floatVersion);

    sscanf(floatVersion, "%f", &floatValue);
    coordinate = floatValue;
}

}
