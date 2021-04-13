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
#include <data/meshes/advanced/primitives/LoopSplit.h>
#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/Utilities.h>

namespace Ultraliser
{

void AdvancedMesh::loopSubdivision(const bool &midPoint)
{
    // Generic
    List nodes;
    Node *node;
    AdvancedVertex *vertex;
    AdvancedEdge *edge, *edge1, *edge2;
    AdvancedTriangle *triangle;
    AdvancedLoopSplit *loopSplit;

    if (!midPoint)
    {
        FOR_EACH_VERTEX(vertex, node)
        {
            loopRelaxOriginal(vertex);
        }
    }

    bool isSelection = false;
    FOR_EACH_TRIANGLE(triangle, node)
    {
        if (IS_VISITED(triangle))
        {
            isSelection = true;
            break;
        }
    }

    if (isSelection)
    {
        FOR_EACH_TRIANGLE(triangle, node)
        {
            if (IS_VISITED(triangle))
            {
                MARK_BIT(triangle->edge1, 3);
                MARK_BIT(triangle->edge2, 3);
                MARK_BIT(triangle->edge3, 3);
            }
        }
    }
    else
    {
        FOR_EACH_EDGE(edge, node)
        {
            MARK_BIT(edge, 3);
        }
    }

    bool detectSharp = false;
    FOR_EACH_EDGE(edge, node)
    {
        if (IS_BIT(edge, 3))
        {
            MARK_BIT(edge->v1, 3);
            MARK_BIT(edge->v2, 3);

            if (!midPoint && IS_SHARPEDGE(edge))
                detectSharp = true;

            nodes.appendHead(new AdvancedLoopSplit(edge, midPoint));
        }
    }

    if (detectSharp)
    {
        LOG_WARNING("loopSubdivision: Crease-preservation is NOT supported!");
    }

    int isOnBoundary;
    FOR_EACH_NODE(nodes, node)
    {
        loopSplit = ((AdvancedLoopSplit *) node->data);

        isOnBoundary = loopSplit->edge->isOnBoundary();

        vertex = splitEdge(loopSplit->edge, &(loopSplit->point));

        edge1 = (AdvancedEdge *) _edges.head()->data;
        edge2 = (AdvancedEdge *) _edges.head()->next()->data;

        MARK_VISIT2(vertex);
        MARK_VISIT2(edge1);

        if (!isOnBoundary)
        {
            MARK_VISIT2(edge2);
        }

        if (loopSplit->edge->t2)
        {
            if (IS_VISITED(loopSplit->edge->t2))
            {
                triangle = (AdvancedTriangle *) _triangles.head()->data;
                MARK_VISIT(triangle);
            }

            if (loopSplit->edge->t1 && IS_VISITED(loopSplit->edge->t1))
            {
                triangle = (AdvancedTriangle *) _triangles.head()->next()->data;
                MARK_VISIT(triangle);
            }
        }
        else if (IS_VISITED(loopSplit->edge->t1))
        {
            triangle = (AdvancedTriangle *)_triangles.head()->data;
            MARK_VISIT(triangle);
        }

        if (IS_SHARPEDGE(loopSplit->edge))
        {
            if (isOnBoundary)
            {
                TAG_SHARPEDGE(edge2);
            }
            else
            {
                TAG_SHARPEDGE((AdvancedEdge *)_edges.head()->next()->next()->data);
            }
        }
    }

    nodes.freeNodes();

    FOR_EACH_EDGE(edge, node)
    {
        if (IS_VISITED2(edge))
        {
            UNMARK_VISIT2(edge);

            if ((IS_VISITED2(edge->v1) && !IS_VISITED2(edge->v2)) ||
               (!IS_VISITED2(edge->v1) && IS_VISITED2(edge->v2)))
            {
                if (edge->swap(1) && (!IS_VISITED2(edge->v1) || !IS_VISITED2(edge->v2)))
                {
                    edge->swap(1);
                }
            }
        }
    }

    FOR_EACH_VERTEX(vertex, node)
    {
        if (!IS_VISITED2(vertex) && !midPoint && IS_BIT(vertex, 3))
        {
            vertex->setValue((AdvancedPoint *) vertex->info);
            delete((AdvancedPoint *) vertex->info);
            vertex->info = nullptr;
        }
    }

    FOR_EACH_VERTEX(vertex, node)
    {
        UNMARK_VISIT2(vertex);
        UNMARK_BIT(vertex, 3);
    }

    if (detectSharp)
    {
        LOG_WARNING("loopSubdivision: Tagged sharp edges have been smoothed!");
    }

    FOR_EACH_EDGE(edge, node)
    {
        UNMARK_BIT(edge, 3);
    }
}

}
