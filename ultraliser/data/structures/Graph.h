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

#ifndef ULTRALISER_DATA_STRUCTURES_GRAPH_H
#define ULTRALISER_DATA_STRUCTURES_GRAPH_H

#include <common/Common.h>
#include <data/structures/List.h>

namespace Ultraliser
{

class GraphEdge;

/**
 * @brief The GraphNode class
 * Base class type for nodes of non-oriented graphs.
 */
class GraphNode
{
public:

    /**
     * @brief edges
     * List of incident edges.
     */
    List edges;

    /**
     * @brief mask
     * Generic 8-bit mask for marking purposes.
     */
    unsigned char mask;

    /**
     * @brief GraphNode
     */
    GraphNode()
    {
        mask=0;
    }

    /**
     * @brief ~GraphNode
     */
    ~GraphNode()
    {
        /// EMPTY
    }

    /**
     * @brief isIsolated
     * Returns TRUE if the node is isolated. O(1).
     * @return
     */
    bool isIsolated()
    {
        return (edges.numberElements() == 0);
    }

    /**
     * @brief getEdge
     * Returns the edge connecting this with 'n'.
     * Returns nullptr if not connected. O(degree).
     * @param n
     * @return
     */
    GraphEdge *getEdge(GraphNode *n);
};


/**
 * @brief The GraphEdge class
 * Base class type for edges of non-oriented graphs.
 */
class GraphEdge
{
public:

    /**
     * @brief n1
     * First edge end point.
     */
    GraphNode *n1;

    /**
     * @brief n2
     * Second edge end point.
     */
    GraphNode *n2;

    /**
     * @brief info
     * Generic info attached to the edge
     */
    void *info;

    /**
     * @brief mask
     * Generic 8-bit mask for marking purposes.
     */
    unsigned char mask;

    /**
     * @brief GraphEdge
     */
    GraphEdge(GraphNode *, GraphNode *);
    ~GraphEdge()
    {
        /// EMPTY
    }

    /**
     * @brief oppositeNode
     * Returns the node oppsite to 'n'. O(1).
     * @param n
     * @return
     */
    GraphNode *oppositeNode(GraphNode *n)
    {
        return (n1 == n) ? (n2) : ((n2 == n) ? (n1) : nullptr);
    }

    /**
     * @brief isUnlinked
     * Returns TRUE if this edge does not connect points. O(1).
     * @return
     */
    bool isUnlinked()
    {
        return (n1 == nullptr);
    }

    /**
     * @brief hasNode
     * Returns TRUE if 'n' is a node of this edge. O(1).
     * @param n
     * @return
     */
    bool hasNode(GraphNode *n)
    {
        return (n1 == n || n2 == n);
    }

    /**
     * @brief makeUnlinked
     * Makes this edge as 'unlinked' from the graph. O(1).
     */
    void makeUnlinked()
    {
        n1=nullptr;
        n2=nullptr;
    }

    /**
     * @brief collapse
     * Edge collapse. O(degree of neighbors).
     * After one (or a series of) collapse, remember to call
     * Graph::deleteUnlinkedElements() to make the graph coherent with its
     * node and edge lists.
     */
    void collapse();

    /**
     * @brief invert
     * Inverts the edge's orientation
     */
    void invert();
};

/**
 * @brief The Graph class
 * Base class type for non oriented graphs
 */
class Graph
{
public:

    /**
     * @brief nodes
     * Nodes of the graph.
     */
    List nodes;

    /**
     * @brief edges
     * Edges of the graph.
     */
    List edges;

    ~Graph();

    /**
     * @brief addNode
     * Adds an existing isolated node to the graph. O(1).
     * @param n
     * @return
     */
    GraphNode *addNode(GraphNode *n)
    {
        nodes.appendHead(n);
        return n;
    }

    /**
     * @brief createEdge
     * Creates a new edge out of a pair of nodes. O(degree of nodes).
     * If the edges already exists, no new edge is created and the existing
     * one is returned. Otherwise the newly created edge is returned.
     * @param n1
     * @param n2
     * @return
     */
    GraphEdge *createEdge(GraphNode *n1, GraphNode *n2);

    /**
     * @brief unlinkEdge
     * Unlinks the edge and updates the graph connectivity accordingly.
     * @param e
     */
    void unlinkEdge(GraphEdge *e);

    /**
     * @brief destroyEdge
     * Removes and deletes the edge. Updates the graph connectivity accordingly.
     * @param e
     */
    void destroyEdge(GraphEdge *e);

    /**
     * @brief unlinkNode
     * Removes a node and updates the graph connectivity accordingly.
     * The unlinked node is returned.
     * @param n
     * @return
     */
    GraphNode *unlinkNode(GraphNode *n);

    /**
     * @brief deleteUnlinkedElements
     * Eliminates isolated nodes and unlinked edges from the lists. O(N).
     * The eliminated elements are deleted too.
     */
    void deleteUnlinkedElements();

    /**
     * @brief isConnected
     * @return
     */
    bool isConnected();
};

}

#endif // ULTRALISER_DATA_STRUCTURES_GRAPH_H
