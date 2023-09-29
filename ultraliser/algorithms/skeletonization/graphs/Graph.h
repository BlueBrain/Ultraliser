/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <common/Headers.hh>
#include <algorithms/skeletonization/graphs/GraphNode.h>

namespace Ultraliser
{

/**
 * @brief GraphComponent
 */
typedef std::vector< size_t > GraphComponent;

/**
 * @brief GraphComponents
 */
typedef std::vector< GraphComponent > GraphComponents;

/**
 * @brief The Graph class
 * Represents a undirected graph using adjacency list representation.
 */
class Graph
{
public:

    /**
     * @brief Graph
     * @param numberNodes
     */
    Graph(const size_t& numberNodes);

    /**
     * @brief Graph
     * @param weighteEdges
     * @param graphNodes
     */
    Graph(SkeletonWeightedEdges& weighteEdges, GraphNodes& graphNodes);
     ~Graph();

public:

    /**
     * @brief addEdge
     * Adds an undirected edge to the graph, based on the indices.
     * @param n1
     * @param n2
     */
    void addEdge(const size_t& n1, const size_t& n2);

    /**
     * @brief getComponents
     * @return
     */
    GraphComponents getComponents();

private:

    /**
     * @brief _addEdges
     * @param edges
     */
    void _addEdges(SkeletonWeightedEdges& edges);

private:

    /**
     * @brief _numberNodes
     */
    size_t _numberNodes;

    /**
     * @brief _adjacencyLists
     * Pointer to an array containing adjacency lists.
     */
    std::list< size_t >* _adjacencyLists;

    /**
     * @brief _makeDSF
     * @param nodeIndex
     * @param visited
     * @param component
     */
    void _makeDSF(size_t nodeIndex, bool* visited, GraphComponent& component);
};

}
