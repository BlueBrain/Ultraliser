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

#include <algorithms/skeletonization/SkeletonWeightedEdge.hh>

namespace Ultraliser
{

/**
 * @brief PathIndices
 */
typedef std::vector< int64_t > PathIndices;

/**
 * @brief PathsIndices
 */
typedef std::vector< PathIndices > PathsIndices;

class ShortestPathFinder
{
public:
    ShortestPathFinder(const SkeletonWeightedEdges &edges, const size_t& numberNodes);
    ~ShortestPathFinder();

    /**
     * @brief findPath
     * @param sourceNodeIndex
     * @param destinationNodeIndex
     * @return
     */
    PathIndices findPath(const size_t& sourceNodeIndex, const size_t& destinationNodeIndex);

private:

    /**
     * @brief _findShortestPathesWithDijkstra
     * @param sourceNodeIndex
     * @return
     */
    PathsIndices _findShortestPathesWithDijkstra(int64_t sourceNodeIndex);

    /**
     * @brief _constructPath
     * @param parent
     * @param i
     * @param path
     */
    void _constructPath(int64_t* parent, int64_t i, PathIndices& path);

    /**
     * @brief _constructGraphAdjacencyMatrix
     * Construct the adjacency matrix representation of the graph from the weighted edges list.
     * @param weightedEdges
     * A list of a given edges with weighted.
     * @param numberNodes
     * The total number of nodes in the graph.
     */
    void _constructGraphAdjacencyMatrix(const SkeletonWeightedEdges& edges);

    /**
     * @brief _computeMinimumDistanceIndex
     * @param distances
     * @param visited
     * @return
     */
    size_t _computeMinimumDistanceIndex(int64_t* distances, const bool* visited);

    void _allocateGraph();

    void _releaseGraph();

private:

    /**
     * @brief _graph
     * The graph adjacency matrix representation that is used to perform the search operation.
     */
    int64_t** _graph = nullptr;

    /**
     * @brief _numberNodes
     * Number of nodes in the graph.
     */
    const size_t _numberNodes;
};

}
