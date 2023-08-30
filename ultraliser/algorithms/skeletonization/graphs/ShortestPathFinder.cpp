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

#include "ShortestPathFinder.h"

namespace Ultraliser
{

ShortestPathFinder::ShortestPathFinder(const SkeletonWeightedEdges &edges,
                                       const size_t& numberNodes)
    : _numberNodes(numberNodes)
{
   // std::cout << "Allocate \n";

    // Allocate the graph
    _allocateGraph();

   // std::cout << "Construct \n";

    // Construct the adjacency matrix
    _constructGraphAdjacencyMatrix(edges);

   // std::cout << "Construct DONE \n";

}

PathIndices ShortestPathFinder::findPath(const size_t& sourceNodeIndex,
                                         const size_t& destinationNodeIndex)
{
    std::cout << "findPath 1\n";
    // Find all the paths from the source node to all the other nodes in the graph
    PathsIndices paths = _findShortestPathesWithDijkstra(sourceNodeIndex);

    std::cout << "findPath 2\n";

    // Return the path that corresponds to the index of the destination node
    return paths[destinationNodeIndex];
}

PathsIndices ShortestPathFinder::_findShortestPathesWithDijkstra( int64_t sourceNodeIndex)
{
    std::cout << "xx1 \n";

    int64_t* distances = new int64_t[_numberNodes];
    bool* visited = new bool[_numberNodes];
    int64_t* parent = new int64_t[_numberNodes];

    std::cout << "xx2 \n";

    for (int i = 0; i < _numberNodes; i++)
    {
        parent[sourceNodeIndex] = -1;
        distances[i] = INT_MAX;
        visited[i] = false;
    }

    std::cout << "xx3 \n";


    distances[sourceNodeIndex] = 0;
    for (size_t i = 0; i < _numberNodes - 1; i++)
    {
        size_t minDistanceIndex = _computeMinimumDistanceIndex(distances, visited);

        // Visited
        visited[minDistanceIndex] = true;

        for (size_t j = 0; j < _numberNodes; j++)
        {
            int64_t currentDistance = distances[minDistanceIndex] + _graph[minDistanceIndex][j];

            if (!visited[j] && _graph[minDistanceIndex][j] && currentDistance < distances[j])
            {
                parent[j] = minDistanceIndex;
                distances[j] = currentDistance;
            }
        }
    }

    std::cout << "xx4 \n";


    PathsIndices paths;
    paths.resize(_numberNodes);

    for (int i = 0; i < _numberNodes; i++)
    {
        paths[i].push_back(sourceNodeIndex);
        _constructPath(parent, i, paths[i]);
    }

    std::cout << "xx5 \n";

    delete [] distances;
    delete [] visited;
    delete [] parent;

    // Returns a list of all the paths
    return paths;
}

void ShortestPathFinder::_constructPath(int64_t* parent, int64_t i, PathIndices& path)
{
    std::cout << i << ",";
    // We have reached the top, return
    if (parent[i] == -1) { return; }

    // Recursively, use the parent array to figure out the path in a reversed order
    _constructPath(parent, parent[i], path);

    // Add the index of the current node from the parent array
    path.push_back(i);
}


void ShortestPathFinder::_constructGraphAdjacencyMatrix(const SkeletonWeightedEdges &edges)
{
   // std::cout << "Num Nodes " << _numberNodes << "\n";

    for (size_t i = 0; i < _numberNodes; ++i)
    {
        for (size_t j = 0; j < _numberNodes; ++j)
        {
            std::cout << _graph[i][j] << " ";
        }
        std::cout << "\n";
    }

    {
        //std::cout << "edges.size() " << i << "\n";



//        _graph[edges[i]->node1WeightedIndex][edges[i]->node2WeightedIndex] = edges[i]->edgeWeight;
  //      _graph[edges[i]->node2WeightedIndex][edges[i]->node1WeightedIndex] = edges[i]->edgeWeight;
    }

    // Construct the graph
    for (size_t i = 0; i < edges.size(); ++i)
    {
        std::cout << "edges.size() " << i << " ";


        std::cout << edges[i]->node1->orderIndex << " " << edges[i]->node2->orderIndex << " " << edges[i]->edgeWeight << "\n";

        _graph[edges[i]->node1->orderIndex - 1][edges[i]->node2->orderIndex - 1] = edges[i]->edgeWeight;
        _graph[edges[i]->node2->orderIndex - 1][edges[i]->node1->orderIndex - 1] = edges[i]->edgeWeight;
    }

    for (size_t i = 0; i < _numberNodes; ++i)
    {
        for (size_t j = 0; j < _numberNodes; ++j)
        {
            std::cout << _graph[i][j] << " ";
        }
        std::cout << "\n";
    }

    std::cout << "DOnme \n";

}

size_t ShortestPathFinder::_computeMinimumDistanceIndex(int64_t* distances, const bool* visited)
{
    int64_t minDistance = INT64_MAX;
    size_t minDistanceIndex = UINT64_MAX;

    for (size_t i = 0; i < _numberNodes; ++i)
    {
        if (!visited[i] && distances[i] <= minDistance)
        {
            minDistance = distances[i];
            minDistanceIndex = i;
        }
    }
    return minDistanceIndex;
}

void ShortestPathFinder::_allocateGraph()
{
    _graph = new int64_t*[_numberNodes];
    for (size_t i = 0; i < _numberNodes; ++i)
    {
        _graph[i] = new int64_t[_numberNodes];
        for (size_t j = 0; j < _numberNodes; ++j)
        {
            _graph[i][j] = 0;
        }
    }
}

void ShortestPathFinder::_releaseGraph()
{
    for (size_t i = 0; i < _numberNodes; ++i)
    {
        delete [] _graph[i];
    }
    delete [] _graph;
}


ShortestPathFinder::~ShortestPathFinder()
{
    _releaseGraph();
}
}
