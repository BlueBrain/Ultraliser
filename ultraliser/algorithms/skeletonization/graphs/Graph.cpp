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

#include "Graph.h"

namespace Ultraliser
{

Graph::Graph(const size_t& numberNodes)
{
    // Get the total number of nodes in the graph
    this->_numberNodes = numberNodes;

    // Allocate the adjacency list
    _adjacencyLists = new std::list< size_t >[numberNodes];
}

Graph::Graph(SkeletonWeightedEdges& weighteEdges, GraphNodes& graphNodes)
{
    // Get the total number of nodes in the graph
    _numberNodes = graphNodes.size();

    // Allocate the adjacency list
    _adjacencyLists = new std::list< size_t >[_numberNodes];

    _addEdges(weighteEdges);
}

Graph::~Graph()
{
    delete[] _adjacencyLists;
}

void Graph::_addEdges(SkeletonWeightedEdges& edges)
{
    for (size_t i =0; i < edges.size(); i++)
    {
        auto n1 = edges[i]->node1;
        auto n2 = edges[i]->node2;

        _adjacencyLists[n1->graphIndex].push_back(n2->graphIndex);
        _adjacencyLists[n2->graphIndex].push_back(n1->graphIndex);
    }
}

void Graph::addEdge(const size_t& n1, const size_t& n2)
{
    _adjacencyLists[n1].push_back(n2);
    _adjacencyLists[n2].push_back(n1);
}

void Graph::_makeDSF(size_t nodeIndex, bool* visited, GraphComponent& component)
{
    // Set the current node to visited
    visited[nodeIndex] = true;

    // Add the node index to the component
    component.push_back(nodeIndex);

    // Recur for all the vertices
    // adjacent to this vertex

    // Apply the operation to all the nodes that are connected, i.e. adjacent to the current node
    std::list< size_t >::iterator i;
    for (i = _adjacencyLists[nodeIndex].begin(); i != _adjacencyLists[nodeIndex].end(); ++i)
    {
        // Ensure that the node is not visited
        if (!visited[*i])
        {
            // Recursively, make a DFS
            _makeDSF(*i, visited, component);
        }
    }
}

GraphComponents Graph::getComponents()
{
    // The components list
    GraphComponents components;

    // Mark all the nodes as not visited
    bool* visited = new bool[_numberNodes];

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _numberNodes; ++i)
        visited[i] = false;

    for (size_t i = 0; i < _numberNodes; ++i)
    {
        if (visited[i] == false)
        {
            // Establish a new component and get all the nodes within it
            GraphComponent component;

            // Make a DFS and get all the nodes connected to the component
            _makeDSF(i, visited, component);

            // Add the GraphComponent to GraphComponents list
            components.push_back(component);
        }
    }

    // Delete the temporary visited array
    delete[] visited;

    // Return the components list
    return components;
}

}
