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
    ~Graph();

    Graph(SkeletonWeightedEdges& weighteEdges, GraphNodes& graphNodes);

    void _addEdges(SkeletonWeightedEdges& edges);

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
