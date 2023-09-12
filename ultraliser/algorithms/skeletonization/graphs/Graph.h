#pragma once

#include <common/Headers.hh>


namespace Ultraliser
{

// Graph class represents a undirected graph
// using adjacency list representation
class Graph
{
    int V; // No. of vertices

    // Pointer to an array containing adjacency lists
    std::list<int>* adj;

    // A function used by DFS
    void DFSUtil(int v, bool visited[]);

public:
    Graph(int V); // Constructor
    ~Graph();

    void addEdge(int v, int w);
    void connectedComponents();
};
}
