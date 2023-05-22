#include "Mesh.h"

#include <cassert>
#include <functional>
#include <unordered_map>

namespace
{
struct Edge
{
    size_t a, b;

    bool operator==(const Edge &other) const
    {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
    }
};
}

namespace std
{
template<>
struct hash<Edge>
{
    size_t operator()(const Edge& edge) const
    {
        auto mi = std::min(edge.a, edge.b);
        auto ma = std::max(edge.a, edge.b);
        size_t hash = 23;
        hash = hash * 31 + mi;
        hash = hash * 31 + ma;
        return hash;
    }
};
}

namespace
{
struct HalfEdge
{
    size_t startVertex;
    size_t twinEdge;
};

class HalfEdgeMesh
{
public:
    HalfEdgeMesh(Ultraliser::Mesh &mesh, std::vector<HalfEdge> halfEdges)
     : _mesh(mesh)
     , _halfEdges(std::move(halfEdges))
    {
    }

    HalfEdge &getHalfEdge(size_t triangleIndex, size_t halfEdge)
    {
        assert(halfEdge < 3);
        assert(triangleIndex < _mesh.getNumberTriangles());

        auto halfEdgeIndex = 3 * triangleIndex + halfEdge;
        assert(halfEdgeIndex < _halfEdges.size());

        return _halfEdges[halfEdgeIndex];
    }

    HalfEdge &getPrevHalfEdge(size_t halfEdgeIndex)
    {
        assert(halfEdgeIndex < _halfEdges.size());

        auto prevHalfEdgeIndex = halfEdgeIndex % 3 == 0? halfEdgeIndex + 2 : halfEdgeIndex - 1;
        assert(prevHalfEdgeIndex < _halfEdges.size());

        return _halfEdges[prevHalfEdgeIndex];
    }

    HalfEdge &getNextHalfEdge(size_t halfEdgeIndex)
    {
        assert(halfEdgeIndex < _halfEdges.size());

        auto nextHalfEdgeIndex = halfEdgeIndex % 2 == 0? halfEdgeIndex - 2 : halfEdgeIndex + 1;
        assert(nextHalfEdgeIndex < _halfEdges.size());

        return _halfEdges[nextHalfEdgeIndex];
    }

    HalfEdge &rotateHalfEdgeCounterClockWise(size_t halfEdgeIndex)
    {
        auto &prev = getPrevHalfEdge(halfEdgeIndex);
        return _halfEdges[prev.twinEdge];
    }

    HalfEdge &rotateHalfEdgeClockWise(size_t halfEdgeIndex)
    {
        assert(halfEdgeIndex < _halfEdges.size());
        return getNextHalfEdge(_halfEdges[halfEdgeIndex].twinEdge);
    }

    const Ultraliser::Triangle &getTriangle(size_t halfEdgeIndex) const
    {
        assert(halfEdgeIndex < _halfEdges.size());

        auto faceIndex = halfEdgeIndex / 3;
        assert(faceIndex < _mesh.getNumberTriangles());

        return _mesh.getTriangles()[faceIndex];
    }

private:
    Ultraliser::Mesh &_mesh;
    std::vector<HalfEdge> _halfEdges;
};

class HalfEdgeMeshFactory
{
public:
    static HalfEdgeMesh create(Ultraliser::Mesh &mesh)
    {
        auto numTriangles = mesh.getNumberTriangles();
        auto triangles = mesh.getTriangles();

        auto halfEdgeCache = std::unordered_map<Edge, size_t>();

        auto halfEdges = std::vector<HalfEdge>();

        for(size_t i = 0; i < numTriangles; ++i)
        {
            auto &triangle = triangles[i];
            for(size_t v = 0; v < 3; ++v)
            {
                auto vertex = triangle[v];
                auto nextVertex = triangle[(v + 1) % 3];

                auto halfEdgeIndex = halfEdges.size();
                auto &halfEdge = halfEdges.emplace_back();
                halfEdge.startVertex = vertex;

                auto edge = Edge{vertex, nextVertex};
                auto it = halfEdgeCache.find(edge);
                if(it == halfEdgeCache.end())
                {
                    halfEdgeCache[edge] = halfEdgeIndex;
                    continue;
                }

                auto twinIndex = it->second;
                halfEdge.twinEdge = twinIndex;
                halfEdges[twinIndex].twinEdge = halfEdgeIndex;
            }
        }

        return HalfEdgeMesh(mesh, std::move(halfEdges));
    }
};
}

namespace Ultraliser
{
void Mesh::collapseEdges(float percentageToCollapse)
{
    auto halfEdgeMesh = HalfEdgeMeshFactory::create(*this);
}
}
