#pragma once

#include <common/Headers.hh>
#include <math/Vector3f.h>

namespace Ultraliser
{

/**
 * @brief The SkeletonNode struct
 * A node in the segmented skeleton.
 */
struct SkeletonNode
{
public:

    /**
     * @brief SkeletonNode
     */
    SkeletonNode()
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief SkeletonNode
     * @param index
     * @param point
     * @param voxel
     */
    SkeletonNode(const size_t& index, const Vector3f& point, const Vector3f& voxel)
    {
        this->index = index;
        this->point = point;
        this->voxel = voxel;
    }

public:

    /**
     * @brief index
     * The index of the node.
     */
    size_t index = 0;

    /**
     * @brief point
     * The Cartesian coordinates of the node.
     */
    Vector3f point;

    /**
     * @brief voxel
     * The corresponding voxel [X, Y, Z] in the volume used to extract the skeleton.
     */
    Vector3f voxel;

    /**
     * @brief radius
     * The radius of the node;
     */
    float radius;

    /**
     * @brief visited
     * If this flag is set, the node must have been visited before.
     * The flag is used to construct the connectivity of the skeleton.
     */
    bool visited = false;

    /**
     * @brief branching
     * If this flag is set, this node is a branching node, i.e. connected to more than two nodes.
     */
    bool branching = false;

    /**
     * @brief terminal
     * If this flag is set, this node is a terminal one, i.e. connected only to a single node.
     */
    bool terminal = false;

    /**
     * @brief insideSoma
     * If this flag is set, the node is located inside the soma and can be safely removed from
     * the graph.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool insideSoma = false;

    /**
     * @brief isSoma
     * If this flag is set, this node is a soma node.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool isSoma = false;

    /**
     * @brief connectedToSoma
     * If this flag is set, the node is connected to the soma, i.e. an initial node of a branch.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool connectedToSoma = false;

    /**
     * @brief edgeNodes
     * A list of all the nodes connected along the edges to this node.
     */
    std::vector< SkeletonNode* > edgeNodes;
};

/**
 * @brief SkeletonNodes
 */
typedef std::vector< SkeletonNode* > SkeletonNodes;
}
