#pragma once

#include <common/Headers.hh>

namespace Ultraliser
{

/**
 * @brief The SkeletonEdge struct
 * An edge composed of two points (nodes or SkeletonNodes in the segmented skeleton.
 */
struct SkeletonEdge
{
public:

    /**
     * @brief SkeletonEdge
     */
    SkeletonEdge()
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief SkeletonEdge
     * @param index
     * @param index0
     * @param index1
     */
    SkeletonEdge(const size_t& index, const size_t& p0Index, const size_t& p1Index)
    {
        this->index = index;
        this->p0Index = p0Index;
        this->p1Index = p1Index;
    }

public:

    /**
     * @brief index
     * The index of the edge.
     */
    size_t index;

    /**
     * @brief p0Index
     * The index of the first point of the edge in the skeleton.
     */
    size_t p0Index;

    /**
     * @brief p1Index
     * The index of the seconf point of the edge in the skeleton.
     */
    size_t p1Index;
};

/**
 * @brief SkeletonEdges
 */
typedef std::vector< SkeletonEdge* > SkeletonEdges;
}
