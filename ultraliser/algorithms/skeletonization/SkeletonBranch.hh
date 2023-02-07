#pragma once

#include <algorithms/skeletonization/SkeletonNode.hh>

namespace Ultraliser
{

/**
 * @brief The SkeletonBranch struct
 * A branch in the segmented skeleton.
 */
struct SkeletonBranch
{
public:

    /**
     * @brief SkeletonBranch
     */
    SkeletonBranch()
    {
        /// EMPTY CONSTRUCTOR
    }

public:

    /**
     * @brief nodes
     * A list of nodes in the skeleton.
     */
    SkeletonNodes nodes;

    /**
     * @brief index
     * The index of the branch in the skeleton.
     */
    size_t index;

    /**
     * @brief parent
     * A list of the parent branches.
     */
    std::vector< SkeletonBranch* > parents;

    /**
     * @brief children
     * A list of the child branches.
     */
    std::vector< SkeletonBranch* > children;

    /**
     * @brief root
     * If this flag is set, this means that it is emanating from the soma
     */
    bool root = false;

    bool valid = true;

    bool duplicate = false;




    SkeletonNode* replacingNode = nullptr;
};

/**
 * @brief SkeletonBranches
 */
typedef std::vector< SkeletonBranch* > SkeletonBranches;
}
