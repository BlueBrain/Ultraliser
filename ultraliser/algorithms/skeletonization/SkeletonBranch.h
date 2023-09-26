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

#include <algorithms/skeletonization/SkeletonNode.hh>
#include <data/common/CommonData.h>

namespace Ultraliser
{


enum BRANCH_CONNECTION_TYPE{
    ROOT,
    TERMINAL,
    ROOT_AND_TERMINAL,
};



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
    SkeletonBranch();
    ~SkeletonBranch();

public:

    /**
     * @brief hasTerminalNodes
     * Checks if the branch has terminal nodes associated with the given node indices.
     * @param node1Index
     * @param node2Index
     * @return
     */
    bool hasTerminalNodes(const size_t& node1Index, const size_t& node2Index);

    /**
     * @brief adjustDirection
     * @param frontNodeIndex
     * @param backNodeIndex
     */
    void adjustDirection(const size_t& frontNodeIndex, const size_t& backNodeIndex);

    /**
     * @brief isLoop
     * Verifies if the branch is a loop or not, i.e. starting and ending at the same node or not.
     * @return True if the branch itself is a loop, otherwise false.
     */
    bool isLoop() const;

    /**
     * @brief setValid
     * Validate the branch.
     */
    void setValid();

    /**
     * @brief setInvalid
     * Invalidate the branch.
     */
    void setInvalid();

    /**
     * @brief isValid
     * Verifies if the branch is valid or not.
     * @return True if the branch is valid, otherwise false.
     */
    bool isValid() const;

    /**
     * @brief setVisited
     * Sets the branch to be visited.
     */
    void setVisited();

    /**
     * @brief setUnvisited
     * Sets the branch to be unvisited for a specific traversal.
     */
    void setUnvisited();

    /**
     * @brief visited
     * Has the branch been visited before or not.
     * @return
     */
    bool visited() const;

    /**
     * @brief setTerminal
     * Sets the branch to be a terminal one, i.e. has no children or the last branch in a skeleton.
     */
    void setTerminal();

    /**
     * @brief unsetTermainal
     * Unsets the branch from being a terminal branch.
     */
    void unsetTermainal();

    /**
     * @brief isTerminal
     * Verifies if the branch is a terminal branch or not.
     * @return True if the branch is terminal, otherwise false.
     */
    bool isTerminal() const;

    /**
     * @brief setRoot
     * Sets the branch to be a root one, i.e. emanating from a soma.
     */
    void setRoot();

    /**
     * @brief unsetRoot
     * Resets the branch to be a non-root branch.
     */
    void unsetRoot();

    /**
     * @brief isRoot
     * Verifies if the branch is root branch, i.e. emanating from a soma, or not.
     * @return True if the branch was set to be root.
     */
    bool isRoot() const;

    /**
     * @brief setDuplicate
     * Sets the branch to be a duplicate, i.e. is not needed due to a merge operation.
     */
    void setDuplicate();

    /**
     * @brief unsetDuplicate
     * Resets the branch to be a non-duplicate branch.
     */
    void unsetDuplicate();

    /**
     * @brief isDuplicate
     * Verifies if the branch is duplicate or not.
     * @return True if the branch is duplicate, otherwise false.
     */
    bool isDuplicate() const;

    /**
     * @brief printTree
     * @param order
     */
    void printTree(size_t order = 0);

    /**
     * @brief computeLength
     * @return
     */
    float computeLength() const;

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
     * @brief t0Connections
     * Connecting branches at terminal 1.
     */
    std::vector< SkeletonBranch* > t1Connections;

    /**
     * @brief t2Connections
     * Connecting branches at termianl 2.
     */
    std::vector< SkeletonBranch* > t2Connections;




private:

    /**
     * @brief _flags
     * A set of flags represented in a BitArray such that every flag is represented by a single
     * bit not a full byte.
     */
    BitArray* _flags;

public:

    bool active = false;
};

/**
 * @brief SkeletonBranches
 */
typedef std::vector< SkeletonBranch* > SkeletonBranches;
}
