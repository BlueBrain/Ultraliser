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




#define ULTRALISER_VALIDAITY_BIT_INDEX 0
#define ULTRALISER_TERMINAL_BIT_INDEX 1

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
        _flags = new BitArray(8);
        _flags->clearAll();
    }

    ~SkeletonBranch() { _flags->~BitArray(); }

    void printTree(size_t order = 0)
    {
        std::cout << std::string(order * 4, '-') << index << "\n";
        for (size_t i = 0; i < children.size(); ++i)
        {
            children[i]->printTree(order + 1);
        }
    }

    /**
     * @brief hasTerminalNodes
     * Checks if the branch has terminal nodes associated with the given node indices.
     * @param node1Index
     * @param node2Index
     * @return
     */
    bool hasTerminalNodes(const size_t& node1Index, const size_t& node2Index)
    {
        const auto& frontNode = nodes.front();
        const auto& backNode = nodes.back();

        if (node1Index == backNode->index && node2Index == frontNode->index)
        {
            return true;
        }

        if (node1Index == frontNode->index && node2Index == backNode->index)
        {
            return true;
        }

        return false;
    }

    /**
     * @brief adjustDirection
     * @param frontNodeIndex
     * @param backNodeIndex
     */
    void adjustDirection(const size_t& frontNodeIndex, const size_t& backNodeIndex)
    {
        if (nodes.front()->index == backNodeIndex && nodes.back()->index == frontNodeIndex)
        {
            std::reverse(this->nodes.begin(), this->nodes.end());
        }
    }

public:

    void setValid() { _flags->setBit(0); }

    void setInvalid() { _flags->clearBit(0); }

    bool isValid() const { return _flags->bit(0); }

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
    BitArray* _flags;

public:

    /**
     * @brief root
     * If this flag is set, this means that it is emanating from the soma
     */
    bool root = false;

    bool terminal = false;

    bool valid = true;

    bool duplicate = false;

    bool visited = false;

    bool active = false;
};

/**
 * @brief SkeletonBranches
 */
typedef std::vector< SkeletonBranch* > SkeletonBranches;
}
