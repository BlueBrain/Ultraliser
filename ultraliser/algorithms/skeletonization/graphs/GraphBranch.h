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

#include <algorithms/skeletonization/SkeletonWeightedEdge.hh>

namespace Ultraliser
{

/**
 * @brief The GraphBranch class
 */
class GraphBranch
{
public:

    /**
     * @brief GraphBranch
     * @param index
     */
    GraphBranch(const int64_t& index)
    {
        this->index = index;
    }

    /**
     * @brief printTree
     * @param order
     */
    void printTree(size_t order = 0)
    {
        std::cout << std::string(order * 2, '-') << skeletonIndex << std::endl;
        for (size_t i = 0; i < children.size(); ++i)
        {
            children[i]->printTree(order + 1);
        }
    }

    /**
     * @brief printStatus
     */
    void printStatus()
    {
        if (this->isRoot)
        {
            std::cout << "Branch " << this->skeletonIndex << " is root" << std::endl;
        }
    }

    /**
     * @brief hasTerminalNodes
     * @param node1Index
     * @param node2Index
     * @return
     */
    bool hasTerminalNodes(const size_t& node1Index, const size_t& node2Index)
    {
        if (node1Index == firstNodeIndex && node2Index == lastNodeIndex)
        {
            return true;
        }

        if (node1Index == lastNodeIndex && node2Index == firstNodeIndex)
        {
            return true;
        }

        return false;
    }

public:

    /**
     * @brief firstNodeIndex
     */
    size_t firstNodeIndex;

    /**
     * @brief lastNodeIndex
     */
    size_t lastNodeIndex;

    /**
     * @brief firstNodeSkeletonIndex
     */
    size_t firstNodeSkeletonIndex;

    /**
     * @brief lastNodeSkeletonIndex
     */
    size_t lastNodeSkeletonIndex;

    /**
     * @brief index
     */
    int64_t index;

    /**
     * @brief skeletonIndex
     */
    size_t skeletonIndex;

    /**
     * @brief parent
     */
    GraphBranch* parent;

    /**
     * @brief children
     */
    std::vector< GraphBranch* > children;

    /**
     * @brief isRoot
     */
    bool isRoot = false;

    /**
     * @brief active
     */
    bool active = false;
};

/**
 * @brief GraphBranches
 */
typedef std::vector< GraphBranch* > GraphBranches;

}
