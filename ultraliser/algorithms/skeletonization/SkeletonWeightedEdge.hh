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

#include <common/Headers.hh>
#include <algorithms/skeletonization/SkeletonNode.hh>
#include <algorithms/skeletonization/SkeletonBranch.h>

namespace Ultraliser
{

/**
 * @brief The SkeletonWeightedEdge class
 */
struct SkeletonWeightedEdge
{
public:

    SkeletonWeightedEdge(SkeletonBranch* branch)
    {
        this->branch = branch;
        this->node1 = branch->nodes.front();
        this->node2 = branch->nodes.back();
        this->edgeWeight = branch->nodes.size();
    }

    bool hasTerminalGraphNodes(const int64_t node1Index, const int64_t node2Index)
    {
        if (node1->graphIndex == node1Index && node2->graphIndex == node2Index)
        {
            return true;
        }

        if (node1->graphIndex == node2Index && node2->graphIndex == node1Index)
        {
            return true;
        }

        return false;
    }

    bool visited = false;

public:

    /**
     * @brief branch
     */
    SkeletonBranch* branch;

    /**
     * @brief node1
     */
    SkeletonNode* node1;

    /**
     * @brief node2
     */
    SkeletonNode* node2;

    /**
     * @brief edgeWeight
     * The weight of the edge. This weight could be positive or negative.
     */
    int64_t edgeWeight;
};

/**
 * @brief WeightedEdges
 * List of WeightedEdge elements.
 */
typedef std::vector< SkeletonWeightedEdge* > SkeletonWeightedEdges;

}
