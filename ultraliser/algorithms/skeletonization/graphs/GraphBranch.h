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

class GraphBranch
{


public:
    GraphBranch(int64_t index)
    {
        this->index = index;
    }

    void printTree(size_t order = 0)
    {
        std::cout << std::string(order * 4, '-') << skeletonIndex << "\n";
        for (size_t i = 0; i < children.size(); ++i)
        {
            children[i]->printTree(order + 1);
        }
    }

    void printStatus()
    {
        if (this->isRoot)
        {
            std::cout << "Branch " << this->skeletonIndex << " is Root \n";
        }
    }

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

    size_t firstNodeIndex;
    size_t lastNodeIndex;

    size_t firstNodeSkeletonIndex;
    size_t lastNodeSkeletonIndex;

    int64_t index;
    size_t skeletonIndex;

    GraphBranch* parent;
    std::vector< GraphBranch* > children;

    bool isRoot = false;
    bool active = false;



};

typedef std::vector< GraphBranch* > GraphBranches;

}
