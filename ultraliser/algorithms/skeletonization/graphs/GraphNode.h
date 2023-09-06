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

class GraphNode
{
public:

    GraphNode(int64_t index, size_t skeletonIndex = 0)
    {
        this->index = index;
        this->skeletonIndex = skeletonIndex;
    }

    int64_t index = -1;
    size_t skeletonIndex;

    std::vector<GraphNode*> parents;
    std::vector<GraphNode*> children;

    bool isNodeInChildren(int64_t childNodeIndex)
    {
        for (size_t i = 0; i < children.size(); ++i)
        {
            if (childNodeIndex == children[i]->index)
            {
                return true;
            }
        }
        return false;
    }
};

}
