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

namespace Ultraliser
{

/**
 * @brief The SkeletonEdge struct
 * An edge composed of two nodes (SkeletonNodes in the segmented skeleton).
 */
struct SkeletonEdge
{
public:

    /**
     * @brief SkeletonEdge
     */
    SkeletonEdge()
    {
        this->node1 = nullptr;
        this->node2 = nullptr;
    }

    /**
     * @brief SkeletonEdge
     * @param index
     * @param index0
     * @param index1
     */
    SkeletonEdge(const int64_t& index, SkeletonNode* node1, SkeletonNode* node2)
    {
        this->index = index;
        this->node1 = node1;
        this->node2 = node2;
    }

public:

    /**
     * @brief index
     * The index of the edge. The index will be -1 to indicate that it is uninitialized.
     */
    int64_t index = -1;

    /**
     * @brief node1
     * A pointer to the first node of the edge.
     */
    SkeletonNode* node1;

    /**
     * @brief node2
     * A pointer to the second node of the edge.
     */
    SkeletonNode* node2;
};

/**
 * @brief SkeletonEdges
 */
typedef std::vector< SkeletonEdge* > SkeletonEdges;
}
