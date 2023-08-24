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

namespace Ultraliser
{

/**
 * @brief The WeightedEdge class
 */
struct WeightedEdge
{
public:

    /**
     * @brief WeightedEdge
     * Constructor
     * @param node1Index
     * The index of the first node of the edge.
     * @param node2Index
     * The index of the second node of the edge.
     * @param edgeWeight
     * The weight of the edge. This weight could be positive or negative.
     */
    WeightedEdge(const size_t& node1Index,
                 const size_t& node2Index,
                 const int64_t& edgeWeight)
    {
        this->node1Index = node1Index;
        this->node2Index = node2Index;
        this->edgeWeight = edgeWeight;
    }

public:

    /**
     * @brief node1Index
     * The index of the first node of the edge.
     */
    size_t node1Index;

    /**
     * @brief node2Index
     * The index of the second node of the edge.
     */
    size_t node2Index;

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
typedef std::vector< WeightedEdge* > WeightedEdges;

}
