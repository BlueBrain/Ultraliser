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
