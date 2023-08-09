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
     * @brief t0Connections
     * Connecting branches at terminal 1.
     */
    std::vector< SkeletonBranch* > t1Connections;

    /**
     * @brief t2Connections
     * Connecting branches at termianl 2.
     */
    std::vector< SkeletonBranch* > t2Connections;

    /**
     * @brief root
     * If this flag is set, this means that it is emanating from the soma
     */
    bool root = false;

    bool terminal = false;

    bool valid = true;

    bool duplicate = false;
};

/**
 * @brief SkeletonBranches
 */
typedef std::vector< SkeletonBranch* > SkeletonBranches;
}
