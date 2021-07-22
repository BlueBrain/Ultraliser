/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero <juanjose.garcia@epfl.ch>
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

#ifndef ULTRALISER_SIM_NODE_H
#define ULTRALISER_SIM_NODE_H

#include <math/Math.h>

namespace Ultraliser
{
namespace sim
{
/**
 * @brief The Node class
 */
class Node
{
public:
    /**
     * @brief Node
     *
     * @param position
     * @param fixed
     */
    Node(Vector3f position_, bool fixed_ = false)
        : position(position_)
        , fixed(fixed_)
        , index(0)
    {
    }

public:
    /**
     * @brief position
     */
    Vector3f position;

    /**
     * @brief velocity
     */
    Vector3f velocity;

    /**
     * @brief force
     */
    Vector3f force;

    /**
     * @brief position fixed
     */
    bool fixed;

    /**
     * @brief index
     */
    uint32_t index;
};

/**
 * @brief NodePtr
 */
typedef Node* NodePtr;

/**
 * @brief Nodes
 */
typedef std::vector<NodePtr> Nodes;

}  // namespace sim

}  // namespace Ultraliser

#endif  // ULTRALISER_SIM_NODE_H
