/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Juan Jose Garcia Cantero <juanjose.garcia@epfl.ch>
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

#ifndef ULTRALISER_ALGORITHMS_SIMULATION_NODE_H
#define ULTRALISER_ALGORITHMS_SIMULATION_NODE_H

#include <math/Math.h>

namespace Ultraliser
{
namespace Simulation
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
    Node(Vector3f pos, float mass = 1.f, bool fix = false)
        : _initPosition(pos)
        , position(pos)
        , mass(mass)
        , fixed(fix)
        , index(0)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief initPosition
     * Gets node initial position
     */
    Vector3f initPosition() const
    {
        return _initPosition;
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
     * @brief mass
     */
    float mass;

    /**
     * @brief position fixed
     */
    bool fixed;

    /**
     * @brief index
     */
    uint32_t index;

private:

    /**
     * @brief _initPosition
     * Node initial position
     */
    Vector3f _initPosition;

};

/**
 * @brief NodePtr
 */
typedef Node* NodePtr;

/**
 * @brief Nodes
 */
typedef std::vector< NodePtr > Nodes;

}
}

#endif  // ULTRALISER_ALGORITHMS_SIMULATION_NODE_H
