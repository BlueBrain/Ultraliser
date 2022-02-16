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

#ifndef ULTRALISER_ALGORITHMS_SIMULATION_SPRING_H
#define ULTRALISER_ALGORITHMS_SIMULATION_SPRING_H

#include "Node.h"

namespace Ultraliser
{
namespace Simulation
{

/**
 * @brief The Spring class
 */
class Spring
{
public:

    /**
     * @brief Spring
     * Constructor
     *
     * @param node0
     * The first node of the spring.
     * @param node1
     * The second node of the spring.
     * @param stiffness
     * The stiffness of the spring.
     * @param restLength
     * The resting length of the spring.
     */
    Spring(NodePtr node0, NodePtr node1, float stiffness, float restLength = 0.f);

    /**
     * @brief length
     * Gets the actual length of the spring.
     */
    float length() const;

    /**
     * @brief computeForce
     * Compute spring force based on nodes positions and velocities
     */
    void computeForce();

public:

    /**
     * @brief node0
     * The first node of the spring.
     */
    NodePtr node0;

    /**
     * @brief node1
     * The second node of the spring.
     */
    NodePtr node1;

    /**
     * @brief spring stiffness
     * The stiffness of the spring.
     */
    float stiffness;

    /**
     * @brief restLength
     * The resting length of the spring.
     */
    float restLength;

    /**
     * @brief spring force
     * The current force applied to the spring.
     */
    Vector3f force;
};

/**
 * @brief SpringPtr
 */
typedef Spring* SpringPtr;

/**
 * @brief Springs
 */
typedef std::vector<SpringPtr> Springs;

/**
 * @brief The Spring Hash class
 */
class SpringHash
{
public:

    /**
     * @brief operator ()
     * @param spring
     * @return
     */
    size_t operator()(const SpringPtr spring) const;
};

/**
 * @brief The Spring Equal class
 */
class SpringEqual
{
public:

    /**
     * @brief operator ()
     * Verifies if the springs are equal or not.
     * @param spring0
     * The first spring.
     * @param spring1
     * The second spring.
     * @return
     * True if the two given springs are equal, otherwise false.
     */
    bool operator()(const SpringPtr spring0, const SpringPtr spring1) const;
};

/**
 * @brief UniqueSprings
 */
typedef std::unordered_set<SpringPtr, SpringHash, SpringEqual> UniqueSprings;

}
}

#endif  // ULTRALISER_ALGORITHMS_SIMULATION_SPRING_H
