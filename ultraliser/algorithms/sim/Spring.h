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

#ifndef ULTRALISER_SIM_SPRING_H
#define ULTRALISER_SIM_SPRING_H

#include "Node.h"

namespace Ultraliser
{
namespace sim
{
/**
 * @brief The Spring class
 */
class Spring
{
public:
    /**
     * @brief Spring
     * Constructor.
     * @param position
     * @param fixed
     */
    Spring(NodePtr node0,
           NodePtr node1,
           float stiffness,
           float restLength = .0f);

    /**
     * @brief length
     */
    float length() const;

public:
    /**
     * @brief node0
     */
    NodePtr node0;

    /**
     * @brief node1
     */
    NodePtr node1;

    /**
     * @brief spring stiffness
     */
    float stiffness;

    /**
     * @brief spring rest length
     */
    float restLength;

    /**
     * @brief spring force
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
    size_t operator()(const SpringPtr spring) const;
};

/**
 * @brief The Spring Equal class
 */
class SpringEqual
{
public:
    bool operator()(const SpringPtr spring0, const SpringPtr spring1) const;
};

/**
 * @brief UniqueSprings
 */
typedef std::unordered_set<SpringPtr, SpringHash, SpringEqual> UniqueSprings;

}  // namespace sim

}  // namespace Ultraliser

#endif  // ULTRALISER_SIM_SPRING_H
