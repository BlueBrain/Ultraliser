/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#pragma once

#include "Node.h"

namespace Ultraliser
{
namespace Simulation
{

/**
 * @brief The Tetrahedron class
 */
class Tetrahedron
{
public:

    /**
     * @brief Tetrahedron
     * Constructor
     *
     * @param node0
     * The first node of the tetrahedron.
     * @param node1
     * The second node of the tetrahedron.
     * @param node2
     * The third node of the tetrahedron.
     * @param node3
     * The fourth node of the tetrahedron.
     */
    Tetrahedron(NodePtr node0, NodePtr node1, NodePtr node2, NodePtr node3);

    /**
     * @brief volume
     * Gets the actual volume of the tetrahedron.
     */
    float volume() const;

    /**
     * @brief initVolume
     * Gets the initial volume of the tetrahedron.
     */
    float initVolume() const;

public:

    /**
     * @brief node0
     * The first node of the tetrahedron.
     */
    NodePtr node0;

    /**
     * @brief node1
     * The second node of the tetrahedron.
     */
    NodePtr node1;

    /**
     * @brief node2
     * The third node of the tetrahedron.
     */
    NodePtr node2;

    /**
     * @brief node3
     * The fourth node of the tetrahedron.
     */
    NodePtr node3;

private:

    /**
     * @brief _initVolume
     * The initial volume of the Tetrahedron.
     */
    float _initVolume;
};

/**
 * @brief TetrahedronPtr
 */
typedef Tetrahedron* TetrahedronPtr;

/**
 * @brief Tetrahedra
 */
typedef std::vector<TetrahedronPtr> Tetrahedra;

}
}
