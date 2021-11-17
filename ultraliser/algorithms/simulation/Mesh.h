/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#ifndef ULTRALISER_ALGORITHMS_SIMULATION_MESH_H
#define ULTRALISER_ALGORITHMS_SIMULATION_MESH_H

#include "Spring.h"

namespace Ultraliser
{
namespace Simulation
{
/**
 * @brief The Mesh class
 */
class Mesh
{
public:

    /**
     * @brief Mesh
     * Constructor
     *
     * @param stiff
     * Spring stiffness
     * @param damp
     * Spring damping
     */
    Mesh(float stiff = 1000.f, float damp = 1.f)
        : stiffness(stiff)
        , damping(damp){};

    /**
     * @brief Mesh
     * Constructor
     *
     * @param nodes
     * Mesh nodes
     * @param springs
     * Mesh springs
     */
    Mesh(Nodes nodes, Springs springs, float stiffnes = 1000.f, float damping = 1.f)
        : Mesh(stiffnes, damping)
    {
        this->nodes = nodes;
        this->springs = springs;
    };

    /**
     * @brief ~Mesh
     * Destructor
     */
    ~Mesh()
    {
        for (auto node : nodes)
            delete node;

        for (auto spring : springs)
            delete spring;
    }

public:
    /**
     * @brief nodes
     */
    Nodes nodes;

    /**
     * @brief springs
     */
    Springs springs;

    /**
     * @brief stiffness
     */
    float stiffness;

    /**
     * @brief damping
     */
    float damping;
};

/**
 * @brief MeshPtr
 */
typedef Mesh* MeshPtr;

/**
 * @brief Meshes
 */
typedef std::vector<MeshPtr> Meshes;

}
}

#endif  // ULTRALISER_ALGORITHMS_SIMULATION_MESH_H
