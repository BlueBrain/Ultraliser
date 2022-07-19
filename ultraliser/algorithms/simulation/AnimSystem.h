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


#pragma once

#include "Mesh.h"

namespace Ultraliser
{
namespace Simulation
{

/**
 * @brief The AnimSystem class
 */
class AnimSystem
{
public:

    /**
     * @brief AnimSystem
     * Constructor
     *
     * @param dt
     * Simulation time-step.
     */
    AnimSystem(float dt);

    /**
     * @brief anim
     *
     * @param mesh
     */
    void animate(MeshPtr mesh, size_t iterations = 1);

    /**
     * @brief anim
     *
     * @param meshes
     */
    void animate(Meshes meshes, size_t iterations = 1);

    /**
     * @brief setZeroForce
     *
     * @param nodes
     */
    void setZeroForce(Nodes& nodes);

    /**
     * @brief setZeroForce
     *
     * @param mesh
     */
    void setZeroForce(MeshPtr mesh);

    /**
     * @brief setZeroForce
     *
     * @param meshes
     */
    void setZeroForce(Meshes meshes);

private:

    /**
     * @brief _dt
     * Animation time-step.
     */
    float _dt;
};

}
}
