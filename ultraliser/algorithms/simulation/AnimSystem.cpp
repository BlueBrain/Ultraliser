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

#include "AnimSystem.h"

namespace Ultraliser
{
namespace Simulation
{

AnimSystem::AnimSystem(float dt)
    : _dt(dt)
{
    /// EMPTY CONSTRUCTOR
}

void AnimSystem::animate(MeshPtr mesh, uint64_t iterations)
{
    for (uint64_t it = 0; it < iterations; it++)
    {
        // Reset node forces to zero value
        setZeroForce(mesh);

        // Parallel compute of per spring force
        OMP_PARALLEL_FOR
        for (uint64_t i = 0; i < mesh->springs.size(); ++i)
        {
            mesh->springs[i]->computeForce();
        }

        // Add spring force to nodes out of the parallel loop to avoid concurrency problems
        for (auto spring : mesh->springs)
        {
            spring->node0->force += spring->force;
            spring->node1->force -= spring->force;
        }

        // Update nodes velocity and position applying an symplectic integration scheme
        OMP_PARALLEL_FOR
        for (uint64_t i = 0; i < mesh->nodes.size(); ++i)
        {
            auto node = mesh->nodes[i];
            if (!node->fixed)
            {
                node->velocity = node->velocity + node->force * _dt;
                node->position = node->position + node->velocity * _dt;
            }
        }
    }
}

void AnimSystem::animate(Meshes meshes, uint64_t iterations)
{
    for (auto mesh : meshes)
    {
        animate(mesh, iterations);
    }
}

void AnimSystem::setZeroForce(Nodes& nodes)
{
    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->force = Vector3f::ZERO;
    }
}

void AnimSystem::setZeroForce(MeshPtr mesh)
{
    setZeroForce(mesh->nodes);
}

void AnimSystem::setZeroForce(Meshes meshes)
{
    for (auto mesh : meshes)
    {
        setZeroForce(mesh);
    }
}

}
}
