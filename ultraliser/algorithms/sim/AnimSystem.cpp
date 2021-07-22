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

#include "AnimSystem.h"

namespace Ultraliser
{
namespace sim
{
AnimSystem::AnimSystem(float dt) : _dt(dt) {}

void AnimSystem::anim(MeshPtr mesh, uint64_t iterations)
{
    double kd = mesh->damping;
    for (uint64_t it = 0; it < iterations; ++it)
    {
        // Reset node forces to zero value
        setZeroForce(mesh);

        // Parallel compute of per spring force
        OMP_PARALLEL_FOR
        for (uint64_t i = 0; i < mesh->springs.size(); ++i)
        {
            auto spring = mesh->springs[i];
            double ks = spring->stiffness;

            Vector3f d = spring->node1->position - spring->node0->position;
            double r = spring->restLength;
            double l = d.abs();
            Vector3f v = spring->node1->velocity - spring->node0->velocity;
            Vector3f f0 = Vector3f::ZERO;
            if (l > 0.0)
            {
                f0 =
                    (ks * (l / r - 1.0) + kd * (v.dot(v, d) / (l * r))) * d / l;
            }
            spring->force = f0;
        }

        // Add spring force to nodes
        for (auto spring : mesh->springs)
        {
            spring->node0->force += spring->force;
            spring->node1->force -= spring->force;
        }

        // Update nodes velocity and position
        OMP_PARALLEL_FOR
        for (uint64_t i = 0; i < mesh->nodes.size(); ++i)
        {
            auto node = mesh->nodes[i];
            if (!node->fixed)
            {
                Vector3f v = node->velocity + node->force * _dt;
                Vector3f x = node->position + v * _dt;
                node->velocity = v;
                node->position = x;
            }
        }
    }
}

void AnimSystem::anim(Meshes meshes, uint64_t iterations)
{
    for (auto mesh : meshes)
    {
        anim(mesh, iterations);
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

void AnimSystem::setZeroForce(MeshPtr mesh) { setZeroForce(mesh->nodes); }

void AnimSystem::setZeroForce(Meshes meshes)
{
    for (auto mesh : meshes)
    {
        setZeroForce(mesh);
    }
}

}  // namespace sim

}  // namespace Ultraliser
