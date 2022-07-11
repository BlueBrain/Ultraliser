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

#include "Spring.h"

namespace Ultraliser
{
namespace Simulation
{

Spring::Spring(NodePtr n0, NodePtr n1, float stiff, float rLength)
    : node0(n0)
    , node1(n1)
    , stiffness(stiff)
    , restLength(rLength)
{
    if (restLength == 0.f)
    {
        restLength = length();
    }
}

float Spring::length() const
{
    // Compute the current length of the spring
    return (node0->position - node1->position).abs();
}

void Spring::computeForce()
{
    // Compute damping constant from stiffness constant
    double kd = stiffness * 0.01f;

    // Compute spring nodes distance
    Vector3f d = node1->position - node0->position;

    // Compute nodes distance norm
    double l = d.abs();

    // Compute nodes velocities difference to apply damping force
    Vector3f v = node1->velocity - node0->velocity;

    // Compute per spring force as the sum of stress and damping components
    force = Vector3f::ZERO;
    if (l > 0.0f)
    {
        force = (stiffness * (l / restLength - 1.0f) +
                 kd * (v.dot(v, d) / (l * restLength))) * d / l;
    }
}

size_t SpringHash::operator()(const SpringPtr spring) const
{
    uint64_t id0 = spring->node0->index;
    uint64_t id1 = spring->node1->index;

    if (id1 > id0)
        std::swap(id1, id0);

    return std::hash<unsigned int>{}(id0) ^ std::hash<unsigned int>{}(id1);
}

bool SpringEqual::operator()(const SpringPtr spring0, const SpringPtr spring1) const
{
    uint64_t id00 = spring0->node0->index;
    uint64_t id01 = spring0->node1->index;
    uint64_t id10 = spring1->node0->index;
    uint64_t id11 = spring1->node1->index;

    if (id01 > id00)
        std::swap(id00, id01);
    if (id11 > id10)
        std::swap(id10, id11);

    return (id00 == id10) && (id01 == id11);
}

}
}
