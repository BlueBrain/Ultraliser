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

#include "Spring.h"

namespace Ultraliser
{
namespace sim
{
Spring::Spring(NodePtr node0, NodePtr node1, float stiffness, float restLength_)
    : node0(node0)
    , node1(node1)
    , stiffness(stiffness)
    , restLength(restLength_)
{
    if (restLength == .0f)
    {
        restLength = length();
    }
}

float Spring::length() const
{
    return (node0->position - node1->position).abs();
}

size_t SpringHash::operator()(const SpringPtr spring) const
{
    uint64_t id0 = spring->node0->index;
    uint64_t id1 = spring->node1->index;
    if (id1 > id0) std::swap(id1, id0);
    return std::hash<unsigned int>{}(id0) ^ std::hash<unsigned int>{}(id1);
}

bool SpringEqual::operator()(const SpringPtr spring0,
                             const SpringPtr spring1) const
{
    uint64_t id00 = spring0->node0->index;
    uint64_t id01 = spring0->node1->index;
    uint64_t id10 = spring1->node0->index;
    uint64_t id11 = spring1->node1->index;
    if (id01 > id00) std::swap(id00, id01);
    if (id11 > id10) std::swap(id10, id11);
    return (id00 == id10) && (id01 == id11);
}

}  // namespace sim

}  // namespace Ultraliser
