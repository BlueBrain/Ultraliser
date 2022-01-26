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

#include "Tetrahedron.h"

namespace Ultraliser
{
namespace Simulation
{

Tetrahedron::Tetrahedron(NodePtr n0, NodePtr n1, NodePtr n2, NodePtr n3)
    : node0(n0)
    , node1(n1)
    , node2(n2)
    , node3(n3)
{
    _initVolume = volume();
}

float Tetrahedron::volume() const
{
    auto x0 = node0->position;
    auto x1 = node1->position;
    auto x2 = node2->position;
    auto x3 = node3->position;
    Matrix3f basis (x1 - x0, x2 - x0, x3 - x0, true);

    return std::abs(basis.determinant() / 6.0);
}

float Tetrahedron::initVolume() const
{
    return _initVolume;
}

}
}
