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

#include "Mesh.h"

namespace Ultraliser
{
namespace Simulation
{

Mesh::Mesh(){}

Mesh::Mesh(Nodes nodes, Nodes surfaceNodes, Springs springs, Tetrahedra tetrahedra)
    : nodes(nodes)
    , surfaceNodes(surfaceNodes)
    , springs(springs)
    , tetrahedra(tetrahedra){}

Mesh::~Mesh()
{
    for (auto node : nodes)
        delete node;
    for (auto node : surfaceNodes)
        delete node;
    for (auto spring : springs)
        delete spring;
    for (auto tetrahedron : tetrahedra)
        delete tetrahedron;
}

void Mesh::computeStiffnessMatrix(float stiffness, float poissonRatio, float dt)
{
    stiffnessMatrix = new StiffnessMatrix(nodes, tetrahedra, stiffness, poissonRatio, dt);
}

}
}
