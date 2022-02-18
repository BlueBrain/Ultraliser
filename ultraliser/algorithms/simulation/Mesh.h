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

#ifndef ULTRALISER_ALGORITHMS_SIMULATION_MESH_H
#define ULTRALISER_ALGORITHMS_SIMULATION_MESH_H

#include "StiffnessMatrix.h"

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
     */ 
    Mesh();

    /**
     * @brief Mesh
     * Constructor
     *
     * @param nodes
     * Mesh nodes
     * @param surfaceNodes
     * Mesh surface nodes
     * @param springs
     * Mesh springs
     * @param tetrahedra
     * Mesh tetrahedra
     */
    Mesh(Nodes nodes, Nodes surfaceNodes, Springs springs, Tetrahedra tetrahedra);

    /**
     * @brief ~Mesh
     * Destructor
     */
    ~Mesh();
    

    void computeStiffnessMatrix( float stiffness = 10000.0,
                                 float poissonRatio = 0.3,
                                 float dt = 0.01);

public:
    /**
     * @brief nodes
     */
    Nodes nodes;

    /**
     * @brief surface nodes
     */
    Nodes surfaceNodes;

    /**
     * @brief springs
     */
    Springs springs;

    /**
     * @brief tetrahedra
     */
    Tetrahedra tetrahedra;

    /**
     * @brief stiffnessMatrix
     */
    StiffnessMatrixPtr stiffnessMatrix;

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
