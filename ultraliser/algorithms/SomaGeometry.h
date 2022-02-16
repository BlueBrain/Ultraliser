/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_SOMA_GEOMETRY_H
#define ULTRALISER_SOMA_GEOMETRY_H

#include <data/meshes/simple/primitives/Primitives.h>
#include <data/morphologies/NeuronMorphology.h>
#include <algorithms/simulation/AnimSystem.h>

namespace Ultraliser
{

/**
 * @brief The SomaGeometry class
 */
class SomaGeometry
{
public:

    /**
     * @brief SomaGeometry
     * Constructor
     *
     * @param morphology
     * The given neuron moprhology to generate the somatic mesh from.
     * @param stiffness
     * Simulation stiffness factor. Takes values between 0 and 1.
     * @param dt
     * Simulation time increment.
     * @param numIterations
     * Nunber of simulation iterations.
     */
    SomaGeometry(NeuronMorphology* morphology, float stiffness = 1.0f, float dt = 0.01f, 
                 uint32_t numIterations = 8000);

private:

    /**
     * @brief _insert
     * Insert spring in the unique springs set.
     *
     * @param spring
     * The given spring to be added.
     * @param springs
     * A set of springs where the given spring will be added.
     */
    void _insert(Simulation::SpringPtr spring, Simulation::UniqueSprings& springs);

    /**
     * @brief _loadIcosphereGeometry
     * Load geoemtry from icosphere data.
     *
     * @return MeshPtr
     * A pointer to the simulation mesh.
     */
    Simulation::MeshPtr _loadIcosphereGeometry();

    /**
     * @brief _nodesToVertices
     * Dumps nodes' information to vertices and triangles.
     */
    void _nodesToVertices(Simulation::Nodes& nodes);

    /**
     * @brief _generateIcosphereSprings
     * Add springs based in the coarse icospehere geometry.
     *
     * @param mesh
     * Simulation mesh.
     * @param stiffness
     * Simulation stiffness
     * @return 
     * Icosphere springs.
     */
     Simulation::Springs _generateIcosphereSprings(Simulation::MeshPtr mesh, float stiffness);

    /**
     * @brief _generatePullSprings
     * Adds springs incharge to pull in the fisrSections directions.
     *
     * @param mesh
     * Simulation mesh.
     * @param stiffness
     * Simulation stiffness.
     * @param somaSamples
     * A list of the soma samples.
     * @param restLengthThreshold
     * The threshold length of the springs at the resting state.
     * @return 
     * Pull springs.
     */
    Simulation::Springs _generatePullSprings(Simulation::MeshPtr mesh,
                                             float stiffness,
                                             const Samples& somaSamples,
                                             float restLengthThreshold = 0.2);

    /**
     * @brief _fixCenterNode
     * Fixes the center nodes of the initial icosphere based on the threshold.
     *
     * @param nodes
     * Ico-sphere nodes list.
     * @param threshold
     * A given threshold to fix the nodes at.
     */
    void _fixCenterNodes(Simulation::Nodes& nodes, float threshold);

public:

    /**
     * @brief vertices
     * A list of all the vertices in the soma geometry.
     */
    Vertex* vertices;

    /**
     * @brief numVertices
     * Number of vertices of the section geometry.
     */
    uint64_t numVertices;

    /**
     * @brief triangles
     * A list of all the triangles in the soma geometry.
     */
    Triangle* triangles;

    /**
     * @brief numTriangles
     * Number of triangles of the section geometry.
     */
    uint64_t numTriangles;

private:

    /**
     * @brief _somaCenter
     * The center of the soma. Typically this is the origin, unless otherwise specified.
     */
    Vector3f _somaCenter;

    /**
     * @brief _somaRadius;
     * The average radius of the soma.
     */
    float _somaRadius;
};

}  // namespace Ultraliser

#endif  // ULTRALISER_SOMA_GEOMETRY_H
