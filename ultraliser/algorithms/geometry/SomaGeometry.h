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

#pragma once

#include <data/meshes/simple/primitives/Primitives.h>
#include <data/morphologies/NeuronMorphology.h>
#include <data/morphologies/AstrocyteMorphology.h>

#ifdef ULTRALISER_USE_EIGEN3
#include <algorithms/simulation/AnimSystem.h>
#else 
#include <algorithms/simulation/Mesh.h>
#endif

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
     * @param alphaRadius
     * Multiplier of the soma radius.
     * @param stiffness
     * Simulation stiffness factor. Takes values over 0.
     * @param poissonRatio
     * Simulation poisson ratio. Takes values between (0, 0.5].
     * @param dt
     * Simulation time increment.
     * @param numIterations
     * Nunber of simulation iterations.
     */
    SomaGeometry(const NeuronMorphology* morphology, const float& alphaRadius = 0.75f, 
                 const float& stiffness = 10000.0f, const float& poissonRatio = 0.2f, 
                 const float& dt = 0.01f, const uint32_t& numIterations = 200);


    SomaGeometry(const AstrocyteMorphology* morphology, const float& alphaRadius = 0.75f, 
                 const float& stiffness = 10000.0f, const float& poissonRatio = 0.2f, 
                 const float& dt = 0.01f, const uint32_t& numIterations = 200);

private:

    /**
     * @brief _somaGeometryGeneration
     * Generates the soma geoemtry.
     * 
     * @param somaSamples
     * Soma samples.
     * @param alphaRadius
     * Multiplier of the soma radius.
     * @param stiffness
     * Simulation stiffness factor. Takes values over 0.
     * @param poissonRatio
     * Simulation poisson ratio. Takes values between (0, 0.5].
     * @param dt
     * Simulation time increment.
     * @param numIterations
     * Nunber of simulation iterations.
     */
    void _somaGeometryGeneration(const Samples& somaSamples, const float& stiffness,
                                 const float& poissonRatio, const float& dt,
                                 const uint32_t& numIterations);

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
     * @brief _computeNeuritesNodes
     * Compute the icosphere nodes linkd to the neurites starts 
     *
     * @param mesh
     * Simulation mesh.
     * @param neuritesNodes
     * Vector of nodes linked to each neurite start.
     * @param surfacePositions
     * Position of the surfce linked to each neurite start
     * @param somaSamples
     * A list of the soma samples.
     */
    void _computeNeuritesNodes(Simulation::MeshPtr mesh,
                               std::vector<Simulation::Nodes>& neuritesNodes,
                               std::vector<Vector3f>& surfacePositions,
                               const Samples& somaSamples);

    /**
     * @brief _pullNeuritesNodes
     * Reposition the nodes linked to the neurites starts in their direction multiplied by alpha
     *
     * @param neuritesNodes
     * Vector of nodes linked to each neurite start.
     * @param surfacePositions
     * Position of the surfce linked to each neurite start
     * @param somaSamples
     * A list of the soma samples.
     * @param alpha
     * Direction multiplicator parameter
     */
    void _pullNeuritesNodes(std::vector<Simulation::Nodes>& neuritesNodes,
                            std::vector<Vector3f>& surfacePositions,
                            const Samples& somaSamples, float alpha = 0.01f);

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
    size_t numVertices;

    /**
     * @brief triangles
     * A list of all the triangles in the soma geometry.
     */
    Triangle* triangles;

    /**
     * @brief numTriangles
     * Number of triangles of the section geometry.
     */
    size_t numTriangles;

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

}
