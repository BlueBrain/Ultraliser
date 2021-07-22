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

#ifndef ULTRALISER_SOMA_GEOMETRY_H
#define ULTRALISER_SOMA_GEOMETRY_H

#include <data/meshes/simple/primitives/Primitives.h>
#include <data/morphologies/NeuronMorphology.h>

#include "sim/AnimSystem.h"

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
     * Constructor.
     * @param morphology
     * neuron moprhology to generate soma mesh from.
     */
    SomaGeometry(NeuronMorphology* morphology);

private:
    /**
     * @brief _insert
     * Insert spring in the unique springs set.
     * @param spring
     * @param springs
     * @return
     */
    void _insert(sim::SpringPtr spring, sim::UniqueSprings& springs);

    /**
     * @brief _loadIcosphereGeometry
     * Load geoemtry from icosphere data.
     * @return MeshPtr
     */
    sim::MeshPtr _loadIcosphereGeometry();

    /**
     * @brief _nodesToVertices
     * Dump nodes information to vertices and triangles.
     * @return
     */
    void _nodesToVertices(sim::Nodes& nodes);

    /**
     * @brief _generateIcosphereSprings
     * Add springs based in the coarse icospehere geometry.
     * @param mesh
     * @param stiffness
     * @return
     */
    void _generateIcosphereSprings(sim::MeshPtr mesh, float stiffness);

    /**
     * @brief _generatePullSprings
     * Add springs incharge to pull in the fisrSections directions.
     * @param mesh
     * @param stiffness
     * @param firstSecions
     * @param restLengthThreshold
     * @return
     */
    void _generatePullSprings(sim::MeshPtr mesh,
                              float stiffness,
                              const Sections& firstSections,
                              float restLengthThreshold = 0.2);

    /**
     * @brief _fixCenterNode
     * Fixes the center nodes of the initial icosphere based on the threshold.
     * @param nodes
     * @param threshold
     * @return
     */
    void _fixCenterNodes(sim::Nodes& nodes, float threshold);

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
     */
    Vector3f _somaCenter;

    /**
     * @brief _somaRadius;
     */
    float _somaRadius;
};

}  // namespace Ultraliser

#endif  // ULTRALISER_SOMA_GEOMETRY_H
