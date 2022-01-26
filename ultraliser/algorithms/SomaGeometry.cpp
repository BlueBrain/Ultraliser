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

#include "SomaGeometry.h"

#include "IcosphereGeometry.hh"

namespace Ultraliser
{
SomaGeometry::SomaGeometry(NeuronMorphology* morphology, float stiffness, float poissonRatio,
                           float dt, uint32_t numIterations)
    : numVertices(0)
    , numTriangles(0)
{
    if (morphology != nullptr)
    {
        _somaCenter = morphology->getSomaCenter();

        // For the moment, use the min soma radius to see the effect of the pulling forces on
        // the ico-sphere till further notice
        _somaRadius = morphology->getSomaMinRadius();

        auto mesh = _loadIcosphereGeometry();
        mesh->computeStiffnessMatrix(stiffness, poissonRatio, dt);

        // Compute the simulation nodes linked to neurites starts
        std::vector<Simulation::Nodes> neuritesNodes;
        std::vector<Vector3f> surfacePositions;
        auto neurites = morphology->getSomaSamples();
        _computeNeuritesNodes(mesh, neuritesNodes, surfacePositions, neurites);
        
        
        Simulation::AnimSystem animSystem(dt);

        // Perform the icosphere animation adapting the stiffness of the pull springs
        for (uint32_t i = 1; i <= numIterations; ++i)
        {
            float alpha = (float)i/ numIterations;
            _pullNeuritesNodes(neuritesNodes, surfacePositions, neurites, alpha);
            animSystem.animate(mesh);
        }

        _nodesToVertices(mesh->nodes);
    }
    else
    {
        LOG_ERROR("The input morphology is EMPTY!");
    }
}

Simulation::MeshPtr SomaGeometry::_loadIcosphereGeometry()
{
    auto mesh = new Simulation::Mesh();
    Simulation::Nodes& nodes = mesh->nodes;
    nodes.resize(IcosphereVerticesSize);

    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < nodes.size(); ++i)
    {
        Vector3f position(IcosphereVertices[i * 3 + 0],
                          IcosphereVertices[i * 3 + 1],
                          IcosphereVertices[i * 3 + 2]);
        nodes[i] = new Simulation::Node(position * _somaRadius + _somaCenter);
        nodes[i]->index = i;
    }

    Simulation::Tetrahedra& tetrahedra = mesh->tetrahedra;
    tetrahedra.resize(IcosphereTetsIndicesSize);

    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < tetrahedra.size(); ++i)
    {
        auto node0 = nodes[IcosphereTetsIndices[i * 4 + 0]];
        auto node1 = nodes[IcosphereTetsIndices[i * 4 + 1]];
        auto node2 = nodes[IcosphereTetsIndices[i * 4 + 2]];
        auto node3 = nodes[IcosphereTetsIndices[i * 4 + 3]];
        tetrahedra[i] = new Simulation::Tetrahedron(node0, node1, node2, node3);
    }

    return mesh;
}

void SomaGeometry::_nodesToVertices(Simulation::Nodes& nodes)
{
    // Get surface mesh nodes
    std::unordered_set<Simulation::NodePtr> uniqueNodes;
    for (uint32_t i = 0; i < IcosphereTrianglesIndicesSize; ++i)
    {
        uniqueNodes.insert(nodes[IcosphereTrianglesIndices[i * 3]]);
        uniqueNodes.insert(nodes[IcosphereTrianglesIndices[i * 3 + 1]]);
        uniqueNodes.insert(nodes[IcosphereTrianglesIndices[i * 3 + 2]]);
    }
    std::vector<Simulation::NodePtr> vecNodes(uniqueNodes.begin(), uniqueNodes.end());
    numVertices = vecNodes.size();
    vertices = new Vertex[numVertices];

    // Reassign nodes indices
    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < numVertices; ++i)
    {
        auto node = vecNodes[i];
        node->index = i;
        vertices[i] = node->position;
    }

    // Get surface triangles
    numTriangles = IcosphereTrianglesIndicesSize;
    triangles = new Triangle[numTriangles];

    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < numTriangles; ++i)
    {
        triangles[i] = Triangle(nodes[IcosphereTrianglesIndices[i * 3 + 0]]->index,
                                nodes[IcosphereTrianglesIndices[i * 3 + 1]]->index,
                                nodes[IcosphereTrianglesIndices[i * 3 + 2]]->index);
    }
}

void SomaGeometry::_computeNeuritesNodes(Simulation::MeshPtr mesh,
                                         std::vector<Simulation::Nodes>& neuritesNodes,
                                         std::vector<Vector3f>& surfacePositions,
                                         const Samples& somaSamples)
{
    neuritesNodes.resize(somaSamples.size());
    surfacePositions.resize(somaSamples.size());

    for (uint32_t i = 0; i < somaSamples.size(); ++i)
    {
        auto sample = somaSamples[i];
        Vector3f startPos = sample->getPosition();
        float startRadius = sample->getRadius();

        // Compute icosphere surface position nearer to the neurite start
        Vector3f direction = (startPos - _somaCenter).normalized();
        Vector3f surfacePos = direction * _somaRadius + _somaCenter;
        surfacePositions[i] = surfacePos;

        // Check nodes within the neurite radius and add the nodes to neurite nodes
        for (auto node: mesh->nodes)
        {
            if ((node->position - surfacePos).abs() <= startRadius)
            {
                neuritesNodes[i].push_back(node);
                node->fixed = true;
            }
        }
    }
}

void SomaGeometry::_pullNeuritesNodes(std::vector<Simulation::Nodes>& neuritesNodes,
                                      std::vector<Vector3f>& surfacePositions,
                                      const Samples& somaSamples, float alpha)
{
    // Reposition nodes in the direction of the neurites starts
    for (uint32_t i = 0; i < somaSamples.size(); ++i)
    {
        Vector3f increment = (somaSamples[i]->getPosition() - surfacePositions[i]) * alpha;
        for (auto node: neuritesNodes[i])
            node->position = node->initPosition() + increment; 
    }
}

}  // namespace Ultraliser
