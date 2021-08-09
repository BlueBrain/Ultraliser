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
SomaGeometry::SomaGeometry(NeuronMorphology* morphology)
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
        _fixCenterNodes(mesh->nodes, 0.5);
        _generateIcosphereSprings(mesh, 500);
        _generatePullSprings(mesh, 100.0f, morphology->getFirstSections(), 0.2);

        Simulation::AnimSystem animSystem(0.01f);
        animSystem.animate(mesh, 1000);
        _nodesToVertices(mesh->nodes);
    }
    else
    {
        LOG_ERROR("The input morphology is EMPTY!");
    }
}

void SomaGeometry::_insert(Simulation::SpringPtr spring, Simulation::UniqueSprings& springs)
{
    std::pair<Simulation::UniqueSprings::iterator, bool> insertion = springs.insert(spring);

    if (!insertion.second)
        delete spring;
}

Simulation::MeshPtr SomaGeometry::_loadIcosphereGeometry()
{
    auto mesh = new Simulation::Mesh(1000.f, 0.2f);
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

void SomaGeometry::_generateIcosphereSprings(Simulation::MeshPtr mesh, float stiffness)
{
    Simulation::UniqueSprings uSprings;
    for (uint64_t i = 0; i < IcosphereTetsIndicesSize; ++i)
    {
        Simulation::NodePtr node0 = mesh->nodes[IcosphereTetsIndices[i * 4 + 0]];
        Simulation::NodePtr node1 = mesh->nodes[IcosphereTetsIndices[i * 4 + 1]];
        Simulation::NodePtr node2 = mesh->nodes[IcosphereTetsIndices[i * 4 + 2]];
        Simulation::NodePtr node3 = mesh->nodes[IcosphereTetsIndices[i * 4 + 3]];

        _insert(new Simulation::Spring(node0, node1, stiffness), uSprings);
        _insert(new Simulation::Spring(node0, node2, stiffness), uSprings);
        _insert(new Simulation::Spring(node0, node3, stiffness), uSprings);
        _insert(new Simulation::Spring(node1, node2, stiffness), uSprings);
        _insert(new Simulation::Spring(node1, node3, stiffness), uSprings);
        _insert(new Simulation::Spring(node2, node3, stiffness), uSprings);
    }
    mesh->springs.insert(mesh->springs.end(), uSprings.begin(), uSprings.end());
}

void SomaGeometry::_generatePullSprings(Simulation::MeshPtr mesh,
                                        float stiffness,
                                        const Sections& firstSections,
                                        float restLengthThreshold)
{
    for (auto section : firstSections)
    {
        Vector3f startPos = section->getSamples()[0]->getPosition();
        float startRadius = section->getSamples()[0]->getRadius();

        // Compute icosphere surface position nearer to the neurite start
        Vector3f direction = (startPos - _somaCenter).normalized();
        Vector3f surfacePos = direction * _somaRadius + _somaCenter;

        // Compute displacement from icosphere surface to neurite start
        float diff = (startPos - surfacePos).abs();
        if (diff < 0.1)
        {
            diff = 0.1;
        }
        direction *= diff;

        // Check nodes within the neurite radius and create a spring connecting
        // them to the neurite start
        size_t nodesSize = mesh->nodes.size();
        for (uint64_t i = 0; i < nodesSize; ++i)
        {
            auto node = mesh->nodes[i];
            if ((node->position - surfacePos).abs() < startRadius)
            {
                Vector3f linkPosition = node->position + direction;
                auto linkNode = new Simulation::Node(linkPosition, true);
                mesh->nodes.push_back(linkNode);
                mesh->springs.push_back(new Simulation::Spring(
                    node, linkNode, stiffness, diff * restLengthThreshold));
            }
        }
    }
}

void SomaGeometry::_fixCenterNodes(Simulation::Nodes& nodes, float threshold)
{
    for (auto node : nodes)
    {
        if ((node->position - _somaCenter).abs() < _somaRadius * threshold)
        {
            node->fixed = true;
        }
    }
}

}  // namespace Ultraliser
