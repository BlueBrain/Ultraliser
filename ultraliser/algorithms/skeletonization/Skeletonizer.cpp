/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
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

#include "Skeletonizer.h"
#include <math/Vector.h>
#include "Neighbors.hh"
#include <data/meshes/simple/TriangleOperations.h>
#include <data/meshes/simple/IcoSphere.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(const Mesh *mesh, Volume* volume)
    : _mesh(mesh)
    , _volume(volume)
{
    // Mesh bounding box
    _pMinMesh = volume->getPMin();
    _pMaxMesh = volume->getPMax();
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    // Volume bounding box
    _pMinVolume = Vector3f(0.f);
    _pMaxVolume = Vector3f(volume->getWidth() * 1.f,
                           volume->getHeight() * 1.f,
                           volume->getDepth() * 1.f);
    _boundsVolume = _pMaxVolume;
    _centerVolume = 0.5 * _boundsVolume;

    // Mesh to volume scale factor
    _scaleFactor = _boundsMesh / _boundsVolume;

    _computeShellPoints();
}

#include<ctime>

void Skeletonizer::applyVolumeThinning()
{
    TIMER_SET;
    LOG_STATUS("Thinning Volume");


    // The thinning kernel that will be used to thin the volume

    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();
    LOG_STATUS("Thinning Starting");

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        size_t numberDeletedVoxels = _volume->deleteCandidateVoxels(thinningKernel);

        // Updating the progess bar
       if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
       LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                     initialNumberVoxelsToBeDeleted);

       if (numberDeletedVoxels == 0)
           break;

       loopCounter++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);








//    while(1)
//    {
//        size_t numberDeletedVoxels = 0;

////        // Search for the border voxels, for this iteration
////        std::cout << "\n1 ";
////        std::vector< std::vector< Vec3ui_64 > > borderVoxels = _volume->searchForBorderVoxels();
////        std::cout << "1 \n";

//        std::cout << "\n1 ";
//        auto candidateVoxels = _volume->searchForCandidateVoxelsOne();
//        std::cout << "1 \n";

//        for (size_t direction = 0; direction < 6; direction++)
//        {
//            // Search for the delerable voxels
////            std::vector< Vec3ui_64 > voxelsToBeDeleted =
////                    _volume->searchForDeletableVoxels(
////                        borderVoxels, thinningKernel, direction);

//            _volume->confirmDeletableVoxels(candidateVoxels, thinningKernel, direction);

//            // Delete the voxels
////            for (size_t i = 0; i < voxelsToBeDeleted.size(); ++i)
////            {
////                numberDeletedVoxels++;
////                _volume->clear(voxelsToBeDeleted[i].x(),
////                               voxelsToBeDeleted[i].y(),
////                               voxelsToBeDeleted[i].z());
////            }

//            for (size_t i = 0; i < candidateVoxels.size(); ++i)
//            {
//                if (candidateVoxels[i]->deletable)
//                {
//                    numberDeletedVoxels++;
//                    _volume->clear(candidateVoxels[i]->x,
//                                   candidateVoxels[i]->y,
//                                   candidateVoxels[i]->z);

//                    candidateVoxels[i]->deletable = false;
//                }
//            }
//        }

//         // Clear the border voxels (list of lists)
////        for (size_t i = 0; i < borderVoxels.size(); ++i)
////        {
////            borderVoxels[i].clear();
////            borderVoxels[i].shrink_to_fit();
////        }

//        candidateVoxels.clear();
//        candidateVoxels.shrink_to_fit();

//         // Updating the progess bar
//        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
//        LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
//                      initialNumberVoxelsToBeDeleted);

//        if (numberDeletedVoxels == 0)
//            break;

//        loopCounter++;
//    }
//    LOOP_DONE;
//    LOG_STATS(GET_TIME_SECONDS);
}

std::vector< Vector3f > Skeletonizer::getShellPoints()
{
    return _shellPoints;
}

void Skeletonizer::_computeShellPoints()
{
    // Search for the border voxels (the shell voxels) of the volume
    std::vector< std::vector< Vec3ui_64 > > perSliceSurfaceShell = _volume->searchForBorderVoxels();

    // Concatinate the points in a single list
    for (size_t i = 0; i < perSliceSurfaceShell.size(); ++i)
    {
        for (size_t j = 0; j < perSliceSurfaceShell[i].size(); ++j)
        {
            const auto voxel = perSliceSurfaceShell[i][j];
            _shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSliceSurfaceShell[i].clear();
        perSliceSurfaceShell[i].shrink_to_fit();
    }
    perSliceSurfaceShell.clear();
    perSliceSurfaceShell.shrink_to_fit();

    // Adjust the locations of the shell points taking into consideration the mesh coordinates
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _shellPoints.size(); ++i)
    {
        // Center the shell points (of the volume) at the origin
        _shellPoints[i] -= _centerVolume;

        // Scale to match the dimensions of the mesh
        _shellPoints[i].x() *= _scaleFactor.x();
        _shellPoints[i].y() *= _scaleFactor.y();
        _shellPoints[i].z() *= _scaleFactor.z();

        // Translate to the center of the mesh
        _shellPoints[i] += _centerMesh;
    }
}


void removeEdgeNode(SkeletonNode* node, SkeletonNode* edgeNodeMarkedForRemoval)
{
    auto it = std::find(node->edgeNodes.begin(), node->edgeNodes.end(), edgeNodeMarkedForRemoval);

    if(it != node->edgeNodes.end())
        node->edgeNodes.erase(it);
    node->edgeNodes.shrink_to_fit();
}

void fixTriangle(SkeletonNodes& nodes,
                 SkeletonNode* n1, SkeletonNode* n2, SkeletonNode* n3)
{
    SkeletonNode* centerNode = new SkeletonNode();
    centerNode->point = (n1->point + n2->point + n3->point) / 3.f;
    centerNode->radius = (n1->radius + n2->radius + n3->radius) / 3.f;

    centerNode->index = nodes.back()->index + 1;
    centerNode->branching = true;
    centerNode->terminal = false;

    removeEdgeNode(n1, n2); removeEdgeNode(n1, n3);
    removeEdgeNode(n2, n1); removeEdgeNode(n2, n3);
    removeEdgeNode(n3, n1); removeEdgeNode(n3, n2);

    centerNode->edgeNodes.push_back(n1);
    centerNode->edgeNodes.push_back(n2);
    centerNode->edgeNodes.push_back(n3);

    n1->edgeNodes.push_back(centerNode);
    n2->edgeNodes.push_back(centerNode);
    n3->edgeNodes.push_back(centerNode);

    nodes.push_back(centerNode);
}

bool areConnected(const SkeletonNode* n1, const SkeletonNode* n2)
{
    for (size_t i = 0; i < n1->edgeNodes.size(); ++i)
    {
        if (n1->edgeNodes[i]->index == n2->index)
        {
            return true;
        }
    }

    return false;
}

bool isTriangleNode(const SkeletonNode* n, SkeletonNodes& connectedEdgeNodes)
{
    for (size_t i = 0; i < n->edgeNodes.size(); ++i)
    {
        for (size_t j = 0; j < n->edgeNodes.size(); ++j)
        {
            if (i == j) continue;
            if (n->edgeNodes[i]->index == n->index) continue;
            if (n->edgeNodes[j]->index == n->index) continue;

            if (areConnected(n->edgeNodes[i], n->edgeNodes[j]))
            {
                connectedEdgeNodes.push_back(n->edgeNodes[j]);
                connectedEdgeNodes.push_back(n->edgeNodes[i]);
                return true;
            }
        }
    }

    return false;
}

SkeletonNodes Skeletonizer::constructGraph()
{
    // The graph that will contain the nodes
    SkeletonNodes nodes;

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;

    // Search the filled voxels in the volume
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        for (size_t j = 0; j < _volume->getHeight(); ++j)
        {
            for (size_t k = 0; k < _volume->getDepth(); ++k)
            {
                // If the voxel is filled
                if (_volume->isFilled(i, j, k))
                {
                    // Get the 1D index of the voxel
                    size_t voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(i, j, k);

                    // Get a point representing the center of the voxel (in the volume)
                    Vector3f voxelPosition(i * 1.f, j * 1.f, k * 1.f);

                    // Get a point in the same coordinate space of the mesh
                    Vector3f nodePosition(voxelPosition);
                    nodePosition -= _centerVolume;
                    nodePosition.x() *= _scaleFactor.x();
                    nodePosition.y() *= _scaleFactor.y();
                    nodePosition.z() *= _scaleFactor.z();
                    nodePosition += _centerMesh;

                    // Add the node to the nodes list
                    nodes.push_back(new SkeletonNode(voxelIndex, nodePosition, voxelPosition));

                    // Mapper from voxel to node indices
                    indicesMapper.insert(std::pair< size_t, size_t >(voxelIndex, nodeIndex));

                    // New node
                    nodeIndex++;
                }
            }
        }
    }

    // Compute the approximate radii of all the nodes in the graph, based on the minimum distance
    std::vector< float > nodesRadii;
    nodesRadii.resize(nodes.size());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        float minimumDistance = std::numeric_limits< float >::max();
        for (size_t j = 0; j < _shellPoints.size(); ++j)
        {
            const float distance = (nodes[i]->point - _shellPoints[j]).abs();

            if (distance < minimumDistance)
            {
                minimumDistance = distance;
            }
        }
        nodes[i]->radius = minimumDistance;
        nodesRadii[i] = minimumDistance;
    }

    // Obtain the node with the largest radius, candidate for soma
    const auto iterator = std::max_element(std::begin(nodesRadii), std::end(nodesRadii));
    const auto& index = std::distance(std::begin(nodesRadii), iterator);
    const auto& largestNode = nodes[index];

    // Clear the auxiliary list
    nodesRadii.clear();
    nodesRadii.shrink_to_fit();

    // Construct the graph and connect the nodes
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        // Check if the node has been visited before
        SkeletonNode* node = nodes[i];

        // Count the number of the connected edges to the node
        size_t connectedEdges = 0;

        // Search for the neighbours
        for (size_t l = 0; l < 26; l++)
        {
            size_t idx = node->voxel.x() + VDX[l];
            size_t idy = node->voxel.y() + VDY[l];
            size_t idz = node->voxel.z() + VDZ[l];

            if (_volume->isFilled(idx, idy, idz))
            {
                // Connected edges
                connectedEdges++;

                // Find the index of the voxel
                const auto& vIndex = _volume->mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

                // Find the corresponding index of the node to access the node from the nodes list
                const auto& nIndex = indicesMapper.find(vIndex)->second;

                // Add the node to the edgeNodes, only to be able to access it later
                node->edgeNodes.push_back(nodes[nIndex]);
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;

        if (connectedEdges > 2)
            node->branching = true;

    }

    std::cout << "Trignale Nodes \n";

    // Remove the triangular configurations

    const size_t currentNodesSize = nodes.size();
    for (size_t i = 0; i < currentNodesSize; ++i)
    {
        if (nodes[i]->branching)
        {
            SkeletonNodes sideNodes;
            if (isTriangleNode(nodes[i], sideNodes))
            {
                if (nodes[i]->visited) continue;

                auto& n1 = nodes[i];
                auto& n2 = sideNodes[0];
                auto& n3 = sideNodes[1];

                fixTriangle(nodes, n1, n2, n3);

                if (n1->edgeNodes.size() > 2)
                    n1->branching = true;
                else
                    n1->branching = false;

                if (n2->edgeNodes.size() > 2)
                    n2->branching = true;
                else
                    n2->branching = false;

                if (n3->edgeNodes.size() > 2)
                    n3->branching = true;
                else
                    n3->branching = false;

                n1->visited = true;
                n2->visited = true;
                n3->visited = true;
            }
        }
    }

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        nodes[i]->visited = false;
    }

    std::cout << currentNodesSize << ", " << nodes.size() << "\n";
    std::cout << "Trignale Nodes Done \n";


//    SkeletonNode* somaNode = new SkeletonNode();
//    somaNode->index = nodes.back()->index + 1;
//    somaNode->isSoma = true;
//    nodes.push_back(somaNode);

     // Re-index the samples, for simplicity
     // OMP_PARALLEL_FOR for (size_t i = 0; i < nodes.size(); ++i) { nodes[i]->index = i; }



    return nodes;
}

void Skeletonizer::segmentComponents(SkeletonNodes& nodes)
{
    // SkeletonNode* somaNode = nodes.back();


    // Build the branches from the nodes
    SkeletonBranches branches = _buildBranchesFromNodes(nodes);

    std::cout << "Branches: " << branches.size() << "\n";
    std::cout << "Nodes (Samples): " << nodes.size() << "\n";

    std::fstream stream;
    stream.open("branches.txt", std::ios::out);
    for (size_t i = 0; i < branches.size(); ++i)
    {
        stream << "start\n";

        for (auto& node: branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        stream << "end\n";
    }
    stream.close();

    std::cout << "Branchs: " << branches.size() << "\n";

    return;


















//    _somaMesh = _reconstructSoma(branches);

//    _volume->clear();
//    _volume->surfaceVoxelization(_somaMesh, false, false);
//    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

//    std::vector< size_t > somaVoxels;
//    for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
//    {
//        if (_volume->isFilled(i))
//        {
//            somaVoxels.push_back(i);
//        }
//    }

//    _volume->clear();

//    std::cout << "Soma Voxels " << somaVoxels.size() << "\n";


//    size_t insideSoma = 0;
//    // OMP_PARALLEL_FOR
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        auto& node = nodes[i];

//        size_t key = _volume->mapTo1DIndexWithoutBoundCheck(node->voxel.x(), node->voxel.y(), node->voxel.z());
//        if (std::find(somaVoxels.begin(), somaVoxels.end(), key) != somaVoxels.end())
//        {
//            insideSoma++;

//            // The soma is inside the soma
//            node->insideSoma = true;
//        }
//    }

//    std::cout << "Nodes Inside Soma " << insideSoma << "\n";


//    // OMP_PARALLEL_FOR
//    for (size_t i = 0; i < branches.size(); ++i)
//    {
//        auto& branch = branches[i];

//        size_t countSamplesInsideSoma = 0;

//        for (size_t j = 0; j < branch->nodes.size(); ++j)
//        {
//            if (branch->nodes[j]->insideSoma)
//            {
//                countSamplesInsideSoma++;
//            }
//        }

//        // If the count of the samples located inside the soma is zero, then it is a valid branch
//        if (countSamplesInsideSoma == 0)
//        {
//            branch->root = false;
//            branch->valid = true;
//        }

//        // If all the branch nodes are located inside the soma, then it is not valid
//        else if (countSamplesInsideSoma == branch->nodes.size())
//        {
//            branch->root = false;
//            branch->valid = false;
//        }

//        // Otherwise, it is a branch that is connected to the soma
//        else
//        {
//            std::cout << "I am a root << " << countSamplesInsideSoma << "\n";

//            SkeletonNodes newNodes;

//            // Get the first and last nodes
//            auto& firstNode = branch->nodes.front();
//            auto& lastNode = branch->nodes.back();

//            if (firstNode->insideSoma)
//            {
//                newNodes.push_back(somaNode);
//                for (size_t j = 0; j < branch->nodes.size(); ++j)
//                {
//                    if (!branch->nodes[j]->insideSoma)
//                    {
//                        newNodes.push_back(branch->nodes[j]);
//                    }
//                }
//            }
//            else if (lastNode->insideSoma)
//            {
//                for (size_t j = 0; j < branch->nodes.size(); ++j)
//                {
//                    if (!branch->nodes[j]->insideSoma)
//                    {
//                        newNodes.push_back(branch->nodes[j]);
//                    }
//                }
//                newNodes.push_back(somaNode);

//            }

//            branch->nodes.clear();
//            branch->nodes.shrink_to_fit();
//            branch->nodes = newNodes;

//            branch->root = true;
//            branch->valid = true;
//        }
//    }

//    std::cout << "Detecting Starting Points\n";

//    // Get the starting points of the branches
//    std::vector< Vector3f > rootsStartingPoints;

//    SkeletonBranches possibleRoots;

//    std::cout << branches.size() << " \n";
//    for (size_t i = 0; i < branches.size(); ++i)
//    {
//        const auto& branch = branches[i];
//        std::cout << i << " ";

//        if (branch->root)
//        {
//            possibleRoots.push_back(branch);

//            if (branch->nodes.size() ==0) continue;

//            auto auxNode0 = branch->nodes.front();

//            if (auxNode0->isSoma)
//            {
//                auto& rootStartingNode = branch->nodes[1];
//                somaNode->edgeNodes.push_back(rootStartingNode);
//                rootsStartingPoints.push_back(rootStartingNode->point);
//            }
//            else
//            {
//                auto& rootStartingNode = branch->nodes[branch->nodes.size() - 2];
//                somaNode->edgeNodes.push_back(rootStartingNode);
//                rootsStartingPoints.push_back(rootStartingNode->point);
//            }
//        }
//    }

//    std::cout << "Detecting Starting Points DONE\n";


//    Vector3f center(0.f);
//    float radius = 0;
//    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
//    {
//        center += rootsStartingPoints[i];
//    }

//    center = center / rootsStartingPoints.size();
//    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
//    {
//        radius += rootsStartingPoints[i].distance(center);
//    }
//    radius = radius / rootsStartingPoints.size();

//    somaNode->point = center;
//    somaNode->radius = radius;

//    printf("Soma: %f, [%f, %f, %f], Roots: %ld\n",
//           center.x(), center.y(), center.z(), radius, rootsStartingPoints.size());



    stream.open("branches.txt", std::ios::out);
    for (size_t i = 0; i < branches.size(); ++i)
    {
        stream << "start\n";

        for (auto& node: branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        stream << "end\n";
    }
    stream.close();

    std::cout << "Branchs: " << branches.size() << "\n";



//    _volume->clear();
//    Mesh* sample = new IcoSphere(3);
//    sample->scale(somaNode->radius, somaNode->radius, somaNode->radius);
//    sample->translate(somaNode->point);
//    _volume->surfaceVoxelization(sample);
//    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

    // Double loop to detect the connectivity of the branches
    for (size_t i = 0; i < branches.size(); ++i)
    {
        auto& iBranch = branches[i];

        if (!iBranch->valid)
            continue;

        const auto& iFirstNode = iBranch->nodes.front();
        const auto& iLastNode = iBranch->nodes.back();

        for (size_t j = 0; j < branches.size(); ++j)
        {
            const auto& jBranch = branches[j];

            if (!jBranch->valid)
               continue;

            // Ignore the same branch
            if (iBranch->index == jBranch->index)
                continue;

            const auto& jFirstNode = jBranch->nodes.front();
            const auto& jLastNode = jBranch->nodes.back();

            // Child
            if (iLastNode->index == jFirstNode->index || iLastNode->index == jLastNode->index)
            {
                iBranch->children.push_back(jBranch);
            }

            // Parent
            if (iFirstNode->index == jFirstNode->index || iFirstNode->index == jLastNode->index)
            {
                iBranch->parents.push_back(jBranch);
            }
        }
    }

    std::cout << "Heirarichy \n";

     return;


//    stream.open("/ssd3/scratch/skeletonization-tests/input/output/radiii.txt", std::ios::out);


//    for (size_t i = 0; i < branches.size(); ++i)
//    {


//        if (!branches[i]->valid) continue;

//        if (branches[i]->parents.size() == 0 && branches[i]->children.size() == 0)
//        {
//            for (auto& node: branches[i]->nodes)
//            {
//                stream << node->point.x() << " "
//                       << node->point.y() << " "
//                       << node->point.z() << " "
//                       << node->radius << "\n";

//                Mesh* sample = new IcoSphere(2);
//                sample->scale(node->radius, node->radius, node->radius);
//                sample->translate(node->point);
//                _volume->surfaceVoxelization(sample);
//                sample->~Mesh();
//            }
//       }
//    }
//    stream.close();
//    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);




//    // Re-arraning the branching
//    for (size_t i = 0; i < possibleRoots.size(); ++i)
//    {
//        auto &nodes = possibleRoots[i]->nodes;

//        if (nodes.size() == 0)
//        {
//            std::cout << "Issue \n";
//        }
//        else if (nodes.size() == 1)
//        {
//            std::cout << "Issue2 \n";
//        }
//        else if (nodes.size() > 1)
//        {
//            if (nodes[0]->isSoma)
//            {
//                // That's perfect
//            }
//            else
//            {
//                // Re-arrange the nodes in the branch
//                std::reverse(nodes.begin(), nodes.end());
//            }
//        }
//    }
}



Mesh* Skeletonizer::_reconstructSoma(const SkeletonBranches& branches)
{
    Mesh* somaMesh = new Mesh();

    for (size_t i = 0; i < branches.size(); ++i)
    {
        for (size_t j = 0; j < branches[i]->nodes.size(); ++j)
        {
            auto& node0 = branches[i]->nodes[j];
            if (node0->radius >= 2.0)
            {
                Mesh* sample = new IcoSphere(3);
                sample->scale(node0->radius, node0->radius, node0->radius);
                sample->translate(node0->point);

                sample->map(_shellPoints);

                somaMesh->append(sample);
                sample->~Mesh();
            }
        }
    }

    return somaMesh;
}


SkeletonBranches Skeletonizer::_buildBranchesFromNodes(const SkeletonNodes& nodes)
{
    // A list to collect all the constructed branches from the skeleton graph
    SkeletonBranches branches;

    // Used to index the branch
    size_t branchIndex = 0;

    // Construct the hierarchy to the terminal
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto& node = nodes[i];

        // The node must be branching
        if (node->branching)
        {
            // The node must be visited less number of times than its branching edges
            if (node->iVisit < node->edgeNodes.size())
            {
                // Construct the branch, starting with the edge node
                for (size_t j = 0; j < node->edgeNodes.size(); ++j)
                {
                    // Get a reference to the edge node
                    auto& edgeNode = node->edgeNodes[j];

                    if (edgeNode->iVisit >= edgeNode->edgeNodes.size()) continue;

                    // If the edge node is a terminal
                    if (edgeNode->terminal)
                    {
                        SkeletonBranch* branch = new SkeletonBranch();

                        node->iVisit += 1;
                        branch->nodes.push_back(node);

                        edgeNode->iVisit += 1;
                        branch->nodes.push_back(edgeNode);

                        branch->index = branchIndex;
                        branchIndex++;
                        branches.push_back(branch);
                    }

                    // If the edge node is a branching node
                    else if (edgeNode->branching)
                    {
                        SkeletonBranch* branch = new SkeletonBranch();

                        node->iVisit += 1;
                        branch->nodes.push_back(node);

                        edgeNode->iVisit += 1;
                        branch->nodes.push_back(edgeNode);

                        branch->index = branchIndex;
                        branchIndex++;
                        branches.push_back(branch);
                    }

                    // If the edge node is an intermediate node
                    else
                    {
                        // Ensure that the edge node is not visited before to make a branch
                        if (edgeNode->iVisit < 1)
                        {
                            SkeletonBranch* branch = new SkeletonBranch();

                            node->iVisit += 1;
                            branch->nodes.push_back(node);

                            edgeNode->iVisit += 1;
                            branch->nodes.push_back(edgeNode);

                            // The previous node is the first node
                            SkeletonNode *previousNode = node;

                            // The current node is the edge node
                            SkeletonNode *currentNode = edgeNode;

                            // Ensure that the current node has only two connected edges (or nodes)
                            while (true)
                            {
                                // Get a reference to the connecting nodes to the current node
                                auto edgeNode0 = currentNode->edgeNodes[0];
                                auto edgeNode1 = currentNode->edgeNodes[1];

                                // Ignore the previous node
                                if (edgeNode0->index == previousNode->index)
                                {
                                    previousNode = currentNode;
                                    currentNode = edgeNode1;
                                }
                                else
                                {
                                    previousNode = currentNode;
                                    currentNode = edgeNode0;
                                }

                                currentNode->iVisit += 1;
                                branch->nodes.push_back(currentNode);

                                if (!(currentNode->edgeNodes.size() == 2))
                                    break;
                            }

                            branch->index = branchIndex;
                            branchIndex++;
                            branches.push_back(branch);
                        }

                    }
                }
            }
        }
    }

    return branches;
}

}
