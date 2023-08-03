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

#include "NeuronSkeletonizer.h"
#include <data/meshes/simple/IcoSphere.h>

namespace Ultraliser
{

NeuronSkeletonizer::NeuronSkeletonizer(Volume* volume, const Mesh *mesh)
    : Skeletonizer(volume, mesh)
{

}

void NeuronSkeletonizer::skeletonizeVolume()
{
    LOG_TITLE("Skeletonization");

    // Start the timer
    TIMER_SET;

    applyVolumeThinning();
    constructGraph();
    segmentComponents();

    LOG_STATUS_IMPORTANT("Skeletonization Stats.");
    LOG_STATS(GET_TIME_SECONDS);
}

void NeuronSkeletonizer::skeletonizeVolumeBlockByBlock(const size_t& blockSize,
                                                 const size_t& numberOverlappingVoxels,
                                                 const size_t& numberZeroVoxels)
{
    thinVolumeBlockByBlock(blockSize, numberOverlappingVoxels, numberZeroVoxels);
    constructGraph();
    segmentComponents();
}

SkeletonNode* NeuronSkeletonizer::_addSomaNode()
{
    SkeletonNode* somaNode = new SkeletonNode();
    somaNode->index = _nodes.back()->index + 1;
    somaNode->isSoma = true;
    _nodes.push_back(somaNode);
    return somaNode;
}

void NeuronSkeletonizer::_segmentSomaMesh(SkeletonNode* somaNode)
{
    // The _somaMesh should be a complex geometry containing overlapping spheres that would define
    // its structure
    _somaMesh = new Mesh();

    // For every node in the skeleton, if the radius is greater than 2.0, then this is a candidate
    // for the soma sample
    size_t estimatedSomaRadius = 0;
    Vector3f estimatedSomaCenter;
    size_t numberSamples = 0;

    for (size_t i = 0; i < _branches.size(); ++i)
    {
        for (size_t j = 0; j < _branches[i]->nodes.size(); ++j)
        {
            auto& node0 = _branches[i]->nodes[j];
            if (node0->radius >= 2.0)
            {
                Mesh* sample = new IcoSphere(3);
                sample->scale(node0->radius, node0->radius, node0->radius);
                sample->translate(node0->point);
                sample->map(_shellPoints, false);
                _somaMesh->append(sample);
                sample->~Mesh();

                estimatedSomaRadius += node0->radius;
                estimatedSomaCenter += node0->point;
                numberSamples++;
            }
        }
    }

    estimatedSomaCenter /= numberSamples;
    estimatedSomaRadius /= numberSamples;

    // somaNode->radius = estimatedSomaRadius;
    // somaNode->point = estimatedSomaCenter;
}

void NeuronSkeletonizer::_removeBranchesInsideSoma(SkeletonNode* somaNode)
{
    size_t numRoots = 0;

    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];

        size_t countSamplesInsideSoma = 0;

        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            if (branch->nodes[j]->insideSoma)
            {
                countSamplesInsideSoma++;
            }
        }

        // If the count of the samples located inside the soma is zero, then it is a valid
        // branch, but it is not a root branch indeed
        if (countSamplesInsideSoma == 0)
        {
            branch->root = false;
            branch->valid = true;
        }

        // If all the branch nodes are located inside the soma, then it is not valid
        else if (countSamplesInsideSoma == branch->nodes.size())
        {
            branch->root = false;
            branch->valid = false;
        }

        // Otherwise, it is a branch that is connected to the soma
        else
        {
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

            branch->root = true;
            branch->valid = true;

            numRoots++;
        }
    }

    std::cout << "Roots: " << numRoots << std::endl;
}

void NeuronSkeletonizer::_segmentSomaVolume()
{
    _somaMesh->exportMesh("/data/neuron-meshes/output/soma", true);

    _volume->clear();
    _volume->surfaceVoxelization(_somaMesh, false, false);
    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);
    _volume->project("/data/neuron-meshes/output/soma", true);

    // TODO: This could be parallelized
    std::vector< size_t > somaVoxels;
    for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
    {
        if (_volume->isFilled(i))
        {
            somaVoxels.push_back(i);
        }
    }

    // Find out the nodes that are inside the soma
    size_t insideSoma = 0;
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];
        size_t key = _volume->mapTo1DIndexWithoutBoundCheck(
                    node->voxel.x(), node->voxel.y(), node->voxel.z());
        if (std::find(somaVoxels.begin(), somaVoxels.end(), key) != somaVoxels.end())
        {
            insideSoma++;

            // The soma is inside the soma
            node->insideSoma = true;
        }
    }

    // The volume is safe to be deallocated
    _volume->~Volume();
    _volume = nullptr;
}

void NeuronSkeletonizer::exportIndividualBranches(const std::string& prefix) const
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + TXT_EXTENSION;
    LOG_STATUS("Exporting Neuron Branches: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    LOOP_STARTS("Writing Branches");
    size_t progress = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // If the branch does not have any valid nodes, then don't write it
        // if (!_branches[i]->valid || _branches[i]->nodes.size() == 0) continue;
        // if (!_branches[i]->root) continue;

        if (_branches[i]->nodes.size() == 0) continue;


        LOOP_PROGRESS(progress, _branches.size());
        ++progress;

        // The @start marks a new branch in the file
        stream << "start\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }

        // The @end marks the terminal sample of a branch
        stream << "end\n";
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

void NeuronSkeletonizer::constructGraph()
{
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxels();

    // Assign accurate radii to the nodes of the graph, i.e. inflate the nodes
    _inflateNodes();

    // Connect the nodes to construct the edges of the graph
    _connectNodes(indicesMapper);

    // Remove the triangular configurations
    _removeTriangleLoops();

    // Add a virtual soma node, until the soma is reconstructed later
    auto somaNode = _addSomaNode();

    // Re-index the samples, for simplicity
    OMP_PARALLEL_FOR for (size_t i = 1; i <= _nodes.size(); ++i) { _nodes[i - 1]->index = i; }

    // Reconstruct the sections, or the branches from the nodes
    _buildBranchesFromNodes(_nodes);

    // Segmentthe soma mesh from the branches
    _segmentSomaMesh(somaNode);

    // Segment soma volume
    _segmentSomaVolume();

    // Validate the branches, and remove the branches inside the soma
    // _removeBranchesInsideSoma(somaNode);

    // Identify the possible roots

    // Segment the spines, if any

    // Display a report











    return;



    // The graph that will contain the nodes
    SkeletonNodes nodes;

//    // Every constructed node must have an identifier, or index.
//    size_t nodeIndex = 0;

//    // Map to associate between the indices of the voxels and the nodes in the graph
//    std::map< size_t, size_t > indicesMapper;

//    // Search the filled voxels in the volume
//    for (size_t i = 0; i < _volume->getWidth(); ++i)
//    {
//        for (size_t j = 0; j < _volume->getHeight(); ++j)
//        {
//            for (size_t k = 0; k < _volume->getDepth(); ++k)
//            {
//                // If the voxel is filled
//                if (_volume->isFilled(i, j, k))
//                {
//                    // Get the 1D index of the voxel
//                    size_t voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(i, j, k);

//                    // Get a point representing the center of the voxel (in the volume)
//                    Vector3f voxelPosition(i * 1.f, j * 1.f, k * 1.f);

//                    // Get a point in the same coordinate space of the mesh
//                    Vector3f nodePosition(voxelPosition);
//                    nodePosition -= _centerVolume;
//                    nodePosition.x() *= _scaleFactor.x();
//                    nodePosition.y() *= _scaleFactor.y();
//                    nodePosition.z() *= _scaleFactor.z();
//                    nodePosition += _centerMesh;

//                    // Add the node to the nodes list
//                    nodes.push_back(new SkeletonNode(voxelIndex, nodePosition, voxelPosition));

//                    // Mapper from voxel to node indices
//                    indicesMapper.insert(std::pair< size_t, size_t >(voxelIndex, nodeIndex));

//                    // New node
//                    nodeIndex++;
//                }
//            }
//        }
//    }

//    // Compute the approximate radii of all the nodes in the graph, based on the minimum distance
//    std::vector< float > nodesRadii;
//    nodesRadii.resize(nodes.size());

//    OMP_PARALLEL_FOR
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        float minimumDistance = std::numeric_limits< float >::max();
//        for (size_t j = 0; j < _shellPoints.size(); ++j)
//        {
//            const float distance = (nodes[i]->point - _shellPoints[j]).abs();

//            if (distance < minimumDistance)
//            {
//                minimumDistance = distance;
//            }
//        }
//        nodes[i]->radius = minimumDistance;
//        nodesRadii[i] = minimumDistance;
//    }

//    // Obtain the node with the largest radius, candidate for soma
//    const auto iterator = std::max_element(std::begin(nodesRadii), std::end(nodesRadii));
//    const auto& index = std::distance(std::begin(nodesRadii), iterator);
//    const auto& largestNode = nodes[index];

//    // Clear the auxiliary list
//    nodesRadii.clear();
//    nodesRadii.shrink_to_fit();

//    // Construct the graph and connect the nodes
//    // OMP_PARALLEL_FOR
//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        // Check if the node has been visited before
//        SkeletonNode* node = nodes[i];

//        // Count the number of the connected edges to the node
//        size_t connectedEdges = 0;

//        // Search for the neighbours
//        for (size_t l = 0; l < 26; l++)
//        {
//            size_t idx = node->voxel.x() + VDX[l];
//            size_t idy = node->voxel.y() + VDY[l];
//            size_t idz = node->voxel.z() + VDZ[l];

//            if (_volume->isFilled(idx, idy, idz))
//            {
//                // Connected edges
//                connectedEdges++;

//                // Find the index of the voxel
//                const auto& vIndex = _volume->mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

//                // Find the corresponding index of the node to access the node from the nodes list
//                const auto& nIndex = indicesMapper.find(vIndex)->second;

//                // Add the node to the edgeNodes, only to be able to access it later
//                node->edgeNodes.push_back(nodes[nIndex]);
//            }
//        }

//        if (connectedEdges == 1)
//            node->terminal = true;

//        if (connectedEdges > 2)
//            node->branching = true;

//    }

//    std::cout << "Trignale Nodes \n";

//    // Remove the triangular configurations

//    const size_t currentNodesSize = nodes.size();
//    for (size_t i = 0; i < currentNodesSize; ++i)
//    {
//        if (nodes[i]->branching)
//        {
//            SkeletonNodes sideNodes;
//            if (isTriangleNode(nodes[i], sideNodes))
//            {
//                if (nodes[i]->visited) continue;

//                auto& n1 = nodes[i];
//                auto& n2 = sideNodes[0];
//                auto& n3 = sideNodes[1];

//                fixTriangle(nodes, n1, n2, n3);

//                if (n1->edgeNodes.size() > 2)
//                    n1->branching = true;
//                else
//                    n1->branching = false;

//                if (n2->edgeNodes.size() > 2)
//                    n2->branching = true;
//                else
//                    n2->branching = false;

//                if (n3->edgeNodes.size() > 2)
//                    n3->branching = true;
//                else
//                    n3->branching = false;

//                n1->visited = true;
//                n2->visited = true;
//                n3->visited = true;
//            }
//        }
//    }

//    for (size_t i = 0; i < nodes.size(); ++i)
//    {
//        nodes[i]->visited = false;
//    }

//    std::cout << currentNodesSize << ", " << nodes.size() << "\n";
//    std::cout << "Trignale Nodes Done \n";


//    SkeletonNode* somaNode = new SkeletonNode();
//    somaNode->index = nodes.back()->index + 1;
//    somaNode->isSoma = true;
//    nodes.push_back(somaNode);

//     // Re-index the samples, for simplicity
//     // OMP_PARALLEL_FOR for (size_t i = 0; i < nodes.size(); ++i) { nodes[i]->index = i; }
}


//void Skeletonizer::segmentComponents(SkeletonNodes& nodes)
//{
//    SkeletonNode* somaNode = nodes.back();


//    // Build the branches from the nodes
//    SkeletonBranches branches = _buildBranchesFromNodes(nodes);

//    std::cout << "Branches: " << branches.size() << "\n";
//    std::cout << "Nodes (Samples): " << nodes.size() << "\n";

//    std::fstream stream;
//    stream.open("/abdellah2/scratch/thinning/output/projections/branches.txt", std::ios::out);
//    for (size_t i = 0; i < branches.size(); ++i)
//    {
//        stream << "start\n";

//        for (auto& node: branches[i]->nodes)
//        {
//            stream << node->point.x() << " "
//                   << node->point.y() << " "
//                   << node->point.z() << " "
//                   << node->radius << "\n";
//        }
//        stream << "end\n";
//    }
//    stream.close();

//    std::cout << "Branchs: " << branches.size() << "\n";

//    return;


















////    _somaMesh = _reconstructSoma(branches);

////    _volume->clear();
////    _volume->surfaceVoxelization(_somaMesh, false, false);
////    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

////    std::vector< size_t > somaVoxels;
////    for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
////    {
////        if (_volume->isFilled(i))
////        {
////            somaVoxels.push_back(i);
////        }
////    }

////    _volume->clear();

////    std::cout << "Soma Voxels " << somaVoxels.size() << "\n";


////    size_t insideSoma = 0;
////    // OMP_PARALLEL_FOR
////    for (size_t i = 0; i < nodes.size(); ++i)
////    {
////        auto& node = nodes[i];

////        size_t key = _volume->mapTo1DIndexWithoutBoundCheck(node->voxel.x(), node->voxel.y(), node->voxel.z());
////        if (std::find(somaVoxels.begin(), somaVoxels.end(), key) != somaVoxels.end())
////        {
////            insideSoma++;

////            // The soma is inside the soma
////            node->insideSoma = true;
////        }
////    }

////    std::cout << "Nodes Inside Soma " << insideSoma << "\n";


////    // OMP_PARALLEL_FOR
////    for (size_t i = 0; i < branches.size(); ++i)
////    {
////        auto& branch = branches[i];

////        size_t countSamplesInsideSoma = 0;

////        for (size_t j = 0; j < branch->nodes.size(); ++j)
////        {
////            if (branch->nodes[j]->insideSoma)
////            {
////                countSamplesInsideSoma++;
////            }
////        }

////        // If the count of the samples located inside the soma is zero, then it is a valid branch
////        if (countSamplesInsideSoma == 0)
////        {
////            branch->root = false;
////            branch->valid = true;
////        }

////        // If all the branch nodes are located inside the soma, then it is not valid
////        else if (countSamplesInsideSoma == branch->nodes.size())
////        {
////            branch->root = false;
////            branch->valid = false;
////        }

////        // Otherwise, it is a branch that is connected to the soma
////        else
////        {
////            std::cout << "I am a root << " << countSamplesInsideSoma << "\n";

////            SkeletonNodes newNodes;

////            // Get the first and last nodes
////            auto& firstNode = branch->nodes.front();
////            auto& lastNode = branch->nodes.back();

////            if (firstNode->insideSoma)
////            {
////                newNodes.push_back(somaNode);
////                for (size_t j = 0; j < branch->nodes.size(); ++j)
////                {
////                    if (!branch->nodes[j]->insideSoma)
////                    {
////                        newNodes.push_back(branch->nodes[j]);
////                    }
////                }
////            }
////            else if (lastNode->insideSoma)
////            {
////                for (size_t j = 0; j < branch->nodes.size(); ++j)
////                {
////                    if (!branch->nodes[j]->insideSoma)
////                    {
////                        newNodes.push_back(branch->nodes[j]);
////                    }
////                }
////                newNodes.push_back(somaNode);

////            }

////            branch->nodes.clear();
////            branch->nodes.shrink_to_fit();
////            branch->nodes = newNodes;

////            branch->root = true;
////            branch->valid = true;
////        }
////    }

////    std::cout << "Detecting Starting Points\n";

//////    // Get the starting points of the branches
//////    std::vector< Vector3f > rootsStartingPoints;

//////    SkeletonBranches possibleRoots;

//////    std::cout << branches.size() << " \n";
//////    for (size_t i = 0; i < branches.size(); ++i)
//////    {
//////        const auto& branch = branches[i];
//////        std::cout << i << " ";

//////        if (branch->root)
//////        {
//////            possibleRoots.push_back(branch);

//////            if (branch->nodes.size() ==0) continue;

//////            auto auxNode0 = branch->nodes.front();

//////            if (auxNode0->isSoma)
//////            {
//////                auto& rootStartingNode = branch->nodes[1];
//////                somaNode->edgeNodes.push_back(rootStartingNode);
//////                rootsStartingPoints.push_back(rootStartingNode->point);
//////            }
//////            else
//////            {
//////                auto& rootStartingNode = branch->nodes[branch->nodes.size() - 2];
//////                somaNode->edgeNodes.push_back(rootStartingNode);
//////                rootsStartingPoints.push_back(rootStartingNode->point);
//////            }
//////        }
//////    }

//////    std::cout << "Detecting Starting Points DONE\n";


//////    Vector3f center(0.f);
//////    float radius = 0;
//////    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
//////    {
//////        center += rootsStartingPoints[i];
//////    }

//////    center = center / rootsStartingPoints.size();
//////    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
//////    {
//////        radius += rootsStartingPoints[i].distance(center);
//////    }
//////    radius = radius / rootsStartingPoints.size();

//////    somaNode->point = center;
//////    somaNode->radius = radius;

//////    printf("Soma: %f, [%f, %f, %f], Roots: %ld\n",
//////           center.x(), center.y(), center.z(), radius, rootsStartingPoints.size());



////    stream.open("/abdellah2/scratch/thinning/output/projections/branches.txt", std::ios::out);
////    for (size_t i = 0; i < branches.size(); ++i)
////    {
////        stream << "start\n";

////        for (auto& node: branches[i]->nodes)
////        {
////            stream << node->point.x() << " "
////                   << node->point.y() << " "
////                   << node->point.z() << " "
////                   << node->radius << "\n";
////        }
////        stream << "end\n";
////    }
////    stream.close();

////    std::cout << "Branchs: " << branches.size() << "\n";



//////    _volume->clear();
//////    Mesh* sample = new IcoSphere(3);
//////    sample->scale(somaNode->radius, somaNode->radius, somaNode->radius);
//////    sample->translate(somaNode->point);
//////    _volume->surfaceVoxelization(sample);
//////    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

////    // Double loop to detect the connectivity of the branches
////    for (size_t i = 0; i < branches.size(); ++i)
////    {
////        auto& iBranch = branches[i];

////        if (!iBranch->valid)
////            continue;

////        const auto& iFirstNode = iBranch->nodes.front();
////        const auto& iLastNode = iBranch->nodes.back();

////        for (size_t j = 0; j < branches.size(); ++j)
////        {
////            const auto& jBranch = branches[j];

////            if (!jBranch->valid)
////               continue;

////            // Ignore the same branch
////            if (iBranch->index == jBranch->index)
////                continue;

////            const auto& jFirstNode = jBranch->nodes.front();
////            const auto& jLastNode = jBranch->nodes.back();

////            // Child
////            if (iLastNode->index == jFirstNode->index || iLastNode->index == jLastNode->index)
////            {
////                iBranch->children.push_back(jBranch);
////            }

////            // Parent
////            if (iFirstNode->index == jFirstNode->index || iFirstNode->index == jLastNode->index)
////            {
////                iBranch->parents.push_back(jBranch);
////            }
////        }
////    }

////    std::cout << "Heirarichy \n";

////     return;


//////    stream.open("/ssd3/scratch/skeletonization-tests/input/output/radiii.txt", std::ios::out);


//////    for (size_t i = 0; i < branches.size(); ++i)
//////    {


//////        if (!branches[i]->valid) continue;

//////        if (branches[i]->parents.size() == 0 && branches[i]->children.size() == 0)
//////        {
//////            for (auto& node: branches[i]->nodes)
//////            {
//////                stream << node->point.x() << " "
//////                       << node->point.y() << " "
//////                       << node->point.z() << " "
//////                       << node->radius << "\n";

//////                Mesh* sample = new IcoSphere(2);
//////                sample->scale(node->radius, node->radius, node->radius);
//////                sample->translate(node->point);
//////                _volume->surfaceVoxelization(sample);
//////                sample->~Mesh();
//////            }
//////       }
//////    }
//////    stream.close();
//////    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);




//////    // Re-arraning the branching
//////    for (size_t i = 0; i < possibleRoots.size(); ++i)
//////    {
//////        auto &nodes = possibleRoots[i]->nodes;

//////        if (nodes.size() == 0)
//////        {
//////            std::cout << "Issue \n";
//////        }
//////        else if (nodes.size() == 1)
//////        {
//////            std::cout << "Issue2 \n";
//////        }
//////        else if (nodes.size() > 1)
//////        {
//////            if (nodes[0]->isSoma)
//////            {
//////                // That's perfect
//////            }
//////            else
//////            {
//////                // Re-arrange the nodes in the branch
//////                std::reverse(nodes.begin(), nodes.end());
//////            }
//////        }
//////    }
//}



void Skeletonizer::segmentComponents()
{
    // SkeletonNode* somaNode = nodes.back();

    // Build the branches from the nodes
    _buildBranchesFromNodes(_nodes);

    std::cout << "Branches: " << _branches.size() << "\n";
    std::cout << "Nodes (Samples): " << _nodes.size() << "\n";

    std::fstream stream;
    stream.open("branches.txt", std::ios::out);
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        stream << "start\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        stream << "end\n";
    }
    stream.close();

    std::cout << "Branchs: " << _branches.size() << "\n";

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
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        stream << "start\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        stream << "end\n";
    }
    stream.close();

    std::cout << "Branchs: " << _branches.size() << "\n";



//    _volume->clear();
//    Mesh* sample = new IcoSphere(3);
//    sample->scale(somaNode->radius, somaNode->radius, somaNode->radius);
//    sample->translate(somaNode->point);
//    _volume->surfaceVoxelization(sample);
//    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

    // Double loop to detect the connectivity of the branches
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& iBranch = _branches[i];

        if (!iBranch->valid)
            continue;

        const auto& iFirstNode = iBranch->nodes.front();
        const auto& iLastNode = iBranch->nodes.back();

        for (size_t j = 0; j < _branches.size(); ++j)
        {
            const auto& jBranch = _branches[j];

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



}
