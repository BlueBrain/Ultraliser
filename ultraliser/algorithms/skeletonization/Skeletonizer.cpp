#include "Skeletonizer.h"
#include <math/Vector.h>
#include "Neighbors.hh"
#include <data/meshes/simple/TriangleOperations.h>
#include <data/meshes/simple/IcoSphere.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(const Mesh *mesh,
                           Volume* volume)
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

void Skeletonizer::applyVolumeThinning()
{
    TIMER_SET;
    LOG_STATUS("Thinning Volume");

    // A block of the volume that is scanned every iteration
    int8_t volumeBlock[26];

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        size_t numberDeletedVoxels = 0;

        // Search for the border voxels, for this iteration
        std::vector< std::vector< Vec3ui_64 > > borderVoxels = _volume->searchForBorderVoxels();

        for (size_t direction = 0; direction < 6; direction++)
        {
            // Search for the delerable voxels
            std::vector< Vec3ui_64 > voxelsToBeDeleted =
                    _volume->searchForDeletableVoxels(
                        borderVoxels, thinningKernel, direction, volumeBlock);

            // Delete the voxels
            for (size_t i = 0; i < voxelsToBeDeleted.size(); ++i)
            {
                numberDeletedVoxels++;
                _volume->clear(voxelsToBeDeleted[i].x(),
                               voxelsToBeDeleted[i].y(),
                               voxelsToBeDeleted[i].z());
            }

            // Clear the container of the deleted voxels
            voxelsToBeDeleted.clear();
        }

        // Clear the border voxels (list of lists)
        for (size_t i = 0; i < borderVoxels.size(); ++i)
        {
            borderVoxels[i].clear();
            borderVoxels[i].shrink_to_fit();
        }
        borderVoxels.clear();
        borderVoxels.shrink_to_fit();

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

void Skeletonizer::constructGraph()
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
    OMP_PARALLEL_FOR
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
                connectedEdges++;

                // Find the index of the voxel
                const auto& voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

                // Find the corresponding index of the node to access the node from the nodes list
                const auto& nodeIndex = indicesMapper.find(voxelIndex)->second;

                // Add the node to the edgeNodes, only to be able to access it later
                node->edgeNodes.push_back(nodes[nodeIndex]);
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;
        else if (connectedEdges > 2)
            node->branching = true;

    }

    SkeletonNode* somaNode = new SkeletonNode();
    somaNode->index = nodeIndex++;
    somaNode->isSoma = true;
    nodes.push_back(somaNode);

    // Re-index the samples, for simplicity
    OMP_PARALLEL_FOR for (size_t i = 0; i < nodes.size(); ++i) { nodes[i]->index = i; }

    // Build the branches from the nodes
    SkeletonBranches branches = _buildBranchesFromNodes(nodes);

    std::cout << "Branches: " << branches.size() << "\n";
    std::cout << "Nodes (Samples): " << nodes.size() << "\n";

    _somaMesh = _reconstructSoma(branches);

//    Vector3f pMinSoma, pMaxSoma;
//    _somaMesh->computeBoundingBox(pMinSoma, pMaxSoma);

    _volume->clear();
    _volume->surfaceVoxelization(_somaMesh, false, false);
    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

    std::vector< size_t > somaVoxels;
    for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
    {
        if (_volume->isFilled(i))
        {
            somaVoxels.push_back(i);
        }
    }

    _volume->clear();

    std::cout << "Soma Voxels " << somaVoxels.size() << "\n";


    size_t insideSoma = 0;
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto& node = nodes[i];

        size_t key = _volume->mapTo1DIndexWithoutBoundCheck(node->voxel.x(), node->voxel.y(), node->voxel.z());
        if (std::find(somaVoxels.begin(), somaVoxels.end(), key) != somaVoxels.end())
        // if (std::count(somaVoxels.begin(), somaVoxels.end(), key))
        {
            insideSoma++;

            // The soma is inside the soma
            node->insideSoma = true;
        }
    }

    std::cout << "Nodes Inside Soma " << insideSoma << "\n";


    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < branches.size(); ++i)
    {
        auto& branch = branches[i];

        size_t countSamplesInsideSoma = 0;

        for (size_t j = 0; j < branch->nodes.size(); ++j)
        {
            if (branch->nodes[j]->insideSoma)
            {
                countSamplesInsideSoma++;
            }
        }

        // If the count of the samples located inside the soma is zero, then it is a valid branch
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

            std::cout << "I am a root << " << countSamplesInsideSoma << "\n";

            SkeletonNodes newNodes;

            // Get the first and last nodes
            auto& firstNode = branch->nodes.front();
            auto& lastNode = branch->nodes.back();

            if (firstNode->insideSoma)
            {
                newNodes.push_back(somaNode);
                for (size_t j = 0; j < branch->nodes.size(); ++j)
                {
                    if (!branch->nodes[j]->insideSoma)
                    {
                        newNodes.push_back(branch->nodes[j]);
                    }
                }
            }
            else if (lastNode->insideSoma)
            {
                for (size_t j = 0; j < branch->nodes.size(); ++j)
                {
                    if (!branch->nodes[j]->insideSoma)
                    {
                        newNodes.push_back(branch->nodes[j]);
                    }
                }
                newNodes.push_back(somaNode);

            }

            branch->nodes.clear();
            branch->nodes.shrink_to_fit();
            branch->nodes = newNodes;

            branch->root = true;
            branch->valid = true;
        }

        // std::cout << i << "/" << branches.size() << "\n";
    }

    std::cout << "Detecting Starting Points\n";

    // Get the starting points of the branches
    std::vector< Vector3f > rootsStartingPoints;

    for (size_t i = 0; i < branches.size(); ++i)
    {
        const auto& branch = branches[i];

        if (branch->root)
        {
            std::cout << "Root " << branches[i]->index <<"\n";
            std::cout << branch->nodes.size() << "\n";

            auto auxNode0 = branch->nodes.front();

            if (auxNode0->isSoma)
            {
                auto& rootStartingNode = branch->nodes[1];
                somaNode->edgeNodes.push_back(rootStartingNode);
                rootsStartingPoints.push_back(rootStartingNode->point);
            }
            else
            {
                auto& rootStartingNode = branch->nodes[branch->nodes.size() - 2];
                somaNode->edgeNodes.push_back(rootStartingNode);
                rootsStartingPoints.push_back(rootStartingNode->point);
            }
        }
    }

    std::cout << "Roots Starting Points: " << rootsStartingPoints.size() << "\n";
    Vector3f center(0.f);
    float radius = 0;
    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
    {
        center += rootsStartingPoints[i];
    }

    center = center / rootsStartingPoints.size();
    for (size_t i = 0; i < rootsStartingPoints.size(); i++)
    {
        radius += rootsStartingPoints[i].distance(center);
    }
    radius = radius / rootsStartingPoints.size();

    somaNode->point = center;
    somaNode->radius = radius;

    printf("Soma: %f, [%f, %f, %f], Roots: %ld\n",
           center.x(), center.y(), center.z(), radius, rootsStartingPoints.size());



    std::cout << "Sphere \n";
    _volume->clear();
    Mesh* sample = new IcoSphere(3);
    sample->scale(somaNode->radius, somaNode->radius, somaNode->radius);
    sample->translate(somaNode->point);
    _volume->surfaceVoxelization(sample);
    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);
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
    SkeletonBranches branches;
    size_t branchIndex = 0;

    // Construct the heirarichy to the terminal
    for (size_t i = 0; i < nodes.size(); ++i)
    {
        auto& node = nodes[i];

        // The node must be branching
        if (node->branching)
        {
            // Construct the branch
            for (size_t j = 0; j < node->edgeNodes.size(); ++j)
            {
                auto& edgeNode = node->edgeNodes[j];
                {
                    // If the edgeNode is visited before, then this branch has been reconstructed
                    if (edgeNode->visited)
                        continue;

                    auto branch = _buildBranch(node, edgeNode);
                    branch->index = branchIndex;
                    branchIndex++;

                    branches.push_back(branch);
                }
            }
        }
    }

    return branches;
}

SkeletonBranch* Skeletonizer::_buildBranch(SkeletonNode* firstNode, SkeletonNode* edgeNode)
{
    SkeletonBranch* branch = new SkeletonBranch();

    firstNode->visited = true;
    edgeNode->visited = true;

    branch->nodes.push_back(firstNode);
    branch->nodes.push_back(edgeNode);

    // If the node along the edge is a terminal, then append this branch to the list and return
    if (edgeNode->edgeNodes.size() == 1)
    {
        return branch;
    }

    // If the node along the edge is a branching node, then append this branch to the list and
    // consider the next branch
    else if (edgeNode->edgeNodes.size() > 2)
    {
        return branch;
    }

    // Otherwise, construct the branches from connected segments
    else
    {
        // The previous node is the first node
        SkeletonNode *previousNode = firstNode;

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

            branch->nodes.push_back(currentNode);
            currentNode->visited = true;

            if (!(currentNode->edgeNodes.size() == 2))
                break;
        }
    }

    return branch;

}



}
