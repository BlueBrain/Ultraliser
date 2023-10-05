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

#include<ctime>
#include "Skeletonizer.h"
#include "SkeletonizerUtils.h"
#include <algorithms/skeletonization/thinning/Neighbors.hh>
#include <math/Vector.h>
#include <data/meshes/simple/TriangleOperations.h>
#include <utilities/Range.h>
#include <data/volumes/voxels/NodeVoxel.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(Volume* volume, const bool &useAcceleration)
    : _volume(volume)
    , _useAcceleration(useAcceleration)
{
    // Mesh bounding box
    _pMinMesh = volume->getPMin();
    _pMaxMesh = volume->getPMax();
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    // TODO: Verify the volume point cloud
    // Volume bounding box
    _pMinVolume = Vector3f(0.f);
    _pMaxVolume = Vector3f((volume->getWidth() - 1) * 1.f,
                           (volume->getHeight() - 1) * 1.f,
                           (volume->getDepth() - 1) * 1.f);
    _boundsVolume = _pMaxVolume;
    _centerVolume = 0.5 * _boundsVolume;

    // Mesh to volume scale factor
    _scaleFactor = _boundsMesh / _boundsVolume;
}

void Skeletonizer::initialize()
{
     LOG_TITLE("Neuron Skeletonization");
     LOG_STATUS("Initialization - Building Structures");

     // Start the timer
     TIMER_SET;

     // Compute the shell points either natively or by using the acceleration structures
     if (_useAcceleration)
     {
         // Build the ThinningVoxels acceleration structure from the input solid volume
         // NOTE: We do not rebuild the ThinningVoxels structure!
         auto thinningVoxels = _volume->getThinningVoxelsList(false);

         // Compute the surface shell from the pre-built ThinningVoxels structure
         _computeShellPointsUsingAcceleration(thinningVoxels);
     }
     else
     {
         _computeShellPoints();
     }

     LOG_STATUS_IMPORTANT("Initialization Stats.");
     LOG_STATS(GET_TIME_SECONDS);
}

void Skeletonizer::_scaleShellPoints()
{
    // Initialize the time
    TIMER_SET;

    // TODO: Adjust the voxel slight shift
    // Adjust the locations of the shell points taking into consideration the mesh coordinates
    PROGRESS_SET;
    LOOP_STARTS("Mapping Shell Points");
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

        LOOP_PROGRESS(PROGRESS, _shellPoints.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Skeletonizer::_computeShellPointsUsingAcceleration(ThinningVoxelsUI16List &thinningVoxels)
{
    // Initialize the timer
    TIMER_SET;

    PROGRESS_SET;
    LOOP_STARTS("Computing Shell Points *");
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        auto& voxel = thinningVoxels[i];
        if (_volume->isBorderVoxel(voxel->x, voxel->y, voxel->z))
        {
            voxel->border = true;
        }

        LOOP_PROGRESS(PROGRESS, thinningVoxels.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Add all the obtained voxels in a single list
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        const auto& voxel = thinningVoxels[i];
        if (voxel->border)
        {
            _shellPoints.push_back(Vector3f(voxel->x, voxel->y, voxel->z));
        }
    }

    // Scale the shell points to match the extent of the input data
    _scaleShellPoints();
}

void Skeletonizer::_computeShellPoints()
{
    // Initialize the time
    TIMER_SET;

    // Search for the border voxels (the shell voxels) of the volume
    std::vector< std::vector< Vec3ui_64 > > perSlice = _volume->searchForBorderVoxels();

    // Concatinate the points in a single list
    for (size_t i = 0; i < perSlice.size(); ++i)
    {
        for (size_t j = 0; j < perSlice[i].size(); ++j)
        {
            const auto voxel = perSlice[i][j];
            _shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSlice[i].clear();
        perSlice[i].shrink_to_fit();
    }
    perSlice.clear();
    perSlice.shrink_to_fit();

    // Scale the shell points to match the extent of the input data
    _scaleShellPoints();
}

void Skeletonizer::_applyVolumeThinning()
{
    LOG_STATUS("Volume Thinning");

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    TIMER_SET;
    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        size_t numberDeletedVoxels = _volume->deleteCandidateVoxelsParallel(thinningKernel);

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

void Skeletonizer::_applyVolumeThinningUsingAcceleration()
{
    LOG_STATUS("Volume Thinning *");

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    auto thinningVoxels = _volume->getThinningVoxelsList(false);

    TIMER_SET;
    LOOP_STARTS("Thinning Loop");
    LOOP_PROGRESS(0, 100);
    while(1)
    {
        // Delete the border voxels based on the ThinningVoxels acceleration structure
        size_t numberDeletedVoxels = _volume->deleteBorderVoxelsUsingThinningVoxels(
                    thinningKernel, thinningVoxels);

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

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxels()
{
    if (_useAcceleration)
        return _extractNodesFromVoxelsUsingAcceleration();
    else
        return _extractNodesFromVoxelsUsingSlicing();
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsNaive()
{
    TIMER_SET;
    LOG_STATUS("Mapping Voxels to Nodes");

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;

    // A list of filled voxels to compute the elements in parallel
    std::vector< Vec4ui_64 > indicesFilledVoxels;

    // Search the filled voxels in the volume
    PROGRESS_SET;
    LOOP_STARTS("Detecting Filled Voxels");
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

                    Vec4ui_64 index(i, j, k, voxelIndex);
                    indicesFilledVoxels.push_back(index);

                    // Mapper from voxel to node indices
                    indicesMapper.insert(std::pair< size_t, size_t >(voxelIndex, nodeIndex));

                    // New node
                    nodeIndex++;
                }
            }
        }
        LOOP_PROGRESS(i, _volume->getWidth());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Resize the nodes
    _nodes.resize(indicesFilledVoxels.size());

    PROGRESS_RESET;
    LOOP_STARTS("Building Graph Nodes");
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < indicesFilledVoxels.size(); ++n)
    {
        const size_t i = indicesFilledVoxels[n].x();
        const size_t j = indicesFilledVoxels[n].y();
        const size_t k = indicesFilledVoxels[n].z();
        const size_t voxelIndex = indicesFilledVoxels[n].w();

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
        _nodes[n] = new SkeletonNode(voxelIndex, nodePosition, voxelPosition);

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, indicesFilledVoxels.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clear the auxiliary lists
    indicesFilledVoxels.clear();
    indicesFilledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingSlicing()
{
    LOG_STATUS("Mapping Voxels to Nodes");
    TIMER_SET;

    struct FilledVoxel
    {
        size_t i, j, k, idx;

        FilledVoxel(size_t ii, size_t jj, size_t kk, size_t index)
        { i = ii; j = jj; k = kk; idx = index; }
    };

    typedef std::vector< FilledVoxel > FilledVoxels;

    // Make a per-slice list
    std::vector< FilledVoxels > allFilledVoxels;
    allFilledVoxels.resize(_volume->getWidth());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        // Get a reference to the per-slice list
        auto& perSlice = allFilledVoxels[i];
        for (size_t j = 0; j < _volume->getHeight(); ++j)
        {
            for (size_t k = 0; k < _volume->getDepth(); ++k)
            {
                if (_volume->isFilled(i, j, k))
                {
                    perSlice.push_back(
                        FilledVoxel(i, j, k, _volume->mapTo1DIndexWithoutBoundCheck(i, j, k)));
                }
            }
        }
    }

    // Put them in a single list
    FilledVoxels filledVoxels;
    for (size_t i = 0; i < allFilledVoxels.size(); ++i)
    {
        if (allFilledVoxels[i].size() > 0)
        {
            filledVoxels.insert(filledVoxels.end(),
                                allFilledVoxels[i].begin(), allFilledVoxels[i].end());
            allFilledVoxels[i].clear();
            allFilledVoxels[i].shrink_to_fit();
        }
    }
    allFilledVoxels.clear();
    allFilledVoxels.shrink_to_fit();

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;
    for (size_t i = 0; i < filledVoxels.size(); ++i)
    {
        indicesMapper.insert(std::pair< size_t, size_t >(filledVoxels[i].idx, i));
    }

    // Resize the nodes
    _nodes.resize(filledVoxels.size());

    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < filledVoxels.size(); ++n)
    {
        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(filledVoxels[n].i * 1.f,
                               filledVoxels[n].j * 1.f,
                               filledVoxels[n].k * 1.f);

        // Get a point in the same coordinate space of the mesh
        Vector3f nodePosition(voxelPosition);
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(filledVoxels[n].idx, nodePosition, voxelPosition);

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, filledVoxels.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    filledVoxels.clear();
    filledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingAcceleration()
{
    LOG_STATUS("Mapping Voxels to Nodes *");
    TIMER_SET;

    auto thinningVoxels = _volume->getThinningVoxelsList(false);

    // Construct the NodeVoxels list
    LOOP_STARTS("Constructing Node Voxels");
    PROGRESS_SET;
    NodeVoxelsUI16 nodeVoxels;
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        auto& voxel = thinningVoxels[i];
        if (voxel->active) // The inactive voxels have been deactivated in during the thinning
        {
            NodeVoxelUI16 nodeVoxel;
            nodeVoxel.x = voxel->x;
            nodeVoxel.y = voxel->y;
            nodeVoxel.z = voxel->z;
            nodeVoxel.index = _volume->mapTo1DIndexWithoutBoundCheck(voxel->x, voxel->y, voxel->z);
            nodeVoxels.push_back(nodeVoxel);
        }
    }
    LOG_STATS(GET_TIME_SECONDS);

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;
    for (size_t i = 0; i < nodeVoxels.size(); ++i)
    {
        indicesMapper.insert(std::pair< size_t, size_t >(nodeVoxels[i].index, i));
    }

    // Resize the nodes to the corresponding size of the NodeVoxels list
    _nodes.resize(nodeVoxels.size());

    PROGRESS_RESET;
    LOOP_STARTS("Building Graph Nodes");
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < nodeVoxels.size(); ++n)
    {
        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(nodeVoxels[n].x * 1.f,
                               nodeVoxels[n].y * 1.f,
                               nodeVoxels[n].z * 1.f);

        // Get a point in the same coordinate space of the mesh
        Vector3f nodePosition(voxelPosition);
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(nodeVoxels[n].index, nodePosition, voxelPosition);

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, nodeVoxels.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return indicesMapper;
}

void Skeletonizer::_inflateNodesUsingAcceleration()
{
    // TODO: Adrien to fill this code
    // FOr the moment, use the old functoin until Adrien fills this section
    _inflateNodes();
}

void Skeletonizer::_inflateNodes()
{
    TIMER_SET;
    LOG_STATUS("Inflating Graph Nodes - Mapping to Surface");

    // Compute the approximate radii of all the nodes in the graph, based on the minimum distance
    std::vector< float > nodesRadii;
    nodesRadii.resize(_nodes.size());

    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        float minimumDistance = std::numeric_limits< float >::max();
        for (size_t j = 0; j < _shellPoints.size(); ++j)
        {
            const float distance = (_nodes[i]->point - _shellPoints[j]).abs();
            if (distance < minimumDistance) { minimumDistance = distance; }
        }

        // TODO: Make some logic to detect the actual radius based on the voxel size
        if (minimumDistance > 0.01)
        {
            _nodes[i]->radius = minimumDistance * 1.2;
            nodesRadii[i] = minimumDistance * 1.2;
        }
        else
        {
            _nodes[i]->radius = 0.1;
            nodesRadii[i] = 0.1;
        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, _nodes.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Obtain the node with the largest radius, candidate for soma
    const auto iterator = std::max_element(std::begin(nodesRadii), std::end(nodesRadii));
    const auto& index = std::distance(std::begin(nodesRadii), iterator);
    const auto& largestNode = _nodes[index];

    // Clear the auxiliary list
    nodesRadii.clear();
    nodesRadii.shrink_to_fit();
}

void Skeletonizer::_connectNodes(const std::map< size_t, size_t >& indicesMapper)
{
    TIMER_SET;
    LOG_STATUS("Connecting Graph Nodes");

    /// TODO: Makre sure that there are no race conditions in the code
    LOOP_STARTS("Building & Linking Edges");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        // Check if the node has been visited before
        SkeletonNode* node = _nodes[i];

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
                node->edgeNodes.push_back(_nodes[nIndex]);
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;

        if (connectedEdges > 2)
            node->branching = true;

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, _nodes.size());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void Skeletonizer::_removeTriangleLoops()
{    
    TIMER_SET;
    LOG_STATUS("Removing Triangle Loops");

    const size_t currentNodesSize = _nodes.size();
    size_t progress = 0;
    for (size_t i = 0; i < currentNodesSize; ++i)
    {
        LOOP_PROGRESS(progress, currentNodesSize);
        progress++;

        if (_nodes[i]->branching)
        {
            SkeletonNodes sideNodes;
            if (isTriangleNode(_nodes[i], sideNodes))
            {
                if (_nodes[i]->visited) continue;

                auto& n1 = _nodes[i];
                auto& n2 = sideNodes[0];
                auto& n3 = sideNodes[1];

                // Collapse a triangle into a single node
                collapseTriangleIntoNode(_nodes, n1, n2, n3);

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
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        _nodes[i]->visited = false;
    }
}

void Skeletonizer::constructGraph()
{
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxels();

    // Assign accurate radii to the nodes of the graph, i.e. inflate the nodes
    _inflateNodes();

    // Connect the nodes to construct the edges of the graph
    _connectNodes(indicesMapper);

    // Remove the triangular configurations
    _removeTriangleLoops();

    // Re-index the samples, for simplicity
    OMP_PARALLEL_FOR for (size_t i = 1; i <= _nodes.size(); ++i) { _nodes[i - 1]->index = i; }
}

void Skeletonizer::_buildBranchesFromNodes(const SkeletonNodes& nodes)
{
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
                        _branches.push_back(branch);
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
                        _branches.push_back(branch);
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
                            _branches.push_back(branch);
                        }

                    }
                }
            }
        }
    }
}

std::vector< Vector3f > Skeletonizer::getShellPoints()
{
    return _shellPoints;
}

}
