/***************************************************************************************************
 * Copyright (c) 2016 - 2024
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
#include <algorithms/mcs/DualMarchingCubes.h>
#include <algorithms/skeletonization/graphs/Graphs.h>
#include <algorithms/skeletonization/SkeletonizerUtils.h>

namespace Ultraliser
{
NeuronSkeletonizer::NeuronSkeletonizer(Volume* volume,
                                       const bool &removeSpines,
                                       const bool &useAcceleration,
                                       const bool &debugSkeleton,
                                       const std::string debuggingPrefix)
    : Skeletonizer(volume, useAcceleration, debugSkeleton, debuggingPrefix)
    , _removeSpines(removeSpines)
{
    /// EMPTY CONSTRUCTOR
}

void NeuronSkeletonizer::skeletonizeVolumeToCenterLines()
{
    if (_useAcceleration) _applyVolumeThinningUsingAcceleration(); else _applyVolumeThinning();
}

void NeuronSkeletonizer::constructGraph()
{
    /// Extract the nodes of the skeleton from the center-line "thinned" voxels and return a
    /// mapper that maps the indices of the voxels in the volume and the nodes in the skeleton
    auto indicesMapper = _extractNodesFromVoxels();

    /// DEBUG: Export the nodes file
    if (_debugSkeleton && _debuggingPrefix != NONE) { _exportGraphNodes(_debuggingPrefix); }

    /// Connect the nodes of the skeleton to construct its edges. This operation will not connect
    /// any gaps, it will just connect the nodes extracted from the voxels.
    _connectNodesToBuildEdges(indicesMapper);

    /// Inflate the nodes, i.e. adjust their radii
    _inflateNodes();

    /// DEBUG: Export the inflated nodes file
    if (_debugSkeleton && _debuggingPrefix != NONE)
    { _exportGraphNodes(_debuggingPrefix + "-inflated"); }

    /// Add a virtual soma node, until the soma is reconstructed later
    _addSomaNode();

    /// Re-index the nodes, for simplicity, i.e. the index of the node represents its location in
    /// the _nodes list
    // OMP_PARALLEL_FOR for (size_t i = 1; i <= _nodes.size(); ++i) { _nodes[i - 1]->index = i; }

    /// Segmentthe soma mesh from the branches
    _segmentSomaMesh();

    /// Identify the somatic nodes in the skeleton
    _identifySomaticNodes();

    /// DEBUG: Export the somatic branches
    if (_debugSkeleton && _debuggingPrefix != NONE)
    { _exportBranches(_debuggingPrefix + "-1", SkeletonBranch::SOMATIC); }

    // Verify graph connectivity
    _verifyGraphConnectivityToClosestPartition(_edges);

    /// In the old approach, we needed to remove the triangle loops, but thanks to Foni's algorithm
    /// these loops are removed automatically during the path construction phase
    // _removeTriangleLoops();

    /// Reconstruct the sections "or branches" from the nodes using the edges data
    _buildBranchesFromNodes(_nodes);

    /// DEBUG: Export all the branches
    if (_debugSkeleton && _debuggingPrefix != NONE)
    { _exportBranches(_debuggingPrefix, SkeletonBranch::ALL); }

    /// Validate the branches, and remove the branches inside the soma, i.e. consider them to
    /// be invalid
    _removeBranchesInsideSoma();

    /// DEBUG: Export the somatic branches
    if (_debugSkeleton && _debuggingPrefix != NONE)
    { _exportBranches(_debuggingPrefix + "-2", SkeletonBranch::SOMATIC); }

    /// DEBUG: Export the valid branches
    if (_debugSkeleton && _debuggingPrefix != NONE)
    { _exportBranches(_debuggingPrefix, SkeletonBranch::VALID); }

    /// Connect all the skeleton branches since the roots have been identified after the
    /// soma segmentation
    _connectBranches();
}

void NeuronSkeletonizer::segmentComponents()
{
    /// Initially, and before constructing the graph, remove the loops between two branching points
    _filterLoopsBetweenTwoBranchingPoints();

    /// Remove the loops at a single branching point, i.e. starting and ending at the same node.
    _filterLoopsAtSingleBranchingPoint();

    /// Reduce the skeleton into a list of SkeletonWeightedEdge's
    SkeletonWeightedEdges weighteEdges = _reduceSkeletonToWeightedEdges();

    /// Get a list of all the branching/terminal nodes within the skeleton from the
    /// SkeletonWeightedEdges list
    SkeletonNodes skeletonBranchingNodes = _selectBranchingNodesFromWeightedEdges(weighteEdges);

    /// Get the soma node index within the weighted graph
    int64_t somaNodeIndex = _getSomaIndexFromGraphNodes(skeletonBranchingNodes);

    // Construct the graph nodes list
    GraphNodes graphNodes = _constructGraphNodesFromSkeletonNodes(skeletonBranchingNodes);

    /// After having the weighted edges and the nodes computed, compute the number of components in
    /// the graph, if the result is more than 1 then then re-connect them to be in a single graph
    auto graph = new Graph(weighteEdges, graphNodes);
    auto components = graph->getComponents();
    if (components.size() == 1)
    {
        LOG_SUCCESS("The skeleton graph has 1 component! OK.");
    }
    else
    {
        LOG_WARNING("The skeleton graph has [ %d ] components!", components.size());
    }

    // Find the shortest paths of all the terminals and get a list of the indices of the active edges
    EdgesIndices edgeIndices = _findShortestPathsFromTerminalNodesToSoma(
                weighteEdges, skeletonBranchingNodes, graphNodes, somaNodeIndex);

    // Construct the GraphBranches from the GraphNodes
    GraphBranches graphBranches = _constructGraphBranchesFromGraphNodes(graphNodes, somaNodeIndex);

    // Construct the hierarchy of the graph
    _constructGraphHierarchy(graphBranches);

    // Construct the hierarchy of the skeleton
    _constructSkeletonHierarchy(graphBranches);

    // Invalidate the inactive branches
    _detectInactiveBranches(weighteEdges, edgeIndices);

    // Adkjust the soma radius
    _adjustSomaRadius();

    // Update all the parents
    _updateParents();

    // Remove the spines from the skeleton
    if (_removeSpines)
    {
        // Remove root spines around the soma
        _removeRootSpines();

        // Remove the branches that have 2 samples along the terminals
        for(size_t i = 0; i < 25; ++i)
        {
            _filterShortTerminalBranches();
            _mergeBranchesAfterFilteringSpines();
        }

        /// DEBUG: Export the somatic branches
        if (_debugSkeleton && _debuggingPrefix != NONE)
        { _exportBranches(_debuggingPrefix, SkeletonBranch::TWO_SAMPLE_AND_VALID); }

        _detectBasePaths();

        for(size_t i = 0; i < 5; ++i)
        {
            // Filter the synapses
            _detectSpines();
            _mergeBranchesAfterFilteringSpines();
        }

        // Remove the branches that have 2 samples along the terminals
        for(size_t i = 0; i < 5; ++i)
        {
            // Filter the terminal branches that have two branches
            _filterShortTerminalBranches();
            _mergeBranchesAfterFilteringSpines();
        }

        // Segment the spines and label (count) them
        segmentSpines();

        /// DEBUG: Export the somatic branches
        if (_debugSkeleton && _debuggingPrefix != NONE)
        { _exportBranches(_debuggingPrefix, SkeletonBranch::SPINE); }

        /// DEBUG: Export the spine locations
        if (_debugSkeleton && _debuggingPrefix != NONE)
        { _exportSpineLocations(_debuggingPrefix); }

        /// DEBUG: Export the spine extents
        if (_debugSkeleton && _debuggingPrefix != NONE)
        { _exportSpineExtents(_debuggingPrefix); }
    }
}

void NeuronSkeletonizer::_detectBasePaths()
{
    for (size_t i = 0; i < _roots.size(); ++i)
    {
        SkeletonBranches terminals;
        getTerminals(_roots[i], terminals);

        for (size_t j = 0; j < terminals.size(); ++j)
        {
            auto& terminal = terminals[j];
            SkeletonBranches path;
            std::vector< size_t > pathIndices;
            constructPathFromBranchToSoma(terminal, path, pathIndices);
        }
    }
}

void NeuronSkeletonizer::_findClosestNodesInTwoPartitions(GraphComponent& partition1,
                                                          GraphComponent& partition2,
                                                          size_t* partition1NodeIndex,
                                                          size_t* partition2NodeIndex,
                                                          float* distance)
{
    // A parameter to store the shotrest distances between the nodes
    float shortestDistance = std::numeric_limits<float>::max();
    size_t indexNode1, indexNode2;

    // Iterate over the nodes of the first partition
    for (size_t i = 0; i < partition1.size(); ++i)
    {
        // Reference to the first node
        auto _node1 = _nodes[partition1[i]];

        // Iterate over the nodes of the second partition
        for (size_t j = 0; j < partition2.size(); ++j)
        {
            // Reference to the second node
            auto _node2 = _nodes[partition2[j]];

            // Calculate the distance
            const auto distance = _node1->point.distance(_node2->point);

            // If the calculated distance is less than the shortest distance, update
            if (distance < shortestDistance)
            {
                shortestDistance = distance;
                indexNode1 = i;
                indexNode2 = j;
            }
        }
    }

    *distance = shortestDistance;
    *partition1NodeIndex = indexNode1;
    *partition2NodeIndex = indexNode2;
}

void NeuronSkeletonizer::_connectPartition(GraphComponents& partitions,
                                           const size_t& partitionIndex,
                                           SkeletonEdges &edges)
{
    // A reference to the primary partition
    auto primaryPartition = partitions[partitionIndex];

    size_t partition1NodeIndex;
    size_t partition2NodeIndex;
    float shortestDistance = std::numeric_limits<float>::max();
    size_t secondaryPartitionIndex;

    for (size_t i = 0; i < partitions.size(); ++i)
    {
        // Skip the same partition
        if (i == partitionIndex) continue;

        // A reference to the secondary partition
        auto secondaryPartition = partitions[i];

        size_t _node1Index, _node2Index;
        float distance;

        _findClosestNodesInTwoPartitions(primaryPartition, secondaryPartition,
                                         &_node1Index, &_node2Index, &distance);

        if (distance < shortestDistance)
        {
            shortestDistance = distance;
            partition1NodeIndex = _node1Index;
            partition2NodeIndex = _node2Index;
            secondaryPartitionIndex = i;
        }
    }

    // Add the missing connectivity information
    auto primaryNode = _nodes[primaryPartition[partition1NodeIndex]];
    auto secondaryNode = _nodes[partitions[secondaryPartitionIndex][partition2NodeIndex]];

    primaryNode->edgeNodes.push_back(secondaryNode);
    secondaryNode->edgeNodes.push_back(primaryNode);

    SkeletonEdge* edge = new SkeletonEdge(edges.size(), primaryNode, secondaryNode);
    edges.push_back(edge);
}


void NeuronSkeletonizer::_verifyGraphConnectivityToClosestPartition(SkeletonEdges &edges)
{
    while (true)
    {
        // Construct the graph
        auto graph = new Graph(edges, _nodes);

        // Get the number of partitions
        auto components = graph->getComponents();

        if (components.size() == 1)
        {
            LOG_SUCCESS("The skeleton graph has 1 component! OK.");
            return;
        }
        else
        {
            LOG_WARNING("The skeleton graph has [ %d ] components! Running the Connectomics Algorithm",
                        components.size());

            // Add the partition
            _connectPartition(components, 0, edges);

            // Updating the branching and terminal nodes
            for (size_t i = 0; i < _nodes.size(); ++i)
            {
                // Check if the node has been visited before
                SkeletonNode* node = _nodes[i];

                if (node->edgeNodes.size() == 1)
                    node->terminal = true;
                else
                    node->terminal = false;

                if (node->edgeNodes.size() > 2)
                    node->branching = true;
                else
                    node->branching = false;
            }
        }
    }
}

void NeuronSkeletonizer::_verifyGraphConnectivityToMainPartition(GraphComponents &components,
                                                                 SkeletonEdges &edges)
{
    // Each component in the GraphComponents is simply a list of nodes
    // Find the largest partition
    size_t primaryPartitionIndex = 0;
    size_t numberNodesPrimaryPartition = 0;
    for (size_t i = 0; i < components.size(); ++i)
    {
        if (components[i].size() > numberNodesPrimaryPartition)
        {
            primaryPartitionIndex = i;
            numberNodesPrimaryPartition = components[i].size();
        }
    }

    // Get a reference to the primary partition
    GraphComponent primaryPartition = components[primaryPartitionIndex];

    // Construct a list of the secondary partitions
    GraphComponents secondaryPartitions;
    for (size_t i = 0; i < components.size(); ++i)
    {
        if (i == primaryPartitionIndex) continue;

        secondaryPartitions.push_back(components[i]);
    }


    // Find the connections between each secondary partition and the parimary partition
    for (size_t i = 0; i < secondaryPartitions.size(); ++i)
    {
        auto secondaryPartition = secondaryPartitions[i];

        // Find the indices of the connecting nodes
        size_t closestPrimaryNodeIndex;
        size_t closesetSecondaryNodeIndex;
        float shortestDistance = 1e32;

        for (size_t j = 0; j < secondaryPartition.size(); ++j)
        {
            auto secondaryNodeIndex = secondaryPartition[j];
            auto secondaryNode = _nodes[secondaryNodeIndex];

            for (size_t k = 0; k < primaryPartition.size(); ++k)
            {
                auto primaryNodeIndex = primaryPartition[k];
                auto primaryNode = _nodes[primaryNodeIndex];

                const auto distance = primaryNode->point.distance(secondaryNode->point);
                if (distance < shortestDistance)
                {
                    shortestDistance = distance;
                    closestPrimaryNodeIndex = primaryNodeIndex;
                    closesetSecondaryNodeIndex = secondaryNodeIndex;
                }
            }
        }

        // The primary and secondary nodes are connected
        auto primaryNode = _nodes[closestPrimaryNodeIndex];
        auto secondaryNode = _nodes[closesetSecondaryNodeIndex];

        primaryNode->edgeNodes.push_back(secondaryNode);
        secondaryNode->edgeNodes.push_back(primaryNode);

        SkeletonEdge* edge = new SkeletonEdge(edges.size(), primaryNode, secondaryNode);
        edges.push_back(edge);
    }

    // Updating the branching and terminal nodes
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        // Check if the node has been visited before
        SkeletonNode* node = _nodes[i];

        if (node->edgeNodes.size() == 1)
            node->terminal = true;
        else
            node->terminal = false;

        if (node->edgeNodes.size() > 2)
            node->branching = true;
        else
            node->branching = false;
    }
}

void NeuronSkeletonizer::_verifyGraphConnectivity(SkeletonEdges& edges)
{
    // Since we have all the nodes and the edges, we can verify if the graph is conected or not
    auto graph = new Graph(edges, _nodes);

    auto components = graph->getComponents();

    if (components.size() == 1)
    {
        LOG_SUCCESS("The skeleton graph has 1 component! OK.");
    }
    else
    {
        LOG_WARNING("The skeleton graph has [ %d ] components!", components.size());

        // Verify the graph connectivity to the main partition
        _verifyGraphConnectivityToMainPartition(components, edges);

        auto newGraph = new Graph(edges, _nodes);

        auto newComponents = newGraph->getComponents();

        if (newComponents.size() == 1)
        {
            LOG_SUCCESS("The skeleton graph has now 1 component! OK.");
        }
        else
        {
            LOG_WARNING("The skeleton graph has still [ %d ] components!", newComponents.size());
        }
    }
}

void NeuronSkeletonizer::_addSomaNode()
{
    _somaNode = new SkeletonNode();
    _somaNode->index = _nodes.back()->index + 1;

    // By default, the soma node is actually the soma
    _somaNode->isSoma = true;

    // The somatic node is considered inside the soma in the processing
    _somaNode->insideSoma = true;

    // Initially, we set the soma node to some value that does not make any conflict
    _somaNode->radius = 0.1;
    _nodes.push_back(_somaNode);
}

void NeuronSkeletonizer::_segmentSomaMesh()
{
    LOG_STATUS("Segmenting Soma");

    // The _somaProxyMesh should be a complex geometry containing overlapping spheres that would define
    // its structure
    _somaProxyMesh = new Mesh();

    // For every node in the skeleton, if the radius is greater than 2.0, then this is a candidate
    // for the soma sample
    Vector3f estimatedSomaCenter(0.f);
    size_t numberSamples = 0;

    // Collecting a subset of the samples, not all of them are needed
    std::map <size_t, float> interSomaticNodes;

    size_t interSomaticNodesCount = 0;
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];

        /// TODO: This is a magic value, it works now, but we need to find an optimum value based
        /// on some statistical analysis.
        if (node->radius >= 2.0)
        {
            interSomaticNodes.insert({node->index, node->radius});
            interSomaticNodesCount++;
        }
    }

    // Sort the nodes by radius
    std::vector< std::pair< size_t , float > > pairsVector = sortIndexRadiusMap(interSomaticNodes);

    // Reverse
    std::reverse(pairsVector.begin(), pairsVector.end());

    /// TODO: Fix me
    size_t numberSelectedNodes = interSomaticNodesCount;

    TIMER_SET;
    LOOP_STARTS("Detecting Soma Nodes");
    for (size_t i = 0; i < numberSelectedNodes; ++i)
    {
        auto& node0 = _nodes[pairsVector[i].first];

        Mesh* sample = new IcoSphere(3);
        sample->scale(node0->radius, node0->radius, node0->radius);
        sample->translate(node0->point);
        sample->map(_shellPoints, false);

        _somaProxyMesh->append(sample);
        sample->~Mesh();

        estimatedSomaCenter += node0->point;
        numberSamples++;

       LOOP_PROGRESS(i, numberSelectedNodes);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Normalize
    estimatedSomaCenter /= numberSamples;

    // Update the location of the soma point, the radius will be updated after detecting root arbors
    _somaNode->point = estimatedSomaCenter;
    _somaNode->radius = 0.5; /// TODO: Fix me
}

void NeuronSkeletonizer::_identifySomaticNodes()
{
    TIMER_SET;
    LOG_STATUS("Identifying Somatic Nodes");

    // Clear the volume to start using the same grid for filling the soma
    _volume->clear();

    // Rasterize the soma mesh
    LOOP_STARTS("Surface Voxelization");
    _volume->surfaceVoxelization(_somaProxyMesh, false, false);
    LOG_STATS(GET_TIME_SECONDS);

    // This is a list of 1D indices of all the voxels that are inside the soma
    std::vector< size_t > insideSomaVoxels;

    if (true)
    {
#ifdef ULTRALISER_DEBUG
        // Compute the active region bounds for validation
        Bounds3D_ui64 activeRegionBounds = _volume->getActiveRegionBounds();
        std::cout << "Active Bounds : "
                  << activeRegionBounds.x1 << ", " << activeRegionBounds.x2 << ","
                  << activeRegionBounds.y1 << ", " << activeRegionBounds.y2 << ", "
                  << activeRegionBounds.z1 << ", " << activeRegionBounds.z2 << "\n";
#endif

        // Compute the soma mesh bounding box (in coordinate space)
        Vector3f pMinSomaMesh, pMaxSomaMesh;
        _somaProxyMesh->computeBoundingBox(pMinSomaMesh, pMaxSomaMesh);

        // Compute the soma volume bounding box ()
        auto somaVolumeBounds = _volume->getROIBounds(pMinSomaMesh, pMaxSomaMesh);

#ifdef ULTRALISER_DEBUG
        std::cout << "Soma Bounds * : "
                  << somaVolumeBounds.x1 << "," << somaVolumeBounds.x2 << ","
                  << somaVolumeBounds.y1 << ", " << somaVolumeBounds.y2 << ", "
                  << somaVolumeBounds.z1 << ", " << somaVolumeBounds.z2 << "\n";
#endif

        // Verify the bounds
        auto x1 = static_cast< size_t >(somaVolumeBounds.x1);
        auto x2 = static_cast< size_t >(somaVolumeBounds.x2);
        if (x2 >= _volume->getWidth()) x2 = _volume->getWidth() - 1;

        auto y1 = static_cast< size_t >(somaVolumeBounds.y1);
        auto y2 = static_cast< size_t >(somaVolumeBounds.y2);
        if (y2 >= _volume->getHeight()) y2 = _volume->getHeight() - 1;

        auto z1 = static_cast< size_t >(somaVolumeBounds.z1);
        auto z2 = static_cast< size_t >(somaVolumeBounds.z2);
        if (z2 >= _volume->getDepth()) z2 = _volume->getDepth() - 1;

        // Apply the solid voxelization only to the selected ROI to avoid flood-filling large slices
        TIMER_RESET;
        LOOP_STARTS("Solid Voxelization ROI *");
        _volume->solidVoxelizationROI(Volume::SOLID_VOXELIZATION_AXIS::X, x1, x2, y1, y2, z1, z2,
                                      false);
        LOG_STATS(GET_TIME_SECONDS);

        // Compute the voxels that are inside the soma
        TIMER_RESET;
        LOOP_STARTS("Finding Somatic Voxels");
        for (size_t i = x1; i <= x2; ++i)
        {
            for (size_t j = y1; j <= y2; ++j)
            {
                for (size_t k = z1; k <= z2; ++k)
                {
                    const size_t& index = _volume->mapTo1DIndexWithoutBoundCheck(i, j, k);
                    if (_volume->isFilled(index))
                    {
                        insideSomaVoxels.push_back(index);
                    }
                }
            }
        }
        LOG_STATS(GET_TIME_SECONDS);
    }

    // This branch is used for debugging just in case any optimizations fail
    else
    {
        // Apply solid voxelization on the entire slice
        TIMER_RESET;
        LOOP_STARTS("Solid Voxelization");
        _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::X);
        LOG_STATS(GET_TIME_SECONDS);

        // Compute the voxels that are inside the soma
        TIMER_RESET;
        LOOP_STARTS("Finding Somatic Voxels");
        for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
        {
            if (_volume->isFilled(i)) { insideSomaVoxels.push_back(i); }
        }
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Find out the nodes that are inside the soma
    TIMER_RESET;
    LOOP_STARTS("Mapping Voxels to Nodes");
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];
        // size_t key = _volume->mapTo1DIndexWithoutBoundCheck(
        // node->voxel.x(), node->voxel.y(), node->voxel.z());
        size_t key = node->voxelIndex;
        if (std::find(insideSomaVoxels.begin(), insideSomaVoxels.end(), key)
                != insideSomaVoxels.end())
        {
            // The node is inside the soma
            node->insideSoma = true;
        }
    }
    LOG_STATS(GET_TIME_SECONDS);

    /// TODO: The volume is safe to be deallocated
    _volume->~Volume();
    _volume = nullptr;
}

void NeuronSkeletonizer::_reconstructSomaMeshFromProxy()
{
    // Get relaxed bounding box to build the volume
    Vector3f pMinInput, pMaxInput;
    _somaProxyMesh->computeBoundingBox(pMinInput, pMaxInput);
    const auto& meshBoundingBox = pMaxInput - pMinInput;

    // Construct the volume
    Volume* somaVolume = new Volume(pMinInput, pMaxInput, 128, 0.1);

    // Voxelize the proxy mesh
    somaVolume->surfaceVoxelization(_somaProxyMesh, false, false, 1.0);
    somaVolume->solidVoxelization();

    // Use DMC to reconstruct a mesh
    _somaMesh = DualMarchingCubes::generateMeshFromVolume(somaVolume);

    // Free the volume
    delete somaVolume;

    // Optimize the mesh
    _somaMesh->smoothSurface(10);
    _somaMesh->optimizeAdaptively(5, 5, 0.05, 5.0);
}

void NeuronSkeletonizer::_removeBranchesInsideSoma()
{
    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];

        // Get the first and last nodes
        auto& firstNode = branch->nodes.front();
        auto& lastNode = branch->nodes.back();

        // If the first and last nodes of the branch are inside the soma, then it is invalid
        if (branch->nodes.front()->insideSoma && branch->nodes.back()->insideSoma)
        {
            // TODO: What could be other possible cases!
            branch->setInvalid();
            branch->setInsideSoma();
            branch->unsetRoot();
        }
        else
        {
            // Count the number of samples inside the soma
            size_t countSamplesInsideSoma = 0;
            for (size_t j = 0; j < branch->nodes.size(); ++j)
            {
                if (branch->nodes[j]->insideSoma) { countSamplesInsideSoma++; }
            }

            // If the count of the samples located inside the soma is zero, then it is a valid
            // branch, but it is not a root branch
            if (countSamplesInsideSoma == 0)
            {
                branch->unsetRoot();
                branch->setValid();
            }

            // If all the branch nodes are located inside the soma, then it is not valid
            else if (countSamplesInsideSoma == branch->nodes.size())
            {
                branch->unsetRoot();
                branch->setInsideSoma();
                branch->setInvalid();
            }

            // Otherwise, it is a branch that is connected to the soma, partially in the soma
            else
            {
                // If the first node is inside the soma, then annotate the branch
                if (firstNode->insideSoma)
                {
                    SkeletonNodes newNodes;
                    newNodes.push_back(_somaNode);
                    for (size_t k = 0; k < branch->nodes.size(); ++k)
                    {
                        if (!branch->nodes[k]->insideSoma)
                        {
                            newNodes.push_back(branch->nodes[k]);
                        }
                    }

                    branch->nodes.clear();
                    branch->nodes.shrink_to_fit();
                    branch->nodes = newNodes;

                    branch->setRoot();
                    branch->setValid();

                    // Add this branch to the roots
                    _roots.push_back(branch);
                }

                // If the last node is inside the soma, then annotate the branch
                else if (lastNode->insideSoma)
                {
                    SkeletonNodes newNodes;
                    for (size_t k = 0; k < branch->nodes.size(); ++k)
                    {
                        if (!branch->nodes[k]->insideSoma)
                        {
                            newNodes.push_back(branch->nodes[k]);
                        }
                    }
                    newNodes.push_back(_somaNode);

                    branch->nodes.clear();
                    branch->nodes.shrink_to_fit();
                    branch->nodes = newNodes;

                    // Reverse the nodes, to make the first node in front
                    std::reverse(std::begin(branch->nodes), std::end(branch->nodes));

                    branch->setRoot();
                    branch->setValid();

                    // Add this branch to the roots
                    _roots.push_back(branch);
                }
                else
                {
                    LOG_WARNING("Undefined case for the branch identification! Possible Errors!");
                }
            }
        }
    }
}

void NeuronSkeletonizer::_filterLoopsBetweenTwoBranchingPoints()
{
    TIMER_SET;
    LOG_STATUS("Filtering Loops Between Two Branches");
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // Reference to the iBranch
        auto& iBranch = _branches[i];

        // The iBranch must be valid to be processed
        if (iBranch->isValid())
        {
            for (size_t j = 0; j < _branches.size(); ++j)
            {
                // Reference to the jBranch
                auto& jBranch = _branches[j];

                // If the same branch, continue
                if (i == j) continue;

                // The jBranch must be valid as well to be processed
                if (jBranch->isValid())
                {
                    // If the jBranch has the same terminal nodes as iBranch
                    if (iBranch->nodes.front()->index == jBranch->nodes.front()->index &&
                        iBranch->nodes.back()->index  == jBranch->nodes.back()->index  ||
                        iBranch->nodes.front()->index == jBranch->nodes.back()->index  &&
                        iBranch->nodes.back()->index  == jBranch->nodes.front()->index)
                    {
                        // Select the shortest path
                        if (iBranch->nodes.size() < jBranch->nodes.size())
                        {
                            iBranch->setValid();
                            jBranch->setInvalid();
                        }
                        else
                        {
                            iBranch->setInvalid();
                            jBranch->setValid();
                        }
                    }
                }
            }
        }
        LOOP_PROGRESS(i, _branches.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void NeuronSkeletonizer::_filterLoopsAtSingleBranchingPoint()
{
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // Reference to the iBranch
        auto& iBranch = _branches[i];

        // If forming a self-loop, i.e. starts and ends at the same node, it is invalid
        if (iBranch->isLoop())
        {
            iBranch->setInvalid();
        }
    }
}

SkeletonWeightedEdges NeuronSkeletonizer::_reduceSkeletonToWeightedEdges()
{
    TIMER_SET;
    LOG_STATUS("Created Simplified Weighted Graph from Skeleton");

    /// The "thinned" skeleton has a list of edges, but indeed we would like to simplify it to
    /// avoid spending hours traversing it. Therefore, we created a "weighted skeleton", where
    /// each edge in this skeleton is a connection between two branching/terminal nodes and the
    /// weight represents the number of samples "or samples" between the two branching/terminal
    /// nodes. The "weighted edges" must be constructed only for valid branches.
    SkeletonWeightedEdges edges;
    const auto branchesCount = _branches.size();
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // Reference to the current branch
        auto& branch = _branches[i];

        // The branch must be valid to be able to have a valid graph
        if (branch->isValid())
        {
            // Reset the traversal state of the "terminal" nodes of the branch
            branch->nodes.front()->visited = false;
            branch->nodes.back()->visited = false;

            // Create a weighted edge and append it to the list, where the weight is indicated by
            // the number of nodes "or samples" in the branch
            SkeletonWeightedEdge* edge = new SkeletonWeightedEdge(branch);
            edges.push_back(edge);
        }

        LOOP_PROGRESS(i, branchesCount);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Return the resulting edges array that will be used for constructing the graph
    return edges;
}

SkeletonNodes NeuronSkeletonizer::_selectBranchingNodesFromWeightedEdges(
        const SkeletonWeightedEdges& edges)
{
    TIMER_SET;
    LOG_STATUS("Identifying Branching Nodes");

    // Use a new index to label the branching nodes, where the maximum value corresponds to the
    // actual number of the branching nodes in the graph
    int64_t branchingNodeIndex = 0;

    // A list to collect the branching nodes
    SkeletonNodes nodes;
    for (size_t i = 0; i < edges.size(); ++i)
    {
        // The node must be visited once to append it to the @skeletonBranchingNodes list
        auto& edge = edges[i];
        auto node1 = edge->node1;
        auto node2 = edge->node2;

        // First node of the edge
        if (!node1->visited)
        {
            node1->graphIndex = branchingNodeIndex;
            nodes.push_back(node1);
            branchingNodeIndex++;
            node1->visited = true;
        }

        // Second node of the edge
        if (!node2->visited)
        {
            node2->graphIndex = branchingNodeIndex;
            nodes.push_back(node2);
            branchingNodeIndex++;
            node2->visited = true;
        }
        LOOP_PROGRESS(i, edges.size());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return nodes;
}

int64_t NeuronSkeletonizer::_getSomaIndexFromGraphNodes(const SkeletonNodes& nodes) const
{
    int64_t somaNodeIndex = -1;
    for (size_t i = 0; i < nodes.size(); ++i)

    {
        const auto& node = nodes[i];
        if (node->isSoma)
        {
            somaNodeIndex = node->graphIndex;
            break;
        }
    }

    // Return the soma node index
    return somaNodeIndex;
}

GraphNodes NeuronSkeletonizer::_constructGraphNodesFromSkeletonNodes(
        const SkeletonNodes& skeletonNodes)
{
    // A list to conatin all the graph nodes
    GraphNodes graphNodes;

    // Resize it to be allow the parallelization of the list
    graphNodes.resize(skeletonNodes.size());

    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < skeletonNodes.size(); ++i)
    {
        const auto& skeletonNode = skeletonNodes[i];
        graphNodes[i] = (new GraphNode(skeletonNode->graphIndex,
                                       skeletonNode->point,
                                       skeletonNode->index,
                                       skeletonNode->edgeNodes.size()));
    }

    // Return the graph nodes list
    return graphNodes;
}

EdgesIndices NeuronSkeletonizer::_findShortestPathsFromTerminalNodesToSoma(
        SkeletonWeightedEdges& edges,
        SkeletonNodes &skeletonBranchingNodes, GraphNodes &graphNodes,
        const int64_t& somaNodeIndex)
{
    TIMER_SET;
    LOG_STATUS("Identifying Short Paths from Terminals To Soma");

    // Identify the terminal nodes to process the paths in parallel
    SkeletonNodes terminalNodes;
    for (size_t i = 0; i < skeletonBranchingNodes.size(); i++)
    {
        if (skeletonBranchingNodes[i]->terminal)
            terminalNodes.push_back(skeletonBranchingNodes[i]);
    }

    std::vector< EdgesIndices > edgesIndicesList;
    edgesIndicesList.resize(skeletonBranchingNodes.size());

    // Generate the ShortestPathFinder only once for all the path retrival functions
    const ShortestPathFinder pathFinder(edges, skeletonBranchingNodes.size());

    // Search for all the terminal nodes
    size_t numberTerminalNodes = terminalNodes.size();
    std::vector< PathIndices > terminalsToSomaPaths;
    terminalsToSomaPaths.resize(numberTerminalNodes);

    LOOP_STARTS("Searching Terminals-to-soma Paths *");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t iNode = 0; iNode < numberTerminalNodes; ++iNode)
    {
#ifdef REVERSE
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = pathFinder->findPath(terminalNodes[iNode]->graphIndex,
                                                       somaNodeIndex);

        // Reverse the terminal to soma path to have the correct order
        std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());
#else
        // Find the path between the terminal node and the soma node
        terminalsToSomaPaths[iNode] = pathFinder.findPath(
                    somaNodeIndex, terminalNodes[iNode]->graphIndex);
#endif
        LOOP_PROGRESS(PROGRESS, numberTerminalNodes);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    LOOP_STARTS("Updating Graph");
    PROGRESS_RESET;
    for (size_t iNode = 0; iNode < numberTerminalNodes; ++iNode)
    {
        // Get a reference to the EdgesIndices list
        EdgesIndices& perTerminalEdgesIndices = edgesIndicesList[iNode];

#ifdef REVERSE
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = pathFinder->findPath(terminalNodes[iNode]->graphIndex,
                                                       somaNodeIndex);

        // Reverse the terminal to soma path to have the correct order
        std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());
#else
        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = terminalsToSomaPaths[iNode];

        // Find the edges
        for (size_t j = 0; j < terminalToSomaPath.size() - 1; ++j)
        {
            auto currentNodeIndex = terminalToSomaPath[j];
            auto nextNodeIndex = terminalToSomaPath[j + 1];

            // Add the edge indices to the list
            perTerminalEdgesIndices.push_back(EdgeIndex(currentNodeIndex, nextNodeIndex));

            // If the next node index is not in the current node index, then add it
            if (!graphNodes[currentNodeIndex]->isNodeInChildren(nextNodeIndex))
            {
                graphNodes[currentNodeIndex]->children.push_back(graphNodes[nextNodeIndex]);
            }
        }
#endif
        LOOP_PROGRESS(PROGRESS, numberTerminalNodes);
        PROGRESS_UPDATE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clear the terminal nodes list
    terminalNodes.clear();
    terminalNodes.shrink_to_fit();

    // Clear the terminal nodes list
    terminalsToSomaPaths.clear();
    terminalsToSomaPaths.shrink_to_fit();

    // The indices of all the edges that have been traversed
    EdgesIndices edgesIndices;
    PROGRESS_RESET;
    LOOP_STARTS("Composing Path Edges");
    for (size_t i = 0; i < edgesIndicesList.size(); ++i)
    {
        // Get the terminal edge identified per terminal
        auto perTerminalEdgesIndices = edgesIndicesList[i];
        if (perTerminalEdgesIndices.size() > 0)
        {
            // Append them to the edgeIndices list
            edgesIndices.insert(edgesIndices.end(),
                                perTerminalEdgesIndices.begin(), perTerminalEdgesIndices.end());

            // Clean the list per terminal
            perTerminalEdgesIndices.clear();
            perTerminalEdgesIndices.shrink_to_fit();
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Clean the list used to collect the edges in parallel
    edgesIndicesList.clear();
    edgesIndicesList.shrink_to_fit();

    // Return the EdgesIndices list
    return edgesIndices;
}

GraphBranches NeuronSkeletonizer::_constructGraphBranchesFromGraphNodes(
        GraphNodes &graphNodes, const int64_t& somaNodeIndex)
{
    TIMER_SET;
    LOG_STATUS("Constructing Graph Branches");

    // Use a new index to label graph branches
    size_t branchGraphIndex = 0;

    // A list of all the constructed GraphBranches
    GraphBranches graphBranches;

    // Construct the valid branches at the end
    for (size_t i = 0; i < graphNodes.size(); i++)
    {
         LOOP_PROGRESS(i, graphNodes.size());
        if (graphNodes[i]->children.size() > 0)
        {
            // This graph node is always the first node, becuase all the other nodes are children
            const auto& firstNodeIndex = graphNodes[i]->index;
            const auto& firstNodeSkeletonIndex = graphNodes[i]->skeletonIndex;

            for (size_t j = 0; j < graphNodes[i]->children.size(); ++j)
            {
                // This graph node is always the last node, because it is a child node
                const auto& lastNodeIndex = graphNodes[i]->children[j]->index;
                const auto& lastNodeSkeletonIndex = graphNodes[i]->children[j]->skeletonIndex;

                // Search for the branches
                for (size_t k = 0; k < _branches.size(); k++)
                {
                    // Reference to the branch
                    auto& branch = _branches[k];

                    // The branch must be valid
                    if (branch->isValid() &&
                        branch->hasTerminalNodes(firstNodeSkeletonIndex, lastNodeSkeletonIndex))
                    {
                        GraphBranch* gBranch = new GraphBranch(branchGraphIndex);
                        branchGraphIndex++;

                        gBranch->skeletonIndex = branch->index;
                        gBranch->firstNodeIndex = firstNodeIndex;
                        gBranch->firstNodeSkeletonIndex = firstNodeSkeletonIndex;
                        gBranch->lastNodeIndex = lastNodeIndex;
                        gBranch->lastNodeSkeletonIndex = lastNodeSkeletonIndex;
                        if (somaNodeIndex == firstNodeIndex) gBranch->isRoot = true;

                        graphBranches.push_back(gBranch);
                    }
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return graphBranches;
}

void NeuronSkeletonizer::_constructGraphHierarchy(GraphBranches& graphBranches)
{
    TIMER_SET;
    LOG_STATUS("Constructing Graph Hierarchy");

    for (size_t i = 0; i < graphBranches.size(); ++i)
    {
        LOOP_PROGRESS(i, graphBranches.size());
        auto& iBranch = graphBranches[i];

        // Get the last node of the iBranch
        const auto& iLastNodeIndex = iBranch->lastNodeIndex;

        for (size_t j = 0; j < graphBranches.size(); ++j)
        {
            auto& jBranch = graphBranches[j];

            // If the same branch, next branch
            if (iBranch->index == jBranch->index) continue;

            // Get the first node of the jBranch
            const auto& jFirstNodeIndex = jBranch->firstNodeIndex;

            // jBranch is a child of the iBranch if the nodes are the same
            if (iLastNodeIndex == jFirstNodeIndex)
            {
                iBranch->children.push_back(jBranch);
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void NeuronSkeletonizer::_constructSkeletonHierarchy(GraphBranches& graphBranches)
{
    TIMER_SET;
    LOG_STATUS("Constructing Skeleton Hierarchy");

    for(size_t i = 0; i < graphBranches.size(); ++i)
    {
        LOOP_PROGRESS(i, graphBranches.size());

        // Reference to the GraphBranch
        const auto& graphBranch = graphBranches[i];

        // Reference to the corresponding SkeletonBranch
        auto& skeletonBranch = _branches[graphBranch->skeletonIndex];

        // Adjust the direction, i.e. the samples list
        skeletonBranch->adjustDirection(graphBranch->firstNodeSkeletonIndex,
                                        graphBranch->lastNodeSkeletonIndex);

        // Update the tree
        for(size_t j = 0; j < graphBranch->children.size(); ++j)
        {
            skeletonBranch->children.push_back(
                        _branches[graphBranch->children[j]->skeletonIndex]);
            skeletonBranch->logicalChildren.push_back(
                        _branches[graphBranch->children[j]->skeletonIndex]);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void NeuronSkeletonizer::_mergeBranchesAfterFilteringSpines()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // A reference to the branch
        auto& branch = _branches[i];

        // If the branch is spine
        if (branch->isSpine())
        {
            for (size_t j = 0; j < branch->parents.size(); ++j)
            {
                // Get the parent
                auto parent = branch->parents[j];

                std::vector< SkeletonBranch* > newChildren;
                for (size_t k = 0; k < parent->children.size(); ++k)
                {
                    if (parent->children[k]->index == branch->index)
                    {
                        // IGNORE
                    }
                    else
                    {
                        newChildren.push_back(parent->children[k]);
                    }
                }

                parent->children.clear();
                parent->children.shrink_to_fit();

                for (size_t k = 0; k < newChildren.size(); ++k)
                {
                    parent->children.push_back(newChildren[k]);
                }
            }

        }
    }

    /// Set the terminals
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto branch = _branches[i];
        if (branch->isValid() && branch->children.size() == 0)
        {
            branch->setTerminal();
        }
    }
}

void NeuronSkeletonizer::_detectInactiveBranches(SkeletonWeightedEdges& graphEdges,
                                                 EdgesIndices& visitedEdgesIndices)
{
    // For every edge in the graph, verifiy if this edge contains the nodes that are visited
    for (size_t i = 0; i < graphEdges.size(); i++)
    {
        // Get a reference to the edge
        auto& edge = graphEdges[i];

        // If the edge is already set to be visited, there is no need to search
        if (!edge->visited)
        {
            // Search for all the visited edges indices
            for (size_t j = 0; j < visitedEdgesIndices.size(); ++j)
            {
                // Get a reference to a single edge indices pair
                auto& edgeIndices = visitedEdgesIndices[j];

                // Verification
                if (edge->hasTerminalGraphNodes(edgeIndices.first, edgeIndices.second))
                {
                    // Set the edge to be visited
                    edge->visited = true;

                    // There is no need to continue, then break to save performance
                    break;
                }
            }
        }
    }

    // Update the skeleton using the SkeletonWeightedEdges list
    for (size_t i = 0; i < graphEdges.size(); ++i)
    {
        // If the edge is not visited, the invalidate and deactivate it in the skeleton
        if (!graphEdges[i]->visited)
        {
            graphEdges[i]->branch->setInvalid();
            // graphEdges[i]->branch->setInactive();
        }
    }
}

void NeuronSkeletonizer::_adjustSomaRadius()
{
    // Only count the valid roots to normalize the size of the soma
    size_t numberValidRoots = 0;

    // Calculate the actual radius of the soma from the valid roots
    float somaRadius = 0.;
    for(size_t i = 0; i < _branches.size(); ++i)
    {
        // Get the branch
        const auto branch = _branches[i];

        // Ensure that the branch is a root and also a valid one
        if (branch->isRoot() && branch->isValid())
        {
            // The second sample is the first sample of the branch and the first sample is
            // the soma center
            somaRadius += branch->nodes[1]->point.distance(_somaNode->point);
            numberValidRoots++;
        }
    }

    // Normalize and update the soma node
    somaRadius /= numberValidRoots;
    _somaNode->radius = somaRadius;
}

void NeuronSkeletonizer::_updateParent(SkeletonBranch* branch)
{
    for(size_t j = 0; j < branch->children.size(); j++)
    {
        auto& child = branch->children[j];

        // Clear old parents if any
        child->parents.clear();
        child->parents.shrink_to_fit();

        // Add the new parent
        child->parents.push_back(branch);
        child->logicalParents.push_back(branch);

        _updateParent(child);
    }
}

void NeuronSkeletonizer::_updateParents()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];

        if (branch->isValid() && branch->isRoot())
        {
            for(size_t j = 0; j < branch->children.size(); j++)
            {
                auto& child = branch->children[j];

                // Clear old parents if any
                child->parents.clear();
                child->parents.shrink_to_fit();

                // Add the new parent
                child->parents.push_back(branch);
                child->logicalParents.push_back(branch);

                _updateParent(child);
            }
        }
    }
}

void NeuronSkeletonizer::_removeRootSpines()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];
        if (branch->isValid() &&
            branch->isRoot() &&
            branch->isTerminal() &&
            branch->computeLength() < 6.0)
        {
            // This is a spine on the soma
            branch->setSpine();
            branch->setInvalid();
        }
    }
}

size_t NeuronSkeletonizer::_estimateNumberSpineCandidates()
{
    size_t numberSpineCandidates = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];
        if (branch->isValid())
        {
            if (branch->isTerminal())
            {
                if (branch->computeLength() < 6.0)
                {
                    numberSpineCandidates += 1;
                }
            }
        }
    }

    return numberSpineCandidates;
}

void NeuronSkeletonizer::_removeSpinesss()
{

    size_t iterations = 0;
    while (true)
    {
        // Filter terminal candidates that could be spines
        _filterSpineCandidates();
        iterations++;


        auto numberSpinesEstimated = _estimateNumberSpineCandidates();
        if (numberSpinesEstimated == 0 || iterations > 5)
        {
            LOG_SUCCESS("Spines are completely eliminated");
            return;
        }
        else
        {
            LOG_WARNING("Estimated Spines to Eliminate %d", numberSpinesEstimated);
        }
    }
}

void NeuronSkeletonizer::segmentSpines()
{
    // Unset the visited flag for all the spines
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto branch = _branches[i];
        branch->setUnvisited();
    }

    // An index that will be assigned to the spine
    size_t spineIndex = 1;

    /// NOTE: This function must be called after the labeling of the spines.
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto branch = _branches[i];

        // The branch is then labeled as a spine
        if (branch->isSpine())
        {
            // The branch must not be visited before, otherwise it is already segmented
            if (!branch->visited())
            {
                // Do the traversal logical and segment the spine
                _spineRoots.push_back(segmentSpineFromBranch(branch, spineIndex));

                // Next spine
                spineIndex++;
            }
        }
    }

    LOG_SUCCESS("[%ld] spine were segmented!", spineIndex - 1);
}





void NeuronSkeletonizer::reskeletonizeSpines(Mesh* neuronMesh,
                                             const float& voxelsPerMicron,
                                             const float& edgeGap)
{
    if (_spineRoots.size() == 0)
        return;

    for (size_t i = 0; i < _spineRoots.size(); ++i)
    {
        // Get the spine morphology
        auto spine =_spineRoots[i];
        SpineMorphology* spineMorphology = new SpineMorphology(spine);

        // Get the mesh of the spine model
        auto spineModelMesh = spineMorphology->reconstructMesh(voxelsPerMicron, edgeGap);

        // Get relaxed bounding box to build the volume
        Vector3f pMinInput, pMaxInput, inputBB, inputCenter;
        spineMorphology->getBoundingBox(pMinInput, pMaxInput, inputBB, inputCenter);

        // Get the largest dimension
        float largestDimension = inputBB.getLargestDimension();
        size_t resolution = static_cast< size_t >(voxelsPerMicron * largestDimension);

        // Construct the volume
        Volume* volume = new Volume(pMinInput, pMaxInput, resolution, edgeGap, VOLUME_TYPE::BIT, false);

        // Rasterize the morphologies into the volume
        volume->surfaceVoxelization(neuronMesh);
        // volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

        std::stringstream prefixs;
        prefixs << "//home/abdellah/spines/models/spine_" << spine->index;

        volume->project(prefixs.str(), true);

        // Use the DMC algorithm to reconstruct a mesh
        auto originalMesh = DualMarchingCubes::generateMeshFromVolume(volume);
        originalMesh->smoothSurface(10);

        // Map the surface
        spineModelMesh->map(originalMesh);


        // prefixs << "//home/abdellah/spines/models/spine_" << spine->index;
        spineModelMesh->exportMesh(prefixs.str(), true);

        spineMorphology->exportExtents(_debuggingPrefix);
    }

}

void NeuronSkeletonizer::_filterSpineCandidates()
{
    /// NOTE: Spine branches must be terminals and can have branching.
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // The spine is always a valid branch
        auto& branch = _branches[i];
        if (branch->isValid())
        {
            // The spine is always a terminal branch
            if (branch->isTerminal())
            {
                // Compute the branch features to determine if it could account for a spine or not
                auto branchLength = branch->computeLength();
                if (branchLength < 6.0)
                {
                    // This is a spine
                    branch->setSpine();

                    // It is indeed an invalid neuronal branch
                    branch->setInvalid();

                    // Determine the parents after the removal of the spine from the tree
                    for (size_t j = 0; j < branch->parents.size(); ++j)
                    {
                        auto& parent = branch->parents[j];

                        std::vector< SkeletonBranch* > newChildren;
                        for (size_t k = 0; k < parent->children.size(); ++k)
                        {
                            // Exclude the spine branch
                            if (parent->children[k]->index == branch->index) { }
                            else { newChildren.push_back(parent->children[k]); }
                        }

                        // Clear the old children and add the new list
                        parent->children.clear();
                        parent->children.shrink_to_fit();
                        for (size_t k = 0; k < newChildren.size(); ++k)
                        {
                            parent->children.push_back(newChildren[k]);
                        }
                    }
                }
            }
        }
    }

    /// Merge the branches with a single children
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];
        if (branch->isValid())
        {
            if (branch->children.size() == 1)
            {
                // Append the nodes of the child branch to the parent branch, except the commone node
                auto& child = branch->children[0];
                for(size_t j = 1; j < child->nodes.size(); ++j)
                {
                    branch->nodes.push_back(child->nodes[j]);
                }

                // Set the new children to be that of the other one
                branch->children = branch->children[0]->children;

                // Disregard the child branch
                branch->children[0]->setInvalid();
            }
        }
    }

    /// Label the new terminals
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];
        if (branch->isValid() && branch->children.size() == 0)
        {
            std::cout << "Terminal: " << branch->index << std::endl;
            branch->setTerminal();
        }
    }
}

void NeuronSkeletonizer::_filterShortTerminalBranches()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isTerminal())
        {
            if (_branches[i]->isValid())
            {
                if (_branches[i]->nodes.size() == 2)
                {
                    // Set the branch to invalid
                    _branches[i]->setInvalid();

                    // Set the branch to be a spine
                    _branches[i]->setSpine();
                }
            }
        }
    }
}

void NeuronSkeletonizer::_detectSpines()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isTerminal())
        {
            if (_branches[i]->isValid())
            {
                if (_branches[i]->traversalCount == 1)
                {
                    if (_branches[i]->computeLength() < 6.0)
                    {
                        // Set the branch to invalid
                        _branches[i]->setInvalid();

                        // Set the branch to be a spine
                        _branches[i]->setSpine();
                    }
                }

                if (_branches[i]->traversalCount == 2)
                {
                    if (_branches[i]->computeLength() < 4.0)
                    {
                        // Set the branch to invalid
                        _branches[i]->setInvalid();

                        // Set the branch to be a spine
                        _branches[i]->setSpine();
                    }
                }
            }
        }
    }
}

void NeuronSkeletonizer::_connectBranches()
{
    // Identify the connections at the terminals of each branch
    identifyTerminalConnections(_branches);

    // Roots, terminals and others
    confirmTerminalsBranches(_branches);
}

void NeuronSkeletonizer::skeletonizeVolumeBlockByBlock(const size_t& blockSize,
                                                       const size_t& numberOverlappingVoxels,
                                                       const size_t& numberZeroVoxels)
{
    thinVolumeBlockByBlock(blockSize, numberOverlappingVoxels, numberZeroVoxels);
    constructGraph();
    segmentComponents();
}

void NeuronSkeletonizer::exportSomaProxyMesh(const std::string& prefix,
                                             const bool& formatOBJ = false,
                                             const bool& formatPLY = false,
                                             const bool& formatOFF = false,
                                             const bool& formatSTL = false)
{
    _somaProxyMesh->exportMesh(prefix + PROXY_SOMA_MESH_SUFFIX,
                               formatOBJ, formatPLY, formatOFF, formatSTL);
}

void NeuronSkeletonizer::exportSomaMesh(const std::string& prefix,
                                        const bool& formatOBJ = false,
                                        const bool& formatPLY = false,
                                        const bool& formatOFF = false,
                                        const bool& formatSTL = false)
{
    // Reconstruct the soma mesh from the proxy one, only if we are going to export it
    _reconstructSomaMeshFromProxy();
    _somaMesh->exportMesh(prefix + SOMA_MESH_SUFFIX, formatOBJ, formatPLY, formatOFF, formatSTL);
}

void NeuronSkeletonizer::collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                                         int64_t& swcIndex, int64_t branchingNodeSWCIndex)
{
    // Get a reference to the nodes of the current branch
    auto& currentBranchNodes = branch->nodes;

    for (size_t i = 1; i < currentBranchNodes.size(); ++i)
    {
        currentBranchNodes[i]->swcIndex = swcIndex;

        if (i == 1) { currentBranchNodes[i]->prevSampleSWCIndex = branchingNodeSWCIndex;}
        else { currentBranchNodes[i]->prevSampleSWCIndex= swcIndex - 1; }

        swcIndex++;
        swcNodes.push_back(currentBranchNodes[i]);
    }

    const int64_t branchingIndex = swcIndex - 1;
    for (size_t i = 0; i < branch->children.size(); ++i)
    {
        if (branch->children[i]->isValid())
        {
            collectSWCNodes(branch->children[i], swcNodes, swcIndex, branchingIndex);
        }
    }
}

SkeletonNodes NeuronSkeletonizer::constructSWCTable(const bool& resampleSkeleton)
{
    // A table, or list that contains all the nodes in order
    SkeletonNodes swcNodes;

    // A global index that will be used to correctly index the nodes
    int64_t swcIndex = 1;

    // Append the somatic mode
    _somaNode->swcIndex = swcIndex;
    _somaNode->prevSampleSWCIndex = -1;

    swcIndex++;
    swcNodes.push_back(_somaNode);

    // Resample the skeleton
    TIMER_SET;
    if (resampleSkeleton)
    {
        LOOP_STARTS("Resampling Skeleton");
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            auto& branch = _branches[i];

            // Do not resample the root sections
            if (branch->isRoot()) continue;

            // Resample only valid branches
            if (branch->isValid()) { branch->resampleAdaptively(); }
            LOOP_PROGRESS(i, _branches.size());
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Get all the root branches
    TIMER_RESET;
    LOOP_STARTS("Constructing SWC Table");
    const size_t numberBranches = _branches.size();
    for (size_t i = 0; i < numberBranches ; ++i)
    {

        auto& branch = _branches[i];
        if (branch->isRoot() && branch->isValid())
        {
            // The branching index is that of the soma
            collectSWCNodes(branch, swcNodes, swcIndex, 1);
        }
        LOOP_PROGRESS(i, numberBranches);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return swcNodes;
}

void NeuronSkeletonizer::exportSWCFile(const std::string& prefix, const bool& resampleSkeleton)
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + SWC_EXTENSION;
    LOG_STATUS("Exporting Neuron to SWC file: [ %s ]", filePath.c_str());

    auto swcNodes = constructSWCTable(resampleSkeleton);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    auto somaNode = swcNodes[0];
    stream << somaNode->swcIndex << " "
           << SWC_SOMA_STRUCT_IDENTIFIER << " "
           << somaNode->point.x() << " "
           << somaNode->point.y() << " "
           << somaNode->point.z() << " "
           << somaNode->radius << " "
           << "-1" << "\n";

    LOOP_STARTS("Writing SWC Table");
    const size_t numberSWCNodes = swcNodes.size();
    for (size_t i = 1; i < numberSWCNodes; ++i)
    {
        // TODO: Export all the branches as basal dendrites UFN
        auto swcNode = swcNodes[i];
        stream << swcNode->swcIndex << " "
               << SWC_BASAL_DENDRITE_STRUCT_IDENTIFIER << " "
               << swcNode->point.x() << " "
               << swcNode->point.y() << " "
               << swcNode->point.z() << " "
               << swcNode->radius << " "
               << swcNode->prevSampleSWCIndex << "\n";
	LOOP_PROGRESS(i, numberSWCNodes);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

void NeuronSkeletonizer::_exportSomaticBranches(const std::string& prefix) const
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + "-somatic" + BRANCHES_EXTENSION;
    LOG_STATUS("Exporting Somatic Branches: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    LOOP_STARTS("Writing Somatic Branches");
    size_t progress = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isInsideSoma())
        {
            // The @start marks a new branch in the file
            stream << "start " << _branches[i]->index << "\n";

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

        LOOP_PROGRESS(progress, _branches.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

void NeuronSkeletonizer::_exportBranches(const std::string& prefix,
                                         const SkeletonBranch::BRANCH_STATE state)
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix;
    if (state == SkeletonBranch::INVALID)
    {
        filePath += "-invalid";
    }
    else if (state == SkeletonBranch::VALID)
    {
        filePath += "-valid";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE)
    {
        filePath += "-two-sample";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_VALID)
    {
        filePath += "-two-sample-valid";
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_INVALID)
    {
        filePath += "-two-sample-invalid";
    }
    else if (state == SkeletonBranch::SOMATIC)
    {
        filePath += "-somatic";
    }
    else if (state == SkeletonBranch::SPINE)
    {
        filePath += "-spines";
    }
    else
    {
        /// NOTHING
    }

    filePath += BRANCHES_EXTENSION;

    LOG_STATUS("Exporting Neuron Branches: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    // A set of selected branches to write for debugging
    SkeletonBranches toWrite;
    if (state == SkeletonBranch::INVALID)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (!_branches[i]->isValid())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->nodes.size() == 2)
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_VALID)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->nodes.size() == 2 && _branches[i]->isValid())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::TWO_SAMPLE_AND_INVALID)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->nodes.size() == 2 && !_branches[i]->isValid())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::SOMATIC)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->isInsideSoma())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::SPINE)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->isSpine())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else if (state == SkeletonBranch::VALID)
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            if (_branches[i]->isValid())
            {
                toWrite.push_back(_branches[i]);
            }
        }
    }
    else
    {
        for (size_t i = 0; i < _branches.size(); ++i)
        {
            toWrite.push_back(_branches[i]);
        }
    }

    LOOP_STARTS("Writing Branches");
    size_t progress = 0;
    for (size_t i = 0; i < toWrite.size(); ++i)
    {
        // The @start marks a new branch in the file
        stream << "start " << toWrite[i]->index << "\n";

        for (auto& node: toWrite[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        // The @end marks the terminal sample of a branch
        stream << "end\n";

        LOOP_PROGRESS(progress, toWrite.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

void NeuronSkeletonizer::_exportSpineLocations(const std::string& prefix) const
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + "-spine-locations" + NODES_EXTENSION;
    LOG_STATUS("Exporting Spines Locations: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    LOOP_STARTS("Writing Spines Locations");
    size_t progress = 0;
    for (size_t i = 0; i < _spineRoots.size(); ++i)
    {
        auto node = _spineRoots[i]->nodes[0];

        stream << node->point.x() << " "
               << node->point.y() << " "
               << node->point.z() << "\n";

        LOOP_PROGRESS(progress, _spineRoots.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}


void NeuronSkeletonizer::_handleSpines(SkeletonBranch* root) const
{
    // From the spine morphology create a mesh
    // Create a new class called SpineMorphology
    SpineMorphology* spineMorphology = new SpineMorphology(root);
    std::stringstream prefixs;
    prefixs << "/home/abdellah/spines/models/spine_" << root->index;
    spineMorphology->exportBranches(prefixs.str());

    auto mesh = spineMorphology->reconstructMesh(20, true);
    std::stringstream prefixss;
    prefixss << "/home/abdellah/spines/models/spine_" << root->index;
    mesh->exportMesh(prefixss.str(), true);







    // Voxelize the neuro mesh at the selected region

    // Get the cloud of the volume

    // Map the spine mesh onto the cloud

    // Re-skeletonize the mesh and extract a nice skeleton !

}



void NeuronSkeletonizer::_exportSpineExtents(const std::string& prefix) const
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + "-spine-extents" + NODES_EXTENSION;
    LOG_STATUS("Exporting Spines Extents: [ %s ]", filePath.c_str());

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    LOOP_STARTS("Writing Spines Extents");
    size_t progress = 0;
    for (size_t i = 0; i < _spineRoots.size(); ++i)
    {
        Vector3f pMin, pMax, bounds, center;
        getLogicalTreeBoundingBox(_spineRoots[i], pMin, pMax, bounds, center);

        stream << center.x() << " "
               << center.y() << " "
               << center.z() << " "
               << bounds.x() << " "
               << bounds.y() << " "
               << bounds.z() << "\n";

        LOOP_PROGRESS(progress, _spineRoots.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();
}

}
