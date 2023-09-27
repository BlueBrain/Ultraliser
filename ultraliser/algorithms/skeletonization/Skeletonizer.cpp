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

namespace Ultraliser
{
Skeletonizer::Skeletonizer(Volume* volume, const Mesh *mesh)
    : _volume(volume)
    , _mesh(mesh)
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

    _computeShellPoints();
}

void Skeletonizer::applyVolumeThinningToVolume(Volume* volume, const bool& displayProgress)
{
    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();

    if (displayProgress)
    {
        // Parameters to calculate the loop progress
        size_t initialNumberVoxelsToBeDeleted = 0;
        size_t loopCounter = 0;

        TIMER_SET;
        LOG_STATUS("Thinning Volume");
        LOOP_STARTS("Thinning Loop");
        LOOP_PROGRESS(0, 100);
        while(1)
        {
            size_t numberDeletedVoxels = volume->deleteCandidateVoxelsParallel(thinningKernel);

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
    else
    {
        // Parameters to calculate the loop progress
        size_t initialNumberVoxelsToBeDeleted = 0;
        size_t loopCounter = 0;
        while(1)
        {
            size_t numberDeletedVoxels = volume->deleteCandidateVoxelsParallel(thinningKernel);
           if (numberDeletedVoxels == 0)
               break;
        }
    }
}

void Skeletonizer::applyVolumeThinningWithDomainDecomposition()
{
    // Start the timer
    TIMER_SET;

    // Copy the input volume into a reference volume
    // The reference volume will be used to retrieve the new bricks, and the _volume will be
    // used to write the skeletonization result
    Volume* referenceVolume = new Volume(_volume->getWidth(),
                                         _volume->getHeight(),
                                         _volume->getDepth());

    const size_t subdivisions = 4;
    const size_t overlappingVoxels = 5;
    const size_t numberZeroVoxels = 2;

    Ranges xRanges = Range::decomposeToRanges(int64_t(0), _volume->getWidth() - 1, subdivisions);
    Ranges yRanges = Range::decomposeToRanges(int64_t(0), _volume->getHeight() - 1, subdivisions);
    Ranges zRanges = Range::decomposeToRanges(int64_t(0), _volume->getDepth() - 1, subdivisions);

    // Add the overlaps
    Ranges xRangesOverlapping = Range::addTwoSidedOverlaps(xRanges, overlappingVoxels);
    Ranges yRangesOverlapping = Range::addTwoSidedOverlaps(yRanges, overlappingVoxels);
    Ranges zRangesOverlapping = Range::addTwoSidedOverlaps(zRanges, overlappingVoxels);

    LOG_STATUS("Skeletonizing Volume Bricks");
    LOOP_STARTS("Skeletonization");
    int64_t progress = 0;
    for (size_t i = 0; i < xRanges.size(); ++i)
    {
        LOOP_PROGRESS(progress, xRanges.size());
        for (size_t j = 0; j < yRanges.size(); ++j)
        {
            // OMP_PARALLEL_FOR
            for (size_t k = 0; k < zRanges.size(); ++k)
            {
                // Extract the brick from the volume
                auto brick = _volume->extractBoundedBrickFromVolume(
                            xRangesOverlapping[i].i1, xRangesOverlapping[i].i2,
                            yRangesOverlapping[j].i1, yRangesOverlapping[j].i2,
                            zRangesOverlapping[k].i1, zRangesOverlapping[k].i2,
                            numberZeroVoxels, false);

                // Skeletonize the brick
                applyVolumeThinningToVolume(brick, false);

                size_t xOverlapping, yOverlapping, zOverlapping = 0;
                if (i > 0) xOverlapping = overlappingVoxels; else xOverlapping = 0;
                if (j > 0) yOverlapping = overlappingVoxels; else yOverlapping = 0;
                if (k > 0) zOverlapping = overlappingVoxels; else zOverlapping = 0;

                referenceVolume->insertOverlappingBoundedBrickToVolume(
                            brick,
                            xRanges[i].i1, xRanges[i].i2,
                            yRanges[j].i1, yRanges[j].i2,
                            zRanges[k].i1, zRanges[k].i2,
                            xOverlapping, yOverlapping, zOverlapping, numberZeroVoxels,
                            false);

                brick->~Volume();
            }
        }

        progress++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    _volume->insertBrickToVolume(referenceVolume,
                                 0, referenceVolume->getWidth() - 1,
                                 0, referenceVolume->getHeight() - 1,
                                 0, referenceVolume->getDepth() - 1);
    referenceVolume->~Volume();
}

void Skeletonizer::thinVolumeBlockByBlock(const size_t& blockSize,
                                          const size_t& numberOverlappingVoxels,
                                          const size_t& numberZeroVoxels)
{
    // Start the timer
    TIMER_SET;

    // Copy the input volume into a reference volume
    // The reference volume will be used to retrieve the new bricks, and the _volume will be
    // used to write the skeletonization result
    Volume* referenceVolume = new Volume(_volume->getWidth(),
                                         _volume->getHeight(),
                                         _volume->getDepth());

    // Initially, the range of the volume is identified
    Range xRange(0, _volume->getWidth() - 1);
    Range yRange(0, _volume->getHeight() - 1);
    Range zRange(0, _volume->getDepth() - 1);

    // Decompose the range into ranges, based on the blockSize
    Ranges xRanges = xRange.decomposeToBlocks(blockSize);
    Ranges yRanges = yRange.decomposeToBlocks(blockSize);
    Ranges zRanges = zRange.decomposeToBlocks(blockSize);

    for (const auto& range: xRanges)
        range.printRange();
    for (const auto& range: yRanges)
        range.printRange();
    for (const auto& range: zRanges)
        range.printRange();


    // Add the overlaps
    Ranges xRangesOverlapping = Range::addTwoSidedOverlaps(xRanges, numberOverlappingVoxels);
    Ranges yRangesOverlapping = Range::addTwoSidedOverlaps(yRanges, numberOverlappingVoxels);
    Ranges zRangesOverlapping = Range::addTwoSidedOverlaps(zRanges, numberOverlappingVoxels);

    LOG_STATUS("Skeletonizing Volume Bricks");
    LOOP_STARTS("Skeletonization");
    int64_t progress = 0;
    for (size_t i = 0; i < xRanges.size(); ++i)
    {
        LOOP_PROGRESS(progress, xRanges.size());
        for (size_t j = 0; j < yRanges.size(); ++j)
        {
            // OMP_PARALLEL_FOR
            for (size_t k = 0; k < zRanges.size(); ++k)
            {
                // Extract the brick from the volume
                auto brick = _volume->extractBoundedBrickFromVolume(
                            xRangesOverlapping[i].i1, xRangesOverlapping[i].i2,
                            yRangesOverlapping[j].i1, yRangesOverlapping[j].i2,
                            zRangesOverlapping[k].i1, zRangesOverlapping[k].i2,
                            numberZeroVoxels, false);

                // Skeletonize the brick
                applyVolumeThinningToVolume(brick, false);

                size_t xOverlapping, yOverlapping, zOverlapping = 0;
                if (i > 0) xOverlapping = numberZeroVoxels; else xOverlapping = 0;
                if (j > 0) yOverlapping = numberZeroVoxels; else yOverlapping = 0;
                if (k > 0) zOverlapping = numberZeroVoxels; else zOverlapping = 0;

                referenceVolume->insertOverlappingBoundedBrickToVolume(
                            brick,
                            xRanges[i].i1, xRanges[i].i2,
                            yRanges[j].i1, yRanges[j].i2,
                            zRanges[k].i1, zRanges[k].i2,
                            xOverlapping, yOverlapping, zOverlapping, numberZeroVoxels,
                            false);

                brick->~Volume();
            }
        }

        progress++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    _volume->insertBrickToVolume(referenceVolume,
                                 0, referenceVolume->getWidth() - 1,
                                 0, referenceVolume->getHeight() - 1,
                                 0, referenceVolume->getDepth() - 1);
    referenceVolume->~Volume();
}

void Skeletonizer::applyVolumeThinning()
{
    TIMER_SET;
    LOG_STATUS("Volume Thinning");

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

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

    // TODO: Adjust the voxel slight shift
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

    std::fstream stream;
    std::string filePath = "/data/microns-explorer-dataset/Meshes-Input-MICrONS/skeletonization-spines/shell.txt";
    stream.open(filePath, std::ios::out);
    for (size_t i = 1; i < _shellPoints.size(); ++i)
    {
        stream << _shellPoints[i].x() << " "
               << _shellPoints[i].y() << " "
               << _shellPoints[i].z() << "\n";
    }
    stream.close();
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsParallel()
{
    LOG_STATUS("Mapping Voxels to Nodes * Parallel");
    TIMER_SET;

    struct FilledVoxel
    {
        size_t i, j, k;
        size_t voxelIndex;

        FilledVoxel(size_t ii, size_t jj, size_t kk, size_t vIndex)
        { i = ii; j = jj; k = kk; voxelIndex = vIndex; }
    };

    typedef std::vector< FilledVoxel* > FilledVoxels;


    std::vector< FilledVoxels > allFilledVoxels;
    allFilledVoxels.resize(_volume->getWidth());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        auto& perSliceFilledVoxels = allFilledVoxels[i];

        for (size_t j = 0; j < _volume->getHeight(); ++j)
        {
            for (size_t k = 0; k < _volume->getDepth(); ++k)
            {
                // If the voxel is filled
                if (_volume->isFilled(i, j, k))
                {
                    perSliceFilledVoxels.push_back(new FilledVoxel(i, j, k, _volume->mapTo1DIndexWithoutBoundCheck(i, j, k)));
                }
            }
        }
    }

    FilledVoxels filledVoxels;
    for (size_t i = 0; i < allFilledVoxels.size(); ++i)
    {
        if (allFilledVoxels[i].size() > 0)
        {
            filledVoxels.insert(filledVoxels.end(), allFilledVoxels[i].begin(), allFilledVoxels[i].end());
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
        // Mapper from voxel to node indices
        indicesMapper.insert(std::pair< size_t, size_t >(filledVoxels[i]->voxelIndex, i));
    }


    // Resize the nodes
    _nodes.resize(filledVoxels.size());

    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < filledVoxels.size(); ++n)
    {
        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, filledVoxels.size());
        PROGRESS_UPDATE;

        const size_t i = filledVoxels[n]->i;
        const size_t j = filledVoxels[n]->j;
        const size_t k = filledVoxels[n]->k;
        const size_t voxelIndex = filledVoxels[n]->voxelIndex;

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
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    filledVoxels.clear();
    filledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxels()
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
    size_t progress = 0;
    for (size_t i = 0; i < _volume->getWidth(); ++i)
    {
        LOOP_PROGRESS(progress, _volume->getWidth());

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

        progress++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Resize the nodes
    _nodes.resize(indicesFilledVoxels.size());

    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < indicesFilledVoxels.size(); ++n)
    {
        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, indicesFilledVoxels.size());
        PROGRESS_UPDATE;

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
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    indicesFilledVoxels.clear();
    indicesFilledVoxels.shrink_to_fit();

    return indicesMapper;
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
        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, _nodes.size());
        PROGRESS_UPDATE;

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
    // Construct the graph and connect the nodes
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        // LOOP_PROGRESS(loopProgress, _nodes.size());

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
    }
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
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxelsParallel();

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

                            // TODO: Remove
                            // branch->setValid();

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
