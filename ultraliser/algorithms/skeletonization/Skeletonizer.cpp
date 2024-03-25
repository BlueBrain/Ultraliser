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

#include <common/Defines.h>
#include "Skeletonizer.h"
#include "SkeletonizerUtils.h"
#include <algorithms/skeletonization/thinning/Neighbors.hh>
#include <algorithms/utilities/KdTree.h>
#include <math/Vector.h>
#include <data/meshes/simple/TriangleOperations.h>
#include <utilities/Range.h>
#include <data/volumes/voxels/NodeVoxel.h>

namespace Ultraliser
{
Skeletonizer::Skeletonizer(Volume* volume,
                           const bool &useAcceleration,
                           const bool &debugSkeleton,
                           const std::string debuggingPrefix)
    : _volume(volume)
    , _mesh(nullptr)
    , _useAcceleration(useAcceleration)
    , _debugSkeleton(debugSkeleton)
    , _debuggingPrefix(debuggingPrefix)
    , _debug(_debugSkeleton && _debuggingPrefix != NONE)
{
    /// NOTE: The mesh is assigned a nullptr, until further notice

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

Skeletonizer::Skeletonizer(Mesh* mesh,
                           const VoxelizationOptions& options,
                           const bool &useAcceleration,
                           const bool &debugSkeleton,
                           const std::string debuggingPrefix)
    : _volume(nullptr)
    , _mesh(mesh)
    , _voxelizationOptions(options)
    , _useAcceleration(useAcceleration)
    , _debugSkeleton(debugSkeleton)
    , _debuggingPrefix(debuggingPrefix)
    , _debug(_debugSkeleton && _debuggingPrefix != NONE)
{
    /// NOTE: The volume is assigned a nullptr, until further notice

    // Compute the mesh bounding box
    mesh->computeBoundingBox(_pMinMesh, _pMaxMesh);
    _boundsMesh = _pMaxMesh - _pMinMesh;
    _centerMesh = _pMinMesh + 0.5 * _boundsMesh;

    /// NOTE: Compute the volume bounds after the generation of the volume
}

void Skeletonizer::_computeVolumeFromMesh()
{
    // If the mesh is a nullptr, then return, there is nothing to compute
    if (_mesh == nullptr)
    {
        LOG_ERROR("Skeletonizer::_computeVolumeFromMesh(): An empty mesh is given to voxelize!");
        return;
    }

    // The volume must be a nullptr to be able to compute it
    if (_volume == nullptr)
    {
        // Create the volume extent
        _volume = new Volume(_pMinMesh, _pMaxMesh,
                             _voxelizationOptions.volumeResolution,
                             _voxelizationOptions.edgeGapPrecentage,
                             _voxelizationOptions.volumeType,
                             _voxelizationOptions.verbose);

        // Apply surface and solid voxelization to the input neuron mesh
        _volume->surfaceVoxelization(_mesh, _voxelizationOptions.verbose, false, 1.0);
        _volume->solidVoxelization(_voxelizationOptions.voxelizationAxis,
                                   _voxelizationOptions.verbose);

        // Remove the border voxels that span less than half the voxel
        auto bordeVoxels = _volume->searchForBorderVoxels(_voxelizationOptions.verbose);
        for (size_t i = 0; i < bordeVoxels.size(); ++i)
        {
            for (size_t j = 0; j < bordeVoxels[i].size(); ++j)
            {
                _volume->clear(bordeVoxels[i][j].x(), bordeVoxels[i][j].y(), bordeVoxels[i][j].z());
            }
            bordeVoxels[i].clear();
        }
        bordeVoxels.clear();
        _volume->surfaceVoxelization(_mesh, _voxelizationOptions.verbose, false, 0.5);
    }
}

void Skeletonizer::initialize(const bool verbose)
{
    TIMER_SET;
    VERBOSE_LOG(LOG_TITLE("Ultraliser Skeletonization"), verbose);
    VERBOSE_LOG(LOG_SUCCESS("Voxel Size [%f] μm", _volume->getVoxelSize()), verbose);
    VERBOSE_LOG(LOG_STATUS("Initialization - Building Structures"), verbose);

    // Compute the shell points either natively or by using the acceleration structures
    if (_useAcceleration)
    {
        // Build the ThinningVoxels acceleration structure from the input solid volume
        // NOTE: We do not rebuild the ThinningVoxels structure!
        auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

        // Compute the surface shell from the pre-built ThinningVoxels structure
        _computeShellPointsUsingAcceleration(thinningVoxels, verbose);
    }
    else { _computeShellPoints(verbose); }

    VERBOSE_LOG(LOG_STATUS_IMPORTANT("Initialization Stats."), verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

Sections Skeletonizer::getValidSections() const
{
    Sections sections;
    size_t sectionIndex = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isValid())
        {
            Section* section = new Section(sectionIndex++);
            for (size_t j = 0; j < _branches[i]->nodes.size(); ++j)
            {
                auto node = _branches[i]->nodes[j];
                section->addSample(new Sample(node->point, node->radius, j));
            }
            sections.push_back(section);
        }
    }

    return sections;
}

void Skeletonizer::_scaleShellPoints(const bool verbose)
{
    // Initialize the timer
    TIMER_SET;

    // TODO: Adjust the voxel slight shift
    // Adjust the locations of the shell points taking into consideration the mesh coordinates
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Mapping Shell Points"), verbose);
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

        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _shellPoints.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_computeShellPointsUsingAcceleration(ThinningVoxelsUI16List &thinningVoxels,
                                                        const bool verbose)
{
    // Initialize the timer
    TIMER_SET;

    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Computing Shell Points *"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        auto& voxel = thinningVoxels[i];
        if (_volume->isBorderVoxel(voxel->x, voxel->y, voxel->z))
        {
            voxel->border = true;
        }

        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, thinningVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

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
    _scaleShellPoints(verbose);
}

void Skeletonizer::_computeShellPoints(const bool verbose)
{
    // Initialize the time
    TIMER_SET;

    // Search for the border voxels (the shell voxels) of the volume
    std::vector< std::vector< Vec3ui_64 > > perSlice = _volume->searchForBorderVoxels(verbose);

    // Concatinate the points in a single list
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Computing Shell Points"), verbose);
    for (size_t i = 0; i < perSlice.size(); ++i)
    {
        for (size_t j = 0; j < perSlice[i].size(); ++j)
        {
            const auto voxel = perSlice[i][j];
            _shellPoints.push_back(Vector3f(voxel.x(), voxel.y(), voxel.z()));
        }
        perSlice[i].clear();
        perSlice[i].shrink_to_fit();
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, perSlice.size()), verbose);
    }
    perSlice.clear();
    perSlice.shrink_to_fit();
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Scale the shell points to match the extent of the input data
    _scaleShellPoints(verbose);
}

void Skeletonizer::_applyVolumeThinning(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Volume Thinning"), verbose);

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Thinning Loop"), verbose);
    VERBOSE_LOG(LOOP_PROGRESS(0, 100), verbose);
    while(1)
    {
        size_t numberDeletedVoxels = _volume->deleteCandidateVoxelsParallel(thinningKernel);

        // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        VERBOSE_LOG(LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                      initialNumberVoxelsToBeDeleted), verbose);

        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_applyVolumeThinningUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Volume Thinning *"), verbose);

    // The thinning kernel that will be used to thin the volume
    std::unique_ptr< Thinning6Iterations > thinningKernel =
            std::make_unique< Thinning6Iterations >();

    // Parameters to calculate the loop progress
    size_t initialNumberVoxelsToBeDeleted = 0;
    size_t loopCounter = 0;

    auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Thinning Loop"), verbose);
    VERBOSE_LOG(LOOP_PROGRESS(0, 100), verbose);
    while(1)
    {
        // Delete the border voxels based on the ThinningVoxels acceleration structure
        size_t numberDeletedVoxels = _volume->deleteBorderVoxelsUsingThinningVoxels(
                    thinningKernel, thinningVoxels);

        // Updating the progess bar
        if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
        VERBOSE_LOG(LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                                  initialNumberVoxelsToBeDeleted), verbose);

        // No more voxels to be deleted
        if (numberDeletedVoxels == 0)
            break;

        loopCounter++;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::skeletonizeVolumeToCenterLines(const bool verbose)
{
    if (_useAcceleration)
        _applyVolumeThinningUsingAcceleration(verbose);
    else
        _applyVolumeThinning(verbose);
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxels(const bool verbose)
{
    if (_useAcceleration)
        return _extractNodesFromVoxelsUsingAcceleration(verbose);
    else
        return _extractNodesFromVoxelsUsingSlicing(verbose);
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsNaive(const bool verbose)
{
    LOG_STATUS("Mapping Voxels to Nodes");

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;

    // A list of filled voxels to compute the elements in parallel
    std::vector< Vec4ui_64 > indicesFilledVoxels;

    // Search the filled voxels in the volume
    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Detecting Filled Voxels"), verbose);
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
        VERBOSE_LOG(LOOP_PROGRESS(i, _volume->getWidth()), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Resize the nodes
    _nodes.resize(indicesFilledVoxels.size());

    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Building Graph Nodes"), verbose);
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
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, indicesFilledVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Clear the auxiliary lists
    indicesFilledVoxels.clear();
    indicesFilledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingSlicing(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Mapping Voxels to Nodes"), verbose);

    /**
     * @brief The FilledVoxel struct
     */
    struct FilledVoxel
    {
        size_t i, j, k, idx;

        FilledVoxel(size_t ii, size_t jj, size_t kk, size_t index)
        { i = ii; j = jj; k = kk; idx = index; }
    };

    /**
     * @brief FilledVoxels
     */
    typedef std::vector< FilledVoxel > FilledVoxels;

    // Make a per-slice list
    std::vector< FilledVoxels > allFilledVoxels;
    allFilledVoxels.resize(_volume->getWidth());

    TIMER_SET;
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
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, filledVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    filledVoxels.clear();
    filledVoxels.shrink_to_fit();

    return indicesMapper;
}

std::map< size_t, size_t > Skeletonizer::_extractNodesFromVoxelsUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Mapping Voxels to Nodes *"), verbose);

    // Get all the center-line voxels from the volume
    auto thinningVoxels = _volume->getThinningVoxelsList(false, verbose);

    // Construct the NodeVoxels list
    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Constructing Node Voxels"), verbose);
    PROGRESS_SET;
    NodeVoxelsUI16 nodeVoxels;
    for (size_t i = 0; i < thinningVoxels.size(); ++i)
    {
        // Reference to the voxel
        auto& voxel = thinningVoxels[i];

        // The inactive voxels have been deactivated during the thinning
        if (voxel->active)
        {
            // Create a corresponding node to the voxel
            NodeVoxelUI16 nodeVoxel;

            // Get the location based on the voxel XYZ coordinates
            nodeVoxel.x = voxel->x; nodeVoxel.y = voxel->y; nodeVoxel.z = voxel->z;

            // Update the voxel index (used later to map the voxels to nodes)
            nodeVoxel.voxelIndex = _volume->mapTo1DIndexWithoutBoundCheck(
                        voxel->x, voxel->y, voxel->z);

            // Add th node voxel to the list
            nodeVoxels.push_back(nodeVoxel);
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, thinningVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Every constructed node must have an identifier, or index.
    size_t nodeIndex = 0;

    // Map to associate between the indices of the voxels and the nodes in the graph
    std::map< size_t, size_t > indicesMapper;
    for (size_t i = 0; i < nodeVoxels.size(); ++i)
    {
        indicesMapper.insert(std::pair< size_t, size_t >(nodeVoxels[i].voxelIndex, i));
    }

    // Resize the nodes to the corresponding size of the NodeVoxels list
    _nodes.resize(nodeVoxels.size());

    PROGRESS_RESET;
    VERBOSE_LOG(LOOP_STARTS("Building Graph Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t n = 0; n < nodeVoxels.size(); ++n)
    {
        // Get a point representing the center of the voxel (in the volume)
        Vector3f voxelPosition(nodeVoxels[n].x * 1.f,
                               nodeVoxels[n].y * 1.f,
                               nodeVoxels[n].z * 1.f);

        // Get a point in the same coordinate space of the mesh
        /// TODO: Adjust the center of the node based on the actual center of the voxel
        Vector3f nodePosition(voxelPosition);

        // Adjust the location based on the dimensions of the input data
        nodePosition -= _centerVolume;
        nodePosition.x() *= _scaleFactor.x();
        nodePosition.y() *= _scaleFactor.y();
        nodePosition.z() *= _scaleFactor.z();
        nodePosition += _centerMesh;

        // Add the node to the nodes list
        _nodes[n] = new SkeletonNode(n, nodeVoxels[n].voxelIndex, nodePosition, voxelPosition);

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, nodeVoxels.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    return indicesMapper;
}

void Skeletonizer::_inflateNodes(const bool verbose)
{
    if (_useAcceleration)
        _inflateNodesUsingAcceleration(verbose);
    else
        _inflateNodesNatively(verbose);
}

void Skeletonizer::_inflateNodesUsingAcceleration(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Inflating Graph Nodes - Mapping to Surface"), verbose);

    auto kdtree = KdTree::from(_shellPoints);

    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Inflating Nodes"), verbose);
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto &node = *_nodes[i];

        auto nearestPoint = kdtree.findNearestPoint(node.point);
        auto minimumDistance = nearestPoint.distance;

        // TODO: Make some logic to detect the actual radius based on the voxel size
        if (minimumDistance > _volume->getVoxelSize())
        {
            node.radius = minimumDistance * 1.2;
        }
        else
        {
            node.radius = _volume->getVoxelSize() * 0.5;
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_inflateNodesNatively(const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Inflating Graph Nodes - Mapping to Surface"), verbose);

    TIMER_SET;
    PROGRESS_SET;
    VERBOSE_LOG(LOOP_STARTS("Inflating Nodes"), verbose);
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
        }
        else
        {
            _nodes[i]->radius = 0.1;
        }

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_connectNodesToBuildEdges(const std::map< size_t, size_t >& indicesMapper,
                                             const bool verbose)
{
    VERBOSE_LOG(LOG_STATUS("Connecting Graph Nodes"), verbose);

    size_t numberEdges;
    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Building & Linking Edges"), verbose);
    PROGRESS_SET;
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
                // Increment the number of conected edges along this node
                connectedEdges++;

                // Find the index of the voxel
                const auto& vIndex = _volume->mapTo1DIndexWithoutBoundCheck(idx, idy, idz);

                // Find the corresponding index of the node to access the node from the nodes list
                const auto& nIndex = indicesMapper.find(vIndex)->second;

                // Add the node to the edgeNodes, only to be able to access it later
                SkeletonNode* edgeNode = _nodes[nIndex];
                node->edgeNodes.push_back(edgeNode);

                // Construct the edge
                SkeletonEdge* edge = new SkeletonEdge(numberEdges, node, edgeNode);

                // Add the constructed edge to the list
                _edges.push_back(edge);

                // Increment the number of edges
                numberEdges++;
            }
        }

        if (connectedEdges == 1)
            node->terminal = true;

        if (connectedEdges > 2)
            node->branching = true;

        // Update the progress bar
        VERBOSE_LOG(LOOP_PROGRESS(PROGRESS, _nodes.size()), verbose);
        PROGRESS_UPDATE;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);
}

void Skeletonizer::_removeTriangleLoops(const bool verbose)
{
    TIMER_SET;
    VERBOSE_LOG(LOG_STATUS("Removing Triangle Loops"), verbose);

    const size_t currentNodesSize = _nodes.size();
    for (size_t i = 0; i < currentNodesSize; ++i)
    {
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

        VERBOSE_LOG(LOOP_PROGRESS(i, currentNodesSize), verbose);
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Reset the visiting state
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i) { _nodes[i]->visited = false; }
}

void Skeletonizer::constructGraph(const bool verbose)
{
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxels();

    // Assign accurate radii to the nodes of the graph, i.e. inflate the nodes
    _inflateNodes(verbose);

    // Connect the nodes to construct the edges of the graph
    _connectNodesToBuildEdges(indicesMapper);

    // Remove the triangular configurations
    _removeTriangleLoops(verbose);

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

void Skeletonizer::skeletonizeVolumeBlockByBlock(const size_t& blockSize,
                                                 const size_t& numberOverlappingVoxels,
                                                 const size_t& numberZeroVoxels,
                                                 const bool verbose)
{
    thinVolumeBlockByBlock(blockSize, numberOverlappingVoxels, numberZeroVoxels);
    constructGraph(verbose);
    segmentComponents();
}

std::vector< Vector3f > Skeletonizer::getShellPoints()
{
    return _shellPoints;
}

void Skeletonizer::_exportGraphNodes(const std::string prefix, const bool verbose)
{
    // Construct the file path
    std::string filePath = prefix + NODES_EXTENSION;
    VERBOSE_LOG(LOG_STATUS("Exporting Nodes : [ %s ]", filePath.c_str()), verbose);

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    TIMER_SET;
    VERBOSE_LOG(LOOP_STARTS("Writing Nodes"), verbose);
    size_t progress = 0;
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];
        stream << node->point.x() << " "
               << node->point.y() << " "
               << node->point.z() << " "
               << node->radius << NEW_LINE;

        VERBOSE_LOG(LOOP_PROGRESS(progress, _nodes.size()), verbose);
        ++progress;
    }
    VERBOSE_LOG(LOOP_DONE, verbose);
    VERBOSE_LOG(LOG_STATS(GET_TIME_SECONDS), verbose);

    // Close the file
    stream.close();
}

}
