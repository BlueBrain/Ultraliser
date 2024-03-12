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

#include "SpineSkeletonizer.h"

namespace Ultraliser
{

SpineSkeletonizer::SpineSkeletonizer(Volume* spineVolume,
                                     const bool useAcceleration,
                                     const bool debugSkeleton,
                                     const std::string debuggingPrefix)
    : Skeletonizer(spineVolume, useAcceleration, debugSkeleton, debuggingPrefix)
{
    /// EMPTY CONSTRUCTOR
}

void SpineSkeletonizer::run(const bool verbose)
{
    // Initialize
    initialize(verbose);

    skeletonizeVolumeToCenterLines();

    /// Extract the nodes of the skeleton from the center-line "thinned" voxels and return a
    /// mapper that maps the indices of the voxels in the volume and the nodes in the skeleton
    auto indicesMapper = _extractNodesFromVoxels();

    /// Connect the nodes of the skeleton to construct its edges. This operation will not connect
    /// any gaps, it will just connect the nodes extracted from the voxels.
    _connectNodesToBuildEdges(indicesMapper);

    /// Inflate the nodes, i.e. adjust their radii
    _inflateNodes();

    /// Reconstruct the sections "or branches" from the nodes using the edges data
    _buildBranchesFromNodes(_nodes);


    _exportBranches(_debuggingPrefix);


    // Identify the connections at the terminals of each branch
    // identifyTerminalConnections(_branches);

    // Roots, terminals and others
    // confirmTerminalsBranches(_branches);
}

void SpineSkeletonizer::segmentComponents()
{

}


void SpineSkeletonizer::_exportBranches(const std::string& prefix)
{

    // Construct the file path
    std::string filePath = prefix + BRANCHES_EXTENSION;

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    size_t progress = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
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

        ++progress;
    }

    // Close the file
    stream.close();
}

void SpineSkeletonizer::_collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                                         int64_t& swcIndex, int64_t branchingNodeSWCIndex, const bool)
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
            _collectSWCNodes(branch->children[i], swcNodes, swcIndex, branchingIndex);
        }
    }
}

SkeletonNodes SpineSkeletonizer::_constructSWCTable(const bool& resampleSkeleton, const bool)
{
    // A table, or list that contains all the nodes in order
    SkeletonNodes swcNodes;

    // A global index that will be used to correctly index the nodes
    int64_t swcIndex = 1;

    // Fake soma node
    SkeletonNode* fakeSomaNode = new SkeletonNode();

    // Append the somatic mode
    fakeSomaNode->swcIndex = swcIndex;
    fakeSomaNode->prevSampleSWCIndex = -1;

    swcIndex++;
    swcNodes.push_back(fakeSomaNode);

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
            _collectSWCNodes(branch, swcNodes, swcIndex, 1);
        }
        LOOP_PROGRESS(i, numberBranches);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return swcNodes;
}

void SpineSkeletonizer::exportSWCFile(const std::string& prefix, const bool& resampleSkeleton, const bool)
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + SWC_EXTENSION;
    LOG_STATUS("Exporting Spine to SWC file: [ %s ]", filePath.c_str());

    auto swcNodes = _constructSWCTable(resampleSkeleton);

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

}
