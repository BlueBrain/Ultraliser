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

#pragma once

#include <algorithms/skeletonization/Skeletonizer.h>
#include <algorithms/skeletonization/SkeletonWeightedEdge.hh>
#include <algorithms/skeletonization/graphs/GraphNode.h>
#include <algorithms/skeletonization/graphs/GraphBranch.h>
#include <algorithms/skeletonization/graphs/ShortestPathFinder.h>

namespace Ultraliser
{

/**
 * @brief The NeuronSkeletonizer class
 */
class NeuronSkeletonizer : public Skeletonizer
{
public:

    /**
     * @brief NeuronSkeletonizer
     * @param volume
     */
    NeuronSkeletonizer(Volume *volume, const bool &useAcceleration = true);
    ~NeuronSkeletonizer();

    /**
     * @brief skeletonizeVolumeToCenterLines
     */
    void skeletonizeVolumeToCenterLines() override;

    /**
     * @brief constructGraph
     */
    void constructGraph() override;

    /**
     * @brief segmentComponents
     */
    void segmentComponents() override;

    /**
     * @brief skeletonizeVolumeBlockByBlock
     * @param blockSize
     * @param numberOverlappingVoxels
     * @param numberZeroVoxels
     */
    void skeletonizeVolumeBlockByBlock(const size_t& blockSize = 512,
                                       const size_t& numberOverlappingVoxels = 25,
                                       const size_t& numberZeroVoxels = 5) override;

    /**
      * @brief exportIndividualBranches
      * Exports the individual branches of the skeleton, irrespective of their connectivity.
      * This function is used for visual debugging.
      * @param prefix
      * File prefix.
      */
    void exportIndividualBranches(const std::string& prefix) const;

    /**
      * @brief exportSWCFile
      * Export the resulting skeleton into an SWC file.
      * This function is called after the segmentation of all the components from the skeleton.
      * @param prefix
      * File prefix.
      */
    void exportSWCFile(const std::string& prefix);

    /**
     * @brief constructSWCTable
     * To facilitate exporting the skeleton into an SWC file, we construct a linear list of
     * all the SWC in order and then export this table in the @exportSWCFile function.
     * @return
     * A list of all the skeleton nodes, constructed in order, and ready for being exported into
     * SWC file.
     */
    SkeletonNodes constructSWCTable();

    /**
     * @brief collectSWCNodes
     * Traverses the skeleton using depth first search and builds a list of nodes that can be
     * exported into SWC files.
     * @param branch
     * @param swcNodes
     * @param swcIndex
     * @param branchingNodeSWCIndex
     */
    void collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                         int64_t &swcIndex, int64_t branchingNodeSWCIndex);

    /**
     * @brief exportSomaMesh
     * Exports the mesh of the segmented soma.
     * @param prefix
     * File prefix.
     * @param formatOBJ
     * If this flag is set, an OBJ mesh will be exported.
     * @param formatPLY
     * If this flag is set, a PLY mesh will be exported.
     * @param formatOFF
     * If this flag is set, an OFF mesh will be exported.
     * @param formatSTL
     * If this flag is set, an STL mesh will be exported.
     */
    void exportSomaMesh(const std::string& prefix,
                        const bool &formatOBJ,
                        const bool &formatPLY,
                        const bool &formatOFF,
                        const bool &formatSTL);

    /**
     * @brief getSomaMesh
     * Returns a pointer to the segmented soma mesh.
     * @return
     * Returns a pointer to the segmented soma mesh.
     */
    Mesh* getSomaMesh() const { return _somaMesh; }

private:

    /**
     * @brief _addSomaNode
     * Adds a dedicated node for the soma. This node is updated later after the soma is segmented
     * from the mesh and all the arbors are detected.
     */
    void _addSomaNode();

    /**
     * @brief _segmentSomaMesh
     */
    void _segmentSomaMesh();

    /**
     * @brief _identifySomaticNodes
     */
    void _identifySomaticNodes();

    /**
     * @brief _removeBranchesInsideSoma
     */
    void _removeBranchesInsideSoma();

    /**
     * @brief _connectBranches
     */
    void _connectBranches();

    /**
     * @brief _filterLoopsBetweenTwoBranchingPoints
     */
    void _filterLoopsBetweenTwoBranchingPoints();

    /**
     * @brief _filterLoopsAtSingleBranchingPoint
     * Remove any loops at a single branching point. These loops indicate that there is some issue
     * either with the skeleton or the input mesh, for example: two spines are very close and
     * cannot be detected as single objects.
     */
    void _filterLoopsAtSingleBranchingPoint();

    /**
     * @brief _filterSpines
     */
    void _filterSpines();

    /**
     * @brief _updateParent
     * @param branch
     */
    void _updateParent(SkeletonBranch* branch);

    /**
     * @brief _updateParents
     */
    void _updateParents();

    /**
     * @brief _reduceSkeletonToWeightedEdges
     * @return
     */
    SkeletonWeightedEdges _reduceSkeletonToWeightedEdges();

    /**
     * @brief _selectBranchingNodesFromWeightedEdges
     * @param edges
     * @return
     */
    SkeletonNodes _selectBranchingNodesFromWeightedEdges(const SkeletonWeightedEdges& edges);

    /**
     * @brief _getSomaIndexFromGraphNodes
     * @param nodes
     * @return
     */
    int64_t _getSomaIndexFromGraphNodes(const SkeletonNodes& nodes) const;

    /**
     * @brief _constructGraphNodesFromSkeletonNodes
     * @param skeletonNodes
     * @return
     */
    GraphNodes _constructGraphNodesFromSkeletonNodes(const SkeletonNodes& skeletonNodes);

    /**
     * @brief _findShortestPathsFromTerminalNodesToSoma
     * Finds the shortest path from a given node to the soma. This function runs only after the
     * soma is constructed and all the terminal branches are identified.
     * @param edges
     * @param skeletonBranchingNodes
     * @param graphNodes
     * @param somaNodeIndex
     * @return
     */
    EdgesIndices _findShortestPathsFromTerminalNodesToSoma(SkeletonWeightedEdges& edges,
                                                           SkeletonNodes& skeletonBranchingNodes,
                                                           GraphNodes &graphNodes,
                                                           const int64_t& somaNodeIndex);

    /**
     * @brief _constructGraphBranchesFromGraphNodes
     * @param graphNodes
     * @param somaNodeIndex
     * @return
     */
    GraphBranches _constructGraphBranchesFromGraphNodes(
            GraphNodes &graphNodes, const int64_t& somaNodeIndex);

    /**
     * @brief _constructGraphHierarchy
     * @param graphBranches
     */
    void _constructGraphHierarchy(GraphBranches& graphBranches);

    /**
     * @brief _constructSkeletonHierarchy
     * Build the hierarchy of the morphology skeleton, i.e. define the children branches from
     * the parent ones.
     * @param graphBranches
     */
    void _constructSkeletonHierarchy(GraphBranches& graphBranches);

    /**
     * @brief _mergeBranchesWithSingleChild
     * This function merges two connected branches where the parent has only one child.
     * The two branches are marked invalid and a new branch is reconstructed and added to the
     * branches list.
     */
    void _mergeBranchesWithSingleChild();

    /**
     * @brief _detectInactiveBranches
     * Find all the inactive branches in the skeleton that do not contribute to the actual
     * skeleton of the neuron morphology.
     * @param graphEdges
     * @param visitedEdgesIndices
     */
    void _detectInactiveBranches(SkeletonWeightedEdges& graphEdges,
                                 EdgesIndices& visitedEdgesIndices);

    /**
     * @brief _adjustSomaRadius
     */
    void _adjustSomaRadius();

private:

    /**
     * @brief _somaMesh
     * A pointer to the mesh of the segmented soma.
     */
    Mesh* _somaMesh = nullptr;

    /**
     * @brief _somaNode
     * A pointer to the node of the soma.
     */
    SkeletonNode* _somaNode = nullptr;

    /**
     * @brief _roots
     * A list of all the roots in the skeleton, i.e. the branches that emanate from the soma.
     */
    SkeletonBranches _roots;
};

}
