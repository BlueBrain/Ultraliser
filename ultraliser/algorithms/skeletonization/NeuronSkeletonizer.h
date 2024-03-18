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
#include <algorithms/skeletonization/graphs/Graph.h>

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
    NeuronSkeletonizer(Volume *volume,
                       const bool &removeSpines = true,
                       const bool &useAcceleration = true,
                       const bool &debugSkeleton = false,
                       const std::string debuggingPrefix = NONE);
    ~NeuronSkeletonizer();

    /**
     * @brief constructGraph
     * @param verbose
     */
    void constructGraph(const bool verbose = VERBOSE) override;

    /**
     * @brief segmentComponents
     * @param verbose
     */
    void segmentComponents(const bool verbose = VERBOSE) override;

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
      * @param resampleSkeleton
      * If this flag is set, the morphology skeleton will be adaptively resampled to remove
      * useless samples and create an optimum skeleton. False by default.
      * @param verbose
      */
    void exportSWCFile(const std::string& prefix, const bool &resampleSkeleton = false,
                       const bool verbose = VERBOSE);

    /**
     * @brief constructSWCTable
     * To facilitate exporting the skeleton into an SWC file, we construct a linear list of
     * all the SWC in order and then export this table in the @exportSWCFile function.
     * @return
     * A list of all the skeleton nodes, constructed in order, and ready for being exported into
     * SWC file.
     * @param resampleSkeleton
     * If this flag is set, the morphology skeleton will be adaptively resampled to remove
     * useless samples and create an optimum skeleton. False by default.
     * @param verbose
     */
    SkeletonNodes constructSWCTable(const bool &resampleSkeleton = false,
                                    const bool verbose = VERBOSE);

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
     * @brief exportSomaProxyMesh
     * Exports the proxy mesh of the reconstructed soma. This is the initial mesh object that is
     * used to reconstruct the mesh.
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
    void exportSomaProxyMesh(const std::string& prefix,
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

    void segmentSpines(const bool verbose = VERBOSE);

    /**
     * @brief reconstructSpineProxyMorphologies
     * @return
     */
    SpineMorphologies reconstructSpineProxyMorphologies();

    /**
     * @brief reconstructSpineMeshes
     * @param neuronMesh
     * @param voxelsPerMicron
     * @param edgeGap
     * @return
     */
    Meshes reconstructSpineMeshes(Mesh* neuronMesh,
                                  const float& voxelsPerMicron,
                                  const float& edgeGap = 0.1);

    /**
     * @brief exportBranches
     * @param prefix
     * @param state
     * @param verbose
     */
    void exportBranches(const std::string &prefix, const SkeletonBranch::BRANCH_STATE state,
                        const bool verbose = VERBOSE);

    /**
     * @brief exportSpineLocations
     * @param prefix
     * @param verbose
     */
    void exportSpineLocations(const std::string& prefix, const bool verbose = VERBOSE) const;

    /**
     * @brief exportSpineExtents
     * @param prefix
     * @param verbose
     */
    void exportSpineExtents(const std::string& prefix, const bool verbose = VERBOSE) const;



private:

    /**
     * @brief _addSomaNode
     * Adds a dedicated node for the soma. This node is updated later after the soma is segmented
     * from the mesh and all the arbors are detected.
     */
    void _addSomaNode();

    /**
     * @brief _segmentSomaMesh
     * @param verbose
     */
    void _segmentSomaMesh(const bool verbose = VERBOSE);

    /**
     * @brief _identifySomaticNodes
     * @param verbose
     */
    void _identifySomaticNodes(const bool verbose = VERBOSE);

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
    void _filterLoopsBetweenTwoBranchingPoints(const bool verbose = VERBOSE);

    /**
     * @brief _filterLoopsAtSingleBranchingPoint
     * Remove any loops at a single branching point. These loops indicate that there is some issue
     * either with the skeleton or the input mesh, for example: two spines are very close and
     * cannot be detected as single objects.
     */
    void _filterLoopsAtSingleBranchingPoint();

    /**
     * @brief _updateParent
     * @param branch
     */
    void _updateParent(SkeletonBranch* branch);

    /**
     * @brief _updateParents
     */
    void _updateParents(const bool verbose = VERBOSE);

    /**
     * @brief _reduceSkeletonToWeightedEdges
     * @return
     */
    SkeletonWeightedEdges _reduceSkeletonToWeightedEdges(const bool verbose = VERBOSE);

    /**
     * @brief _selectBranchingNodesFromWeightedEdges
     * @param edges
     * @return
     */
    SkeletonNodes _selectBranchingNodesFromWeightedEdges(const SkeletonWeightedEdges& edges,
                                                         const bool verbose = VERBOSE);

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
                                                           const int64_t& somaNodeIndex,
                                                           const bool verbose = VERBOSE);

    /**
     * @brief _constructGraphBranchesFromGraphNodes
     * @param graphNodes
     * @param somaNodeIndex
     * @return
     */
    GraphBranches _constructGraphBranchesFromGraphNodes(GraphNodes &graphNodes,
                                                        const int64_t& somaNodeIndex,
                                                        const bool verbose = VERBOSE);

    /**
     * @brief _constructGraphHierarchy
     * @param graphBranches
     */
    void _constructGraphHierarchy(GraphBranches& graphBranches, const bool verbose = VERBOSE);

    /**
     * @brief _constructSkeletonHierarchy
     * Build the hierarchy of the morphology skeleton, i.e. define the children branches from
     * the parent ones.
     * @param graphBranches
     */
    void _constructSkeletonHierarchy(GraphBranches& graphBranches, const bool verbose = VERBOSE);

    /**
     * @brief _mergeBranchesWithSingleChild
     * This function merges two connected branches where the parent has only one child.
     * The two branches are marked invalid and a new branch is reconstructed and added to the
     * branches list.
     */
    void _mergeBranchesWithSingleChild();


    void _filterSpineCandidates();

    void _removeSpinesss();

    size_t _estimateNumberSpineCandidates();

    void _removeRootSpines();



    void _detectBasePaths();

    void _filterShortTerminalBranches();


    void _detectSpines();





    void _mergeBranchesAfterFilteringSpines();

    /**
     * @brief _detectInactiveBranches
     * Find all the inactive branches in the skeleton that do not contribute to the actual
     * skeleton of the neuron morphology.
     * @param graphEdges
     * @param visitedEdgesIndices
     */
    void _detectInactiveBranches(SkeletonWeightedEdges& graphEdges,
                                 EdgesIndices& visitedEdgesIndices,
                                 const bool verbose = VERBOSE);

    /**
     * @brief _adjustSomaRadius
     */
    void _adjustSomaRadius(const bool verbose = VERBOSE);

    void _reconstructSomaMeshFromProxy(const bool verbose = VERBOSE);


    /**
     * @brief _verifyGraphConnectivity
     */
    void _verifyGraphConnectivity(SkeletonEdges &edges);

    /**
     * @brief _verifyGraphConnectivityToMainPartition
     * @param components
     */
    void _verifyGraphConnectivityToMainPartition(GraphComponents &components, SkeletonEdges& edges);

    void _verifyGraphConnectivityToClosestPartition(SkeletonEdges& edges);


    void _findClosestNodesInTwoPartitions(GraphComponent& partition1, GraphComponent& partition2,
                                          size_t* partition1NodeIndex, size_t* partition2NodeIndex,
                                          float* distance);
    void _connectPartition(GraphComponents& partitions, const size_t &partitionIndex,
                           SkeletonEdges &edges);

    void _exportSomaticNodes(const std::string prefix);






private:
    /**
     * @brief _removeSpines
     */
    const bool _removeSpines;

    /**
     * @brief _somaProxyMesh
     * A pointer to the proxy mesh of the segmented soma.
     */
    Mesh* _somaProxyMesh = nullptr;

    /**
     * @brief _somaMesh
     * A pointer to the final mesh of the segmented soma, after processing the proxy mesh.
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

    /**
     * @brief _spines
     */
    std::vector< SkeletonBranches > _spines;


    SkeletonBranches _spineRoots;
};

}
