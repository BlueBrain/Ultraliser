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
    NeuronSkeletonizer(Volume *volume, const Mesh* mesh);
    ~NeuronSkeletonizer();

    void constructGraph() override;

    /**
     * @brief getSomaMesh
     * Returns a pointer to the soma mesh.
     * @return
     * Returns a pointer to the soma mesh.
     */
    Mesh* getSomaMesh() const { return _somaMesh; }

    void skeletonizeVolume() override;

    /**
     * @brief exportSomaMesh
     * Exports the soma mesh.
     * @param filePrefix
     */
    void exportSomaMesh(const std::string& filePrefix);


    void skeletonizeVolumeBlockByBlock(const size_t& blockSize = 512,
                                       const size_t& numberOverlappingVoxels = 25,
                                       const size_t& numberZeroVoxels = 5) override;


     void exportIndividualBranches(const std::string& prefix) const;
     void exportSWCFile(const std::string& prefix);


     SkeletonNodes constructSWCTable();
     void collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes, int64_t &swcIndex, int64_t branchingNodeSWCIndex);

     void segmentComponents() override;



private:

    SkeletonNode *_addSomaNode();

    void _segmentSomaMesh(SkeletonNode *somaNode);

    void _segmentSomaVolume();

    void _removeBranchesInsideSoma(SkeletonNode *somaNode);

    void _connectBranches();

    void _processBranchesToYieldCyclicGraph();

    void _filterLoopsBetweenTwoBranchingPoints();

    void _filterLoopsAtSingleBranchingPoint();





    SkeletonWeightedEdges _reduceSkeletonToWeightedEdges();

    SkeletonNodes _selectBranchingNodesFromWeightedEdges(const SkeletonWeightedEdges& edges);

    int64_t _getSomaIndexFromGraphNodes(const SkeletonNodes& nodes) const;

    GraphNodes _constructGraphNodesFromSkeletonNodes(const SkeletonNodes& skeletonNodes);

    EdgesIndices _findShortestPathsFromTerminalNodesToSoma(SkeletonWeightedEdges& edges,
                                                           SkeletonNodes& skeletonBranchingNodes,
                                                           GraphNodes &graphNodes,
                                                           const int64_t& somaNodeIndex);

    GraphBranches _constructGraphBranchesFromGraphNodes(
            GraphNodes &graphNodes, const int64_t& somaNodeIndex);

    void _constructGraphHierarchy(GraphBranches& graphBranches);

    void _constructSkeletonHierarchy(GraphBranches& graphBranches);

    void _mergeBranchesWithSingleChild();

    void _detectInactiveBranches(SkeletonWeightedEdges& graphEdges,
                                 EdgesIndices& visitedEdgesIndices);




private:

    Mesh* _somaMesh = nullptr;

    SkeletonNode* _somaNode = nullptr;


    SkeletonBranches _roots;


};
}
