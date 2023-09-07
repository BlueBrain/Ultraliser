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
#include <algorithms/skeletonization/graphs/Graphs.h>

namespace Ultraliser
{

NeuronSkeletonizer::NeuronSkeletonizer(Volume* volume, const Mesh *mesh)
    : Skeletonizer(volume, mesh)
{
    /// EMPTY CONSTRUCTOR
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
    somaNode->insideSoma = true; // The somatic node is considered inside the soma in the processing
    _nodes.push_back(somaNode);

    _somaNode = somaNode;
    return somaNode;
}

void NeuronSkeletonizer::_segmentSomaMesh(SkeletonNode* somaNode)
{
    // The _somaMesh should be a complex geometry containing overlapping spheres that would define
    // its structure
    _somaMesh = new Mesh();

    // For every node in the skeleton, if the radius is greater than 2.0, then this is a candidate
    // for the soma sample
    Vector3f estimatedSomaCenter(0.f);
    size_t numberSamples = 0;

    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node0 = _nodes[i];
        if (node0->radius >= 2.0)
        {
            Mesh* sample = new IcoSphere(3);
            sample->scale(node0->radius, node0->radius, node0->radius);
            sample->translate(node0->point);
            sample->map(_shellPoints, false);
            _somaMesh->append(sample);
            sample->~Mesh();

            estimatedSomaCenter += node0->point;
            numberSamples++;
        }

    }

    // Normalize
    estimatedSomaCenter /= numberSamples;

    // Update the location of the soma point
    somaNode->point = estimatedSomaCenter;
    somaNode->radius = 0.5; /// TODO: Fix me
}

void NeuronSkeletonizer::_removeBranchesInsideSoma(SkeletonNode* somaNode)
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
                branch->unsetRoot();;
                branch->setValid();
            }

            // If all the branch nodes are located inside the soma, then it is not valid
            else if (countSamplesInsideSoma == branch->nodes.size())
            {
                branch->unsetRoot();
                branch->setInvalid(); // branch->valid = false;
            }

            // Otherwise, it is a branch that is connected to the soma, partially in the soma
            else
            {
                // If the first node is inside the soma, then annotate the branch
                if (firstNode->insideSoma)
                {
                    SkeletonNodes newNodes;
                    newNodes.push_back(somaNode);
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
                    newNodes.push_back(somaNode);

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
                    LOG_WARNING("Undefined case for the branch identification!");
                }
            }
        }
    }
}

void NeuronSkeletonizer::exportSomaMesh(const std::string& filePrefix)
{

}

void NeuronSkeletonizer::_segmentSomaVolume()
{
    // _somaMesh->exportMesh("/data/neuron-meshes/output/soma", true);

    _volume->clear();
    _volume->surfaceVoxelization(_somaMesh, false, false);
    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);
    // _volume->project("/data/neuron-meshes/output/soma", true);

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
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];
        size_t key = _volume->mapTo1DIndexWithoutBoundCheck(
                    node->voxel.x(), node->voxel.y(), node->voxel.z());
        if (std::find(somaVoxels.begin(), somaVoxels.end(), key) != somaVoxels.end())
        {
            // The soma is inside the soma
            node->insideSoma = true;
        }
    }

    // The volume is safe to be deallocated
    _volume->~Volume();
    _volume = nullptr;
}

void NeuronSkeletonizer::collectSWCNodes(const SkeletonBranch* branch, SkeletonNodes& swcNodes,
                                         int64_t& swcIndex, int64_t branchingNodeSWCIndex)
{
    // Get a reference to the nodes of the current branch
    auto& currentBranchNodes = branch->nodes;

    for (size_t i = 1; i < currentBranchNodes.size(); ++i)
    {
        currentBranchNodes[i]->swcIndex = swcIndex;

        if (i == 1)
        {
            currentBranchNodes[i]->prevSampleSWCIndex = branchingNodeSWCIndex;
        }
        else
        {
            currentBranchNodes[i]->prevSampleSWCIndex= swcIndex - 1;
        }

        swcIndex++;
        swcNodes.push_back(currentBranchNodes[i]);
    }

    const int64_t branchingIndex = swcIndex - 1;

    for (size_t i = 0; i < branch->children.size(); ++i)
    {
        if (branch->isValid())
        {
            collectSWCNodes(branch->children[i], swcNodes, swcIndex, branchingIndex);
        }
    }
}


SkeletonNodes NeuronSkeletonizer::constructSWCTable()
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

    // Get all the root branches
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];
        if (branch->isRoot() && branch->isValid())
        {
            // The branching index is that of the soma
            collectSWCNodes(branch, swcNodes, swcIndex, 1);
        }
    }

    return swcNodes;
}



void NeuronSkeletonizer::exportSWCFile(const std::string& prefix)
{
    // Start the timer
    TIMER_SET;

    // Construct the file path
    std::string filePath = prefix + SWC_EXTENSION;
    LOG_STATUS("Exporting Neuron to SWC file: [ %s ]", filePath.c_str());

    auto swcNodes = constructSWCTable();

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    auto somaNode = swcNodes[0];
    stream << somaNode->swcIndex << " "
           << "1" << " "
           << somaNode->point.x() << " "
           << somaNode->point.y() << " "
           << somaNode->point.z() << " "
           << somaNode->radius << " "
           << "-1" << "\n";

    for (size_t i = 1; i < swcNodes.size(); ++i)
    {
        auto swcNode = swcNodes[i];
        stream << swcNode->swcIndex << " "
               << "3" << " "
               << swcNode->point.x() << " "
               << swcNode->point.y() << " "
               << swcNode->point.z() << " "
               << swcNode->radius << " "
               << swcNode->prevSampleSWCIndex << "\n";
    }

    // Close the file
    stream.close();
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
        if (!_branches[i]->isValid() || _branches[i]->isDuplicate() || _branches[i]->nodes.size() == 0) continue;
        // if (!_branches[i]->root) continue;

        LOOP_PROGRESS(progress, _branches.size());
        ++progress;

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
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    stream.close();



    filePath = prefix + "-non-somatic-samples" + TXT_EXTENSION;
    stream.open(filePath, std::ios::out);
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];

        if (!node->insideSoma)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
    }
     stream.close();


    filePath = prefix + "-samples" + TXT_EXTENSION;
    stream.open(filePath, std::ios::out);
    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];

        if (node->insideSoma)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
    }
     stream.close();
}

void NeuronSkeletonizer::_filterLoopsBetweenTwoBranchingPoints()
{
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
    }
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
    // A list that should have all the branches of the valid branches of the neuron represented
    // by weighted edges
    SkeletonWeightedEdges edges;
    for (auto& branch: _branches)
    {
        // The branch must be valid to be able to have a valid graph
        if (branch->isValid())
        {
            // Reset the nodes of the branch
            branch->nodes.front()->visited = false;
            branch->nodes.back()->visited = false;

            // Create a weighted edge and append it to the list
            SkeletonWeightedEdge* edge = new SkeletonWeightedEdge(branch);
            edges.push_back(edge);
        }
    }

    // Return the resulting edges array that will be used for constructing the graph
    return edges;
}


SkeletonNodes NeuronSkeletonizer::_selectBranchingNodesFromWeightedEdges(const SkeletonWeightedEdges& edges)
{
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
    }

    return nodes;
}

int64_t NeuronSkeletonizer::_getSomaIndexFromGraphNodes(const SkeletonNodes& nodes) const
{
    int64_t somaNodeIndex = -1;
    for (size_t i = 0; nodes.size(); ++i)
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

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < skeletonNodes.size(); ++i)
    {
        const auto& skeletonNode = skeletonNodes[i];
        graphNodes[i] = (new GraphNode(skeletonNode->graphIndex, skeletonNode->index));
    }

    // Return the graph nodes list
    return graphNodes;
}

EdgesIndices NeuronSkeletonizer::_findShortestPathsFromTerminalNodesToSoma(
        SkeletonWeightedEdges& edges,
        SkeletonNodes &skeletonBranchingNodes, GraphNodes &graphNodes,
        const int64_t& somaNodeIndex)
{
    // The indices of all the edges that have been traversed
    EdgesIndices edgesIndices;

    // Search for all the terminal nodes
    for (size_t i = 0; i < skeletonBranchingNodes.size(); i++)
    {
        // The node must be a terminal
        if (skeletonBranchingNodes[i]->terminal)
        {
            // Construct the ShortestPathFinder solution
            std::unique_ptr< ShortestPathFinder > pathFinder =
                    std::make_unique< ShortestPathFinder >(edges, skeletonBranchingNodes.size());

            // Find the path between the terminal node and the soma node
            auto terminalToSomaPath = pathFinder->findPath(
                        skeletonBranchingNodes[i]->graphIndex, somaNodeIndex);

            // Reverse the terminal to soma path to have the correct order
            std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());

            for (size_t j = 0; j < terminalToSomaPath.size() - 1; ++j)
            {
                auto currentNodeIndex = terminalToSomaPath[j];
                auto nextNodeIndex = terminalToSomaPath[j + 1];

                // Add the edge indices to the list
                edgesIndices.push_back(EdgeIndex(currentNodeIndex, nextNodeIndex));

                // If the next node index is not in the current node index, then add it
                if (!graphNodes[currentNodeIndex]->isNodeInChildren(nextNodeIndex))
                {
                    graphNodes[currentNodeIndex]->children.push_back(graphNodes[nextNodeIndex]);
                }
            }
        }
    }

    // Return the EdgesIndices list
    return edgesIndices;
}

GraphBranches NeuronSkeletonizer::_constructGraphBranchesFromGraphNodes(
        GraphNodes &graphNodes, const int64_t& somaNodeIndex)
{
    // Use a new index to label graph branches
    size_t branchGraphIndex = 0;

    // A list of all the constructed GraphBranches
    GraphBranches graphBranches;

    // Construct the valid branches at the end
    for (size_t i = 0; i < graphNodes.size(); i++)
    {
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

                        std::cout << "Branch: " << branch->index << " ";
                        std::cout << "(" << firstNodeIndex << ", " << lastNodeIndex << ") ";

                        //continue;
                    }

                    // If the branch is valid
                    //if (branch->valid)
                }
            }

            std::cout << "\n";
        }
    }

    return graphBranches;
}

void NeuronSkeletonizer::_constructGraphHierarchy(GraphBranches& graphBranches)
{
    for (size_t i = 0; i < graphBranches.size(); ++i)
    {
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
}

void NeuronSkeletonizer::_constructSkeletonHierarchy(GraphBranches& graphBranches)
{
    for(size_t i = 0; i < graphBranches.size(); ++i)
    {
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
            skeletonBranch->children.push_back(_branches[graphBranch->children[j]->skeletonIndex]);
        }
    }
}

void NeuronSkeletonizer::_mergeBranchesWithSingleChild()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // A reference to the branch
        auto& branch = _branches[i];

        // If the branch is valid and has a single child, merge it with the parent
        if (branch->isValid() && branch->children.size() == 1)
        {
            // Append the SkeletonNodes from the child branch
            for(size_t j = 1; j < branch->children[0]->nodes.size(); ++j)
            {
                branch->nodes.push_back(branch->children[0]->nodes[j]);
            }

            // Invalidate the child branch
            branch->children[0]->setInvalid();
            branch->children[0]->setDuplicate();

            // If the child branch has any children, then update the children of the parent
            if (branch->children[0]->children.size() > 0)
            {
                branch->children = branch->children[0]->children;
            }
        }
    }
}


void NeuronSkeletonizer::_processBranchesToYieldCyclicGraph()
{
    // Initially, and before constructing the graph, remove the loops between two branching points
    _filterLoopsBetweenTwoBranchingPoints();

    // Remove the loops at a single branching point, i.e. starting and ending at the same node.
    _filterLoopsAtSingleBranchingPoint();

    // Reduce the skeleton into a list of SkeletonWeightedEdge's
    SkeletonWeightedEdges weighteEdges = _reduceSkeletonToWeightedEdges();

    // Get a list of all the branching nodes within the skeleton from the SkeletonWeightedEdges list
    SkeletonNodes skeletonBranchingNodes = _selectBranchingNodesFromWeightedEdges(weighteEdges);

    // Get the soma node index within the weighted graph to be used to construct the paths
    int64_t somaNodeIndex = _getSomaIndexFromGraphNodes(skeletonBranchingNodes);

    // Construct the graph nodes list
    GraphNodes graphNodes = _constructGraphNodesFromSkeletonNodes(skeletonBranchingNodes);

    // Find the shortest paths of all the terminals and get a list of the indices of the active edges
    EdgesIndices edgeIndices = _findShortestPathsFromTerminalNodesToSoma(weighteEdges,
                                                                         skeletonBranchingNodes,
                                                                         graphNodes, somaNodeIndex);



    for (size_t i = 0; i < edgeIndices.size(); i++)
    {
           std::cout << "edges: " << edgeIndices[i].first << ", " << edgeIndices[i].second << "\n";
    }

    for (size_t i = 0; i < weighteEdges.size(); i++)
    {
           std::cout << "w Edge: " << weighteEdges[i]->node1->graphIndex << ", " << weighteEdges[i]->node2->graphIndex << "\n";
    }



    // From the edge indices, find the branches that were traversed in the graph and then the skeleton
    for (size_t i = 0; i < weighteEdges.size(); i++)
    {
        auto& edge = weighteEdges[i];
        if (!edge->visited)
        {
            for (size_t j = 0; j < edgeIndices.size(); ++j)
            {
                auto& indices = edgeIndices[j];
                if (edge->hasTerminalGraphNodes(indices.first, indices.second))
                {
                    edge->visited = true;
                    break;
                }
            }
        }
    }


    for (size_t i = 0; i < weighteEdges.size(); ++i)
    {
        if (weighteEdges[i]->visited)
        {
            std::cout << "Visited Edge : " << weighteEdges[i]->branch->index << "\n";
        }
        else
        {
            weighteEdges[i]->branch->setInvalid();
            std::cout << "Useless Edge : " << weighteEdges[i]->branch->index << "\n";
        }
    }




    // Construct the GraphBranches from the GraphNodes
    GraphBranches graphBranches = _constructGraphBranchesFromGraphNodes(graphNodes, somaNodeIndex);

    // Construct the hierarchy of the graph
    _constructGraphHierarchy(graphBranches);

    // Construct the hierarchy of the skeleton
    _constructSkeletonHierarchy(graphBranches);

    // Merge the branches that have a single child
    _mergeBranchesWithSingleChild();





//    std::cout << "Graphs " << graphBranches.size() << "\n";
//    for(size_t i = 0; i < graphBranches.size(); ++i)
//    {
//        auto& gBranch = graphBranches[i];

//        std::cout << "* Branch: " << gBranch->skeletonIndex << "->";
//        for(size_t j = 0; j < edgeIndices.size(); ++j)
//        {
//            auto edgePair = edgeIndices[j];

//            if (gBranch->hasTerminalNodes(edgePair.first, edgePair.second))
//            {
//                std::cout << "<" << edgePair.first << ", " << edgePair.second << ">, ";
//                std::cout << "[" << gBranch->firstNodeIndex << ", " << gBranch->lastNodeIndex << "],,  ";

//                gBranch->active = true;
//                std::cout << "Active: " << gBranch->index << "\n";
//            }


//        }
//        std::cout << "\n";
//    }

//    for(size_t i = 0; i < graphBranches.size(); ++i)
//    {
//        auto& gBranch = graphBranches[i];
//        if (gBranch->active)
//        {
//            std::cout << i << " Active: " << gBranch->skeletonIndex << "\n";
//        }
//        else
//        {
//             std::cout << i << " Inactive: " << gBranch->skeletonIndex << "\n";
//        }
//    }



    // merge branches with a single child



    for(size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isRoot())
            _branches[i]->printTree();
    }
}


void identifyTerminalConnections(SkeletonBranches& branches)
{
    for (size_t i = 0; i < branches.size(); ++i)
    {
        auto& iBranch = branches[i];

        // Invalid branch, ignore
        if (!iBranch->isValid()) continue;

        auto& iBranchT1 = iBranch->nodes.front();
        auto& iBranchT2 = iBranch->nodes.back();

        for (size_t j = 0; j < branches.size(); ++j)
        {
            auto& jBranch = branches[j];

            // Invalid branch, ignore
            if (!jBranch->isValid()) continue;

            // Same branch, next branch
            if (iBranch->index == jBranch->index) continue;

            auto& jBranchT1 = jBranch->nodes.front();
            auto& jBranchT2 = jBranch->nodes.back();

            // A new connection to the T1 terminal
            if (iBranchT1->index == jBranchT1->index)
            {
                // The sample must not be the soma
                if (jBranchT1->isSoma) continue;

                iBranch->t1Connections.push_back(jBranch);
            }

            if (iBranchT1->index == jBranchT2->index)
            {
                // The sample must not be the soma
                if (jBranchT2->isSoma) continue;

                iBranch->t1Connections.push_back(jBranch);
            }

            if (iBranchT2->index == jBranchT1->index)
            {
                // The sample must not be the soma
                if (jBranchT1->isSoma) continue;

                iBranch->t2Connections.push_back(jBranch);
            }

            // A new connection to the T2 terminal
            if (iBranchT2->index == jBranchT2->index)
            {
                // The sample must not be the soma
                if (jBranchT2->isSoma) continue;

                iBranch->t2Connections.push_back(jBranch);
            }
        }
    }
}

void confirmTerminalsBranches(SkeletonBranches& branches)
{
    for (size_t i = 0; i < branches.size(); ++i)
    {
        // Reference to the branch
        auto& branch = branches[i];

        // Invalid branches, next
        if (!branch->isValid()) continue;

        // If root branch, then next
        if (branch->isRoot()) continue;

        // If terminal branch, then adjust the samples
        if (branch->t1Connections.size() == 0 || branch->t2Connections.size() == 0)
        {
            // Update the status
            branch->setTerminal();

            // Reverse the sample if needed
            if (branch->nodes.back()->edgeNodes.size() > 1)
            {
                std::reverse(std::begin(branch->nodes), std::end(branch->nodes));
            }
        }
    }
}

SkeletonBranches _detectChildren(SkeletonBranch* currentBranch, SkeletonBranches& allBranches)
{
    SkeletonBranches childrenBranches;

    // Get a reference to the last node in the root
    auto rootLastNode = currentBranch->nodes.back();

    for (size_t i = 0; i < allBranches.size(); ++i)
    {
        // Reference to the branch
        auto& branch = allBranches[i];

        // Invalid branches, next
        if (!branch->isValid()) continue;

        // If root branch, then next
        if (branch->isRoot()) continue;

        // If the currentBranch and the branch have the same index, then invalid
        if (currentBranch->index == branch->index) continue;

        // Access the nodes of the branch
        auto& branchFirstNode = branch->nodes.front();
        auto& branchLastNode = branch->nodes.back();

        // If the last node of the root is the first node of the branch, then it is a child
        if (rootLastNode->index == branchFirstNode->index)
        {
            childrenBranches.push_back(branch);

            // Add the branch to the children list
            // currentBranch->children.push_back(branch);

            // Add the root to be a parent of the branch
            // branch->parents.push_back(currentBranch);
        }

        // If the last node of the root is the last node of the branch, then it is a
        // reversed child
        if (rootLastNode->index == branchLastNode->index)
        {
            // Reverse the order of the nodes
            std::reverse(std::begin(branch->nodes), std::end(branch->nodes));

            // Add the branch to the children list
            currentBranch->children.push_back(branch);

            // Add the root to be a parent of the branch
            branch->parents.push_back(currentBranch);
        }
    }

    return childrenBranches;
}


void NeuronSkeletonizer::_connectBranches()
{
    // Identify the connections at the terminals of each branch
    identifyTerminalConnections(_branches);

    // Roots, terminals and others
    confirmTerminalsBranches(_branches);
}

void NeuronSkeletonizer::constructGraph()
{
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxels();

    // Assign accurate radii to the nodes of the graph, i.e. inflate the nodes
    _inflateNodes();

    // Connect the nodes to construct the edges of the graph
    _connectNodes(indicesMapper);

    // Add a virtual soma node, until the soma is reconstructed later
    auto somaNode = _addSomaNode();

    // Re-index the samples, for simplicity
    OMP_PARALLEL_FOR for (size_t i = 1; i <= _nodes.size(); ++i) { _nodes[i - 1]->index = i; }

    // Segmentthe soma mesh from the branches
    _segmentSomaMesh(somaNode);

    // Segment soma volume
     _segmentSomaVolume();

     // Remove the triangular configurations, based on the edges
     _removeTriangleLoops();

     // Reconstruct the sections, or the branches from the nodes
     _buildBranchesFromNodes(_nodes);

    // Validate the branches, and remove the branches inside the soma
    _removeBranchesInsideSoma(somaNode);

    // The roots have been identified
    _connectBranches();
}

void NeuronSkeletonizer::segmentComponents()
{
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

    _processBranchesToYieldCyclicGraph();
}
}
