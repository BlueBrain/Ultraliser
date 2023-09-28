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

    // Initially, we set the soma node to some value that does not make any conflict
    somaNode->radius = 0.1;
    _nodes.push_back(somaNode);
    _somaNode = somaNode;
    return somaNode;
}

void NeuronSkeletonizer::_segmentSomaMesh(SkeletonNode* somaNode)
{
    LOG_STATUS("Segmenting Soma");

    // The _somaMesh should be a complex geometry containing overlapping spheres that would define
    // its structure
    _somaMesh = new Mesh();

    // For every node in the skeleton, if the radius is greater than 2.0, then this is a candidate
    // for the soma sample
    Vector3f estimatedSomaCenter(0.f);
    size_t numberSamples = 0;

    // Collecting a subset of the samples, not all of them are needed
    std::map <size_t, float> interSomaticNodes;

    for (size_t i = 0; i < _nodes.size(); ++i)
    {
        auto& node = _nodes[i];

        /// TODO: This is a magic value, it works now, but we need to find an optimum value based
        // on some statistical analysis.

        if (node->radius >= 2.0)
        {
            interSomaticNodes.insert({node->index, node->radius});
        }
    }

    // Sort the nodes by radius
    std::vector< std::pair< size_t , float > > pairsVector = sortIndexRadiusMap(interSomaticNodes);

    // Reverse
    std::reverse(pairsVector.begin(), pairsVector.end());

    size_t numberSelectedNodes = 50;

    TIMER_SET;
    LOOP_STARTS("Detecting Soma Nodes");
    for (size_t i = 0; i < numberSelectedNodes; ++i)
    {
        auto& node0 = _nodes[pairsVector[i].first - 1];

        // std::cout << pairsVector[i].first << ", " << node0->index << ",,, " << pairsVector[i].second << ", " << node0->radius << "\n";


        Mesh* sample = new IcoSphere(3);
        sample->scale(node0->radius, node0->radius, node0->radius);
        sample->translate(node0->point);
        sample->map(_shellPoints, false);

        _somaMesh->append(sample);
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

void NeuronSkeletonizer::exportSomaMesh(const std::string& filePrefix,
                                        const bool& formatOBJ = false,
                                        const bool& formatPLY = false,
                                        const bool& formatOFF = false,
                                        const bool& formatSTL = false)
{
    const std::string somaMeshPrefix = filePrefix + "-soma";
    _somaMesh->exportMesh(somaMeshPrefix, formatOBJ, formatPLY, formatOFF, formatSTL);
}

void NeuronSkeletonizer::_segmentSomaVolume()
{
    TIMER_SET;
    LOG_STATUS("Segmenting Soma Volume");

    _volume->clear();
    _volume->surfaceVoxelization(_somaMesh, false, false);
    _volume->solidVoxelization(Volume::SOLID_VOXELIZATION_AXIS::XYZ);

    // Get a reference to the occupancy ranges in case it is not computed
    auto occupancyRanges = _volume->getOccupancyRanges();

    std::vector< std::vector< size_t > > perSliceSomaVoxels;
    perSliceSomaVoxels.resize(_volume->getWidth());

    OMP_PARALLEL_FOR
    for (size_t i = 0; i < occupancyRanges.size(); ++i)
    {
        for (size_t j = 0; j < occupancyRanges[i].size(); ++j)
        {
            for (size_t k = 0; k < occupancyRanges[i][j].size(); ++k)
            {
                for (size_t w = occupancyRanges[i][j][k].lower; w <= occupancyRanges[i][j][k].upper; w++)
                {
                    if (_volume->isFilled(i, j, w))
                    {
                        perSliceSomaVoxels[i].push_back(_volume->mapTo1DIndexWithoutBoundCheck(i, j, w));
                    }
                }
            }
        }
    }

    std::vector< size_t > somaVoxels;
    for (size_t i = 0; i < perSliceSomaVoxels.size(); ++i)
    {
        somaVoxels.insert(somaVoxels.end(), perSliceSomaVoxels[i].begin(), perSliceSomaVoxels[i].end());
        perSliceSomaVoxels[i].clear();
        perSliceSomaVoxels[i].shrink_to_fit();
    }

    perSliceSomaVoxels.clear();
    perSliceSomaVoxels.shrink_to_fit();

 // TODO: This could be parallelized
//    for (size_t i = 0; i < _volume->getNumberVoxels(); ++i)
//    {
//        if (_volume->isFilled(i))
//        {
//            somaVoxels.push_back(i);
//        }
//    }

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

    LOG_STATS(GET_TIME_SECONDS);
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
    TIMER_SET;
    LOOP_STARTS("Constructing SWC Table");
    const size_t numberBranches = _branches.size();
    for (size_t i = 0; i < numberBranches ; ++i)
    {
        LOOP_PROGRESS(i, numberBranches);

        auto& branch = _branches[i];
        if (branch->isRoot() && branch->isValid())
        {
            // The branching index is that of the soma
            collectSWCNodes(branch, swcNodes, swcIndex, 1);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

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

    LOOP_STARTS("Writing SWC Table");
    const size_t numberSWCNodes = swcNodes.size();
    for (size_t i = 1; i < numberSWCNodes; ++i)
    {
        LOOP_PROGRESS(i, numberSWCNodes);

        auto swcNode = swcNodes[i];
        stream << swcNode->swcIndex << " "
               << "3" << " "
               << swcNode->point.x() << " "
               << swcNode->point.y() << " "
               << swcNode->point.z() << " "
               << swcNode->radius << " "
               << swcNode->prevSampleSWCIndex << "\n";
    }
    LOOP_DONE;

    // Close the file
    stream.close();

    LOG_STATS(GET_TIME_SECONDS);
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
        if (!_branches[i]->isValid() || _branches[i]->nodes.size() == 0) continue;
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
    TIMER_SET;
    LOG_STATUS("Filtering Loops Between Two Branches");
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        LOOP_PROGRESS(i, _branches.size());

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
    LOG_STATUS("Created Weighted Graph from Skeleton");

    // A list that should have all the branches of the valid branches of the neuron represented
    // by weighted edges
    SkeletonWeightedEdges edges;
    const auto branchesCount = _branches.size();
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        auto& branch = _branches[i];

        LOOP_PROGRESS(i, branchesCount);

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
        LOOP_PROGRESS(i, edges.size());

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
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

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

    // OMP_PARALLEL_FOR
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

    // Search for all the terminal nodes
    PROGRESS_SET;
    LOOP_STARTS("Detecting Paths");
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < terminalNodes.size(); i++)
    {
        EdgesIndices& perTerminalEdgesIndices = edgesIndicesList[i];

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, terminalNodes.size());
        PROGRESS_UPDATE;

        // Construct the ShortestPathFinder solution
        std::unique_ptr< ShortestPathFinder > pathFinder =
                std::make_unique< ShortestPathFinder >(edges, skeletonBranchingNodes.size());

        // Find the path between the terminal node and the soma node
        auto terminalToSomaPath = pathFinder->findPath(terminalNodes[i]->graphIndex, somaNodeIndex);

        // Reverse the terminal to soma path to have the correct order
        std::reverse(terminalToSomaPath.begin(), terminalToSomaPath.end());

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

    }
    LOOP_DONE;

    // Initially, log the time for the loop
    LOG_STATS(GET_TIME_SECONDS);


    LOOP_STARTS("Composing Path Edges");
    // Clear the terminal nodes list
    terminalNodes.clear();
    terminalNodes.shrink_to_fit();

    // The indices of all the edges that have been traversed
    EdgesIndices edgesIndices;
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

    // Clean the list used to collect the edges in parallel
    edgesIndicesList.clear();
    edgesIndicesList.shrink_to_fit();

    // Log the time for the whole process after the edge construction
    LOG_STATS(GET_TIME_SECONDS);

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
            skeletonBranch->children.push_back(_branches[graphBranch->children[j]->skeletonIndex]);
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void NeuronSkeletonizer::_mergeBranchesWithSingleChild()
{
    TIMER_SET;
    LOG_STATUS("Merging Branches with Single Child");
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        LOOP_PROGRESS(i, _branches.size());

        // A reference to the branch
        auto& branch = _branches[i];

        // If the branch is valid and has a single child, merge it with the parent
        if (branch->isValid() && branch->active && branch->children.size() == 1)
        {
            // Append the SkeletonNodes from the child branch
            for(size_t j = 1; j < branch->children[0]->nodes.size(); ++j)
            {
                branch->nodes.push_back(branch->children[0]->nodes[j]);
            }

            // Invalidate the child branch
            branch->children[0]->setInvalid();

            // If the child branch has any children, then update the children of the parent
            if (branch->children[0]->children.size() > 0)
            {
                branch->children = branch->children[0]->children;
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
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
            somaRadius += branch->nodes[1]->point.distance(_somaNode->point);;
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

                _updateParent(child);
            }
        }
    }
}

void NeuronSkeletonizer::_filterSpines()
{
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        if (_branches[i]->isTerminal())
        {
            auto length = _branches[i]->computeLength();
            if (length < 6.0)
            {
                _branches[i]->setInvalid();
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

    // Get the soma node index within the weighted graph
    int64_t somaNodeIndex = _getSomaIndexFromGraphNodes(skeletonBranchingNodes);

    // Construct the graph nodes list
    GraphNodes graphNodes = _constructGraphNodesFromSkeletonNodes(skeletonBranchingNodes);

    /// After having the weighted edges and the nodes computed, compute the number of components in
    /// the graph, if the result is more than 1 then then re-connect them to be in a single graph

    auto graph = new Graph(weighteEdges, graphNodes);

    auto components = graph->getComponents();
    std::cout << "Number Components " << components.size() << "\n";

    // Find the shortest paths of all the terminals and get a list of the indices of the active edges
    EdgesIndices edgeIndices = _findShortestPathsFromTerminalNodesToSoma(
                weighteEdges, skeletonBranchingNodes, graphNodes, somaNodeIndex);

    // Construct the GraphBranches from the GraphNodes
    GraphBranches graphBranches = _constructGraphBranchesFromGraphNodes(graphNodes, somaNodeIndex);

    // Construct the hierarchy of the graph
    _constructGraphHierarchy(graphBranches);

    // Construct the hierarchy of the skeleton
    _constructSkeletonHierarchy(graphBranches);

    // merge branches with a single child
    _mergeBranchesWithSingleChild();

    // Invalidate the inactive branches
    _detectInactiveBranches(weighteEdges, edgeIndices);

    // Adkjust the soma radius
    _adjustSomaRadius();

    // Update all the parents
    _updateParents();

    // Filter the synapses
    _filterSpines();
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
    std::map< size_t, size_t > indicesMapper = _extractNodesFromVoxelsParallel();

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
    // _buildBranchesFromNodes(_nodes);

    // std::cout << "Branches: " << _branches.size() << "\n";
    // std::cout << "Nodes (Samples): " << _nodes.size() << "\n";

    _processBranchesToYieldCyclicGraph();
}
}
