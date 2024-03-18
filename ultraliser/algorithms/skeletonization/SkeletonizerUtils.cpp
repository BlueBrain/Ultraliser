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

#include "SkeletonizerUtils.h"

namespace Ultraliser
{

void removeEdgeNode(SkeletonNode* node, SkeletonNode* edgeNodeMarkedForRemoval)
{
    auto it = std::find(node->edgeNodes.begin(), node->edgeNodes.end(), edgeNodeMarkedForRemoval);

    if(it != node->edgeNodes.end())
        node->edgeNodes.erase(it);
    node->edgeNodes.shrink_to_fit();
}

void collapseTriangleIntoNode(SkeletonNodes& nodes,
                              SkeletonNode* n1, SkeletonNode* n2, SkeletonNode* n3)
{
    SkeletonNode* centerNode = new SkeletonNode();
    centerNode->point = (n1->point + n2->point + n3->point) / 3.f;
    centerNode->radius = (n1->radius + n2->radius + n3->radius) / 3.f;

    centerNode->index = nodes.back()->index + 1;
    centerNode->branching = true;
    centerNode->terminal = false;

    removeEdgeNode(n1, n2); removeEdgeNode(n1, n3);
    removeEdgeNode(n2, n1); removeEdgeNode(n2, n3);
    removeEdgeNode(n3, n1); removeEdgeNode(n3, n2);

    centerNode->edgeNodes.push_back(n1);
    centerNode->edgeNodes.push_back(n2);
    centerNode->edgeNodes.push_back(n3);

    n1->edgeNodes.push_back(centerNode);
    n2->edgeNodes.push_back(centerNode);
    n3->edgeNodes.push_back(centerNode);

    if (n1->insideSoma || n2->insideSoma || n3->insideSoma)
    {
        centerNode->insideSoma = true;
    }

    nodes.push_back(centerNode);
}

bool areConnected(const SkeletonNode* n1, const SkeletonNode* n2)
{
    for (size_t i = 0; i < n1->edgeNodes.size(); ++i)
    {
        if (n1->edgeNodes[i]->index == n2->index)
        {
            return true;
        }
    }

    return false;
}

bool isTriangleNode(const SkeletonNode* n, SkeletonNodes& connectedEdgeNodes)
{
    for (size_t i = 0; i < n->edgeNodes.size(); ++i)
    {
        for (size_t j = 0; j < n->edgeNodes.size(); ++j)
        {
            if (i == j) continue;
            if (n->edgeNodes[i]->index == n->index) continue;
            if (n->edgeNodes[j]->index == n->index) continue;

            if (areConnected(n->edgeNodes[i], n->edgeNodes[j]))
            {
                connectedEdgeNodes.push_back(n->edgeNodes[j]);
                connectedEdgeNodes.push_back(n->edgeNodes[i]);
                return true;
            }
        }
    }

    return false;
}

void constructPathFromBranchToSoma(SkeletonBranch* branch,
                                   SkeletonBranches& path,
                                   std::vector< size_t >& pathIndices)
{
    // Add the branch to the path
    path.push_back(branch);
    pathIndices.push_back(branch->index);
    branch->traversalCount += 1;

    // Ensure that the branch has a single parent
    if (branch->parents.size() == 1)
    {
        constructPathFromBranchToSoma(branch->parents[0], path, pathIndices);
    }
}

void getTerminals(SkeletonBranch* root, SkeletonBranches& terminals)
{
    if (root->children.size() == 0)
    {
        terminals.push_back(root);
        return;
    }
    else
    {
        for (size_t i = 0; i < root->children.size(); ++i)
        {
            getTerminals(root->children[i], terminals);
        }
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

void identifyTerminalBranchesForSpine(SkeletonBranches& branches)
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
        }
    }
}

SkeletonNode* identifyRootNodeForSpine(SkeletonBranch* rootBranch, const Vector3f basePoint)
{
    // The skeleton node to be identified later
    SkeletonNode* rootNode = nullptr;

    // The shortest distance between a terminal branch and the base point
    float shortestDistance = std::numeric_limits< float >::max();

    // Get the first and last sample of the branch
    auto firstSample = rootBranch->nodes.front();
    auto lastSample = rootBranch->nodes.back();

    // Compute the distance to the front and back nodes
    auto distanceToFirstSample = firstSample->point.distance(basePoint);
    auto distanceToLastSample = lastSample->point.distance(basePoint);

    if (distanceToFirstSample < shortestDistance)
    {
        shortestDistance = distanceToFirstSample;
        rootNode = firstSample;
    }

    if (distanceToLastSample < shortestDistance)
    {
        shortestDistance = distanceToLastSample;
        rootNode = lastSample;
    }

    return rootNode;
}

SkeletonBranch* identifyRootBranchForSpine(SkeletonBranches& branches, const Vector3f basePoint)
{
    // The skeleton branch to be identified later
    SkeletonBranch* rootBranch = nullptr;

    // The shortest distance between a terminal branch and the base point
    float shortestDistance = std::numeric_limits< float >::max();

    std::cout << branches.size() << "\n";
    // Select the terminal branch from the list of branches
    for (size_t i = 0; i < branches.size(); ++i)
    {
        // std::cout << i << " ";

        // Reference to the branch
        auto& branch = branches[i];

        // The branch must be terminal , otherwise next
        // if (!branch->isTerminal()) continue;

        // Invalid branches, next
        // if (!branch->isValid()) continue;

        // Get the first and last sample of the branch
        auto firstSample = branch->nodes.front();
        auto lastSample = branch->nodes.back();

        // Compute the distance to the front and back nodes
        auto distanceToFirstSample = firstSample->point.distance(basePoint);
        auto distanceToLastSample = lastSample->point.distance(basePoint);

        if (distanceToFirstSample < shortestDistance)
        {
            shortestDistance = distanceToFirstSample;
            rootBranch = branch;

        }

        if (distanceToLastSample < shortestDistance)
        {
            shortestDistance = distanceToLastSample;
            rootBranch = branch;
        }
    }

    return rootBranch;
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

void setIndexToSpineTree(SkeletonBranch *spineBranch, const size_t &spineIndex)
{
    spineBranch->spineIndex = spineIndex;
    spineBranch->setVisited();

    // Verify the rest of the tree
    for (size_t i = 0; i < spineBranch->logicalChildren.size(); ++i)
    {
        setIndexToSpineTree(spineBranch->logicalChildren[i], spineIndex);
    }
}

SkeletonBranch *verifySpineTree(SkeletonBranch *spine,
                                const size_t &spineIndex) {
    SkeletonBranch* nonSpinyParent = nullptr;

    // Get the parent (processed neurons will have indeed a single parent)
    auto parent = spine->logicalParents[0];
    if (parent->isSpine())
    {
        // Set the spine index
        parent->spineIndex = spineIndex;

        // Verify the tree of the parent to be spines or not
        nonSpinyParent = verifySpineTree(parent, spineIndex);
    }

    // Return the non-spiny parent to propagate backwards
    return nonSpinyParent;
}

SkeletonBranch* findEmanatingBranchOfSpine(SkeletonBranch *spine)
{
    // This is the branch, from which the spine emanates
    SkeletonBranch* spineEmanatingBranch = spine;

    // Get the parent (processed neurons will have indeed a single parent)
    auto parent = spine->logicalParents[0];
    if (parent->isSpine())
    {
        // Verify the tree of the parent to be spines or not
        spineEmanatingBranch = findEmanatingBranchOfSpine(parent);
    }

    // Return the non-spiny parent to propagate backwards
    return spineEmanatingBranch;
}

void segmentSpine(SkeletonBranch* spineTerminal, const size_t& spineIndex)
{
    // Set the spine index, and set it visited
    spineTerminal->spineIndex = spineIndex;
    spineTerminal->setVisited();

    // Verify the spine tree, and return the non-spiny parent
    auto nonSpinyParent = verifySpineTree(spineTerminal, spineIndex);

    // Propagate backwards
    if (nonSpinyParent != nullptr)
    {
        for (size_t i = 0; i < nonSpinyParent->logicalChildren.size(); ++i)
        {
            setIndexToSpineTree(nonSpinyParent->logicalChildren[i], spineIndex);
        }
    }
}


SkeletonBranch* segmentSpineFromBranch(SkeletonBranch* spineBranch, const size_t& spineIndex)
{
    SkeletonBranch* spineRoot;

    // If the spine branch does not have any parent
    if (spineBranch->logicalParents.size() == 0)
    {
        // Set the index to the spine branch
        spineBranch->spineIndex = spineIndex;

        // As this branch has no parents and it is a spine, then it is the spine root.
        spineRoot = spineBranch;

        // Traverse the forward tree to label the spine
        for (size_t i = 0; i < spineBranch->logicalChildren.size(); ++i)
        {
            setIndexToSpineTree(spineBranch->logicalChildren[i], spineIndex);
        }
    }

    // If the spine branch has a parent
    else
    {
        // Get this parent
        auto parent = spineBranch->logicalParents[0];

        // If the branch is a spine but the parent is not a spine
        if (!parent->isSpine())
        {
            spineRoot = spineBranch;

            // Traverse the forward tree to label the spine
            for (size_t i = 0; i < spineBranch->logicalChildren.size(); ++i)
            {
                setIndexToSpineTree(spineBranch->logicalChildren[i], spineIndex);
            }
        }

        // If the branch is a spine and the parent is also a spine
        else
        {
            // Find the parent that is not a spine
            auto emanatingBranch = findEmanatingBranchOfSpine(spineBranch);

            spineRoot = emanatingBranch;

            if (emanatingBranch != nullptr)
            {
                // Traverse the forward tree to label the spine
                for (size_t i = 0; i < emanatingBranch->logicalChildren.size(); ++i)
                {
                    setIndexToSpineTree(emanatingBranch->logicalChildren[i], spineIndex);
                }
            }
        }
    }

    // The spine branch is visited
    spineBranch->setVisited();

    return spineRoot;
}

void getTreeBoundingBox(SkeletonBranch* root,
                        Vector3f& pMin, Vector3f& pMax, Vector3f& bounds, Vector3f &center)
{
    // Get the bounding box of the root branch
    Vector3f pMinRoot, pMaxRoot, boundsRoot, centerRoot;
    root->getBoundingBox(pMinRoot, pMaxRoot, boundsRoot, centerRoot);

    for (size_t i = 0; i < root->logicalChildren.size(); ++i)
    {
        // Get the bounding box of the child branch
        Vector3f pMinChild, pMaxChild, boundsChild, centerChild;
        getTreeBoundingBox(root->logicalChildren[i], pMinChild, pMaxChild, boundsChild, centerChild);

        // Expandthe bounding box
        if (pMinChild.x() < pMinRoot.x()) pMinRoot.x() = pMinChild.x();
        if (pMinChild.y() < pMinRoot.y()) pMinRoot.y() = pMinChild.y();
        if (pMinChild.z() < pMinRoot.z()) pMinRoot.z() = pMinChild.z();

        if (pMaxChild.x() > pMaxRoot.x()) pMaxRoot.x() = pMaxChild.x();
        if (pMaxChild.y() > pMaxRoot.y()) pMaxRoot.y() = pMaxChild.y();
        if (pMaxChild.z() > pMaxRoot.z()) pMaxRoot.z() = pMaxChild.z();
    }

    // Calculate the center and bounds at the end
    pMin = pMinRoot;
    pMax = pMaxRoot;
    bounds = pMax - pMin;
    center = pMin + bounds * 0.5f;
}

void getLogicalTreeBoundingBox(SkeletonBranch* root,
                               Vector3f& pMin, Vector3f& pMax, Vector3f& bounds, Vector3f &center)
{
    // Get the bounding box of the root branch
    Vector3f pMinRoot, pMaxRoot, boundsRoot, centerRoot;
    root->getBoundingBox(pMinRoot, pMaxRoot, boundsRoot, centerRoot);

    for (size_t i = 0; i < root->logicalChildren.size(); ++i)
    {
        // Get the bounding box of the child branch
        Vector3f pMinChild, pMaxChild, boundsChild, centerChild;
        getTreeBoundingBox(root->logicalChildren[i], pMinChild, pMaxChild, boundsChild, centerChild);

        // Expandthe bounding box
        if (pMinChild.x() < pMinRoot.x()) pMinRoot.x() = pMinChild.x();
        if (pMinChild.y() < pMinRoot.y()) pMinRoot.y() = pMinChild.y();
        if (pMinChild.z() < pMinRoot.z()) pMinRoot.z() = pMinChild.z();

        if (pMaxChild.x() > pMaxRoot.x()) pMaxRoot.x() = pMaxChild.x();
        if (pMaxChild.y() > pMaxRoot.y()) pMaxRoot.y() = pMaxChild.y();
        if (pMaxChild.z() > pMaxRoot.z()) pMaxRoot.z() = pMaxChild.z();
    }

    // Calculate the center and bounds at the end
    pMin = pMinRoot;
    pMax = pMaxRoot;
    bounds = pMax - pMin;
    center = pMin + bounds * 0.5f;
}

}
