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

#include <algorithms/skeletonization/SkeletonNode.hh>
#include <algorithms/skeletonization/SkeletonEdge.hh>
#include <algorithms/skeletonization/SkeletonBranch.h>

namespace Ultraliser
{

/**
 * @brief removeEdgeNode
 * @param node
 * @param edgeNodeMarkedForRemoval
 */
void removeEdgeNode(SkeletonNode* node, SkeletonNode* edgeNodeMarkedForRemoval);

/**
 * @brief collapseTriangleIntoNode
 * @param nodes
 * @param n1
 * @param n2
 * @param n3
 */
void collapseTriangleIntoNode(SkeletonNodes& nodes,
                              SkeletonNode* n1, SkeletonNode* n2, SkeletonNode* n3);

/**
 * @brief areConnected
 * @param n1
 * @param n2
 * @return
 */
bool areConnected(const SkeletonNode* n1, const SkeletonNode* n2);

/**
 * @brief isTriangleNode
 * @param n
 * @param connectedEdgeNodes
 * @return
 */
bool isTriangleNode(const SkeletonNode* n, SkeletonNodes& connectedEdgeNodes);

/**
 * @brief updateTraversalCountAlongPathFromBranchToSoma
 * @param branch
 * @param path
 * @param pathIndices
 */
void updateTraversalCountAlongPathFromBranchToSoma(SkeletonBranch* branch,
                                                   SkeletonBranches& path,
                                                   std::vector< size_t >& pathIndices);

/**
 * @brief getTerminals
 * @param root
 * @param terminals
 */
void getTerminals(SkeletonBranch* root, SkeletonBranches& terminals);

/**
 * @brief _detectChildren
 * @param currentBranch
 * @param allBranches
 * @return
 */
SkeletonBranches _detectChildren(SkeletonBranch* currentBranch, SkeletonBranches& allBranches);

/**
 * @brief identifyRootNodeForSpine
 * @param rootBranch
 * @param basePoint
 * @return
 */
SkeletonNode* identifyRootNodeForSpine(SkeletonBranch* rootBranch, const Vector3f basePoint);

/**
 * @brief identifyTerminalConnections
 * @param branches
 */
void identifyTerminalConnections(SkeletonBranches& branches);

/**
 * @brief identifyTerminalBranchesForSpine
 * @param branches
 */
void identifyTerminalBranchesForSpine(SkeletonBranches& branches);

/**
 * @brief identifyRootBranchForSpine
 * @param branches
 * @param basePoint
 * @return
 */
SkeletonBranch* identifyRootBranchForSpine(SkeletonBranches& branches,
                                           const Vector3f basePoint);

/**
 * @brief confirmTerminalsBranches
 * @param branches
 */
void confirmTerminalsBranches(SkeletonBranches& branches);

/**
 * @brief segmentSpine
 * @param spineTerminal
 * @param spineIndex
 */
void segmentSpine(SkeletonBranch* spineTerminal, const size_t& spineIndex);

SkeletonBranch* segmentSpineFromBranch(SkeletonBranch* spineBranch, const size_t& spineIndex);

void getLogicalTreeBoundingBox(SkeletonBranch* root,
                               Vector3f& pMin, Vector3f& pMax, Vector3f& bounds, Vector3f &center);

}
