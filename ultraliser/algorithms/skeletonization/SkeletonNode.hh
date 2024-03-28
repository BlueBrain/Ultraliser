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

#include <common/Headers.hh>
#include <math/Vector3f.h>

namespace Ultraliser
{

/**
 * @brief The SkeletonNode struct
 * A node in the segmented skeleton.
 */
struct SkeletonNode
{
public:

    /**
     * @brief SkeletonNode
     */
    SkeletonNode()
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief SkeletonNode
     * @param index
     * @param point
     * @param radius
     */
    SkeletonNode(const size_t& index, const Vector3f& point, const float radius)
    {
        this->index = index;
        this->point = point;
        this->radius = radius;
    }

    /**
     * @brief SkeletonNode
     * @param index
     * @param point
     * @param voxel
     */
    SkeletonNode(const size_t& index, const Vector3f& point, const Vector3f& voxel)
    {
        this->index = index;
        this->point = point;
        this->voxel = voxel;
    }

    /**
     * @brief SkeletonNode
     * @param index
     * @param voxelIndex
     * @param point
     * @param voxel
     */
    SkeletonNode(const size_t& index, const int64_t voxelIndex,
                 const Vector3f& point, const Vector3f& voxel)
    {
        this->index = index;
        this->voxelIndex = voxelIndex;
        this->point = point;
        this->voxel = voxel;
    }

public:

    /**
     * @brief index
     * The index of the node.
     */
    size_t index = 0;

    /**
     * @brief voxelIndex
     * The 1D index of the corresponding voxel in the thinned volume.
     */
    int64_t voxelIndex = -1;

    /**
     * @brief orderIndex
     * An index that is used to order the nodes during the shortest path construction for
     * the simplified graph.
     */
    int64_t graphIndex = -1;

    /**
     * @brief swcIndex
     */
    int64_t swcIndex = -1;

    /**
     * @brief prevSampleSWCIndex
     */
    int64_t prevSampleSWCIndex = -1;

    /**
     * @brief point
     * The Cartesian coordinates of the node.
     */
    Vector3f point;

    /**
     * @brief voxel
     * The corresponding voxel [X, Y, Z] in the volume used to extract the skeleton.
     */
    Vector3f voxel;

    /**
     * @brief radius
     * The radius of the node;
     */
    float radius = 0.f;

    /**
     * @brief visited
     * If this flag is set, the node must have been visited before.
     * The flag is used to construct the connectivity of the skeleton.
     */
    bool visited = false;

    /**
     * @brief iVisit
     */
    size_t iVisit = 0;

    /**
     * @brief branching
     * If this flag is set, this node is a branching node, i.e. connected to more than two nodes.
     */
    bool branching = false;

    /**
     * @brief terminal
     * If this flag is set, this node is a terminal one, i.e. connected only to a single node.
     */
    bool terminal = false;

    /**
     * @brief insideSoma
     * If this flag is set, the node is located inside the soma and can be safely removed from
     * the graph.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool insideSoma = false;

    /**
     * @brief isSoma
     * If this flag is set, this node is a soma node.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool isSoma = false;

    /**
     * @brief connectedToSoma
     * If this flag is set, the node is connected to the soma, i.e. an initial node of a branch.
     * ONLY FOR NEURONS AND ASTROCYTES.
     */
    bool connectedToSoma = false;

    /**
     * @brief edgeNodes
     * A list of all the nodes connected along the edges to this node.
     */
    std::vector< SkeletonNode* > edgeNodes;
};

/**
 * @brief SkeletonNodes
 */
typedef std::vector< SkeletonNode* > SkeletonNodes;
}
