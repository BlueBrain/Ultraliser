/* Copyright (c) 2020, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of SimCrusher
 * <LINK>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#pragma once

#include <cstdint>
#include <vector>

namespace sc
{
class SparseOctreeNode
{
public:
    inline static const uint8_t SLOT_ROOT = 0;
    inline static const uint8_t SLOT_BACK_BOTTOM_LEFT = 1;
    inline static const uint8_t SLOT_BACK_BOTTOM_RIGHT = 2;
    inline static const uint8_t SLOT_BACK_TOP_LEFT = 4;
    inline static const uint8_t SLOT_BACK_TOP_RIGHT = 8;
    inline static const uint8_t SLOT_FRONT_BOTTOM_LEFT = 16;
    inline static const uint8_t SLOT_FRONT_BOTTOM_RIGHT = 32;
    inline static const uint8_t SLOT_FRONT_TOP_LEFT = 64;
    inline static const uint8_t SLOT_FRONT_TOP_RIGHT = 128;

public:
    SparseOctreeNode(const uint8_t value = 0, const uint8_t slotMask = SLOT_ROOT);
    
    /**
     * @brief Transform node into leaf node if it has 8 children
     */
    void compact();

    /**
     * @brief Returns this node children list
     * 
     * @return const std::vector<SparseOctreeNode>& 
     */
    const std::vector<SparseOctreeNode>& getChildren() const;

    /**
     * @brief Returns the number of children on this node
     * 
     * @return size_t 
     */
    size_t getNumChildren() const;

    /**
     * @brief Return the children mask. The children mask is 1 byte flag indicating which children are present.
     * 
     * @return uint8_t 
     */
    uint8_t getChildrenMask() const;
    
    /**
     * @brief Adds a new children node.
     * 
     * @param node 
     */
    void addChildNode(const SparseOctreeNode& node);

    /**
     * @brief Returns the node slot mask. The slot mask is bit flag that indicates the position relative to the 
     * parent node.
     * 
     * @return uint8_t 
     */
    uint8_t getSlotMask() const;

private:
    uint8_t _slotMask;
    uint8_t _childMask;
    std::vector<SparseOctreeNode> _children;

    friend class SparseOctree;
};
}