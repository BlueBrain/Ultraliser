/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Nadir Roman Guerrero < nadir.romanguerrero@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include "Bounds.h"

#include <cstdint>
#include <vector>

namespace Ultraliser
{
class SparseOctreeNode
{
public:
    SparseOctreeNode() = default;

    explicit SparseOctreeNode(uint8_t slotMask);

    /**
     * @brief Tries to compact this node by emptying the children nodes and storing a full children mask.
     *
     * @return true If the node was full and thus compacted, false otherwise.
     */
    bool tryCompact();

    /**
     * @brief Return a specific children of this node
     *
     * @param index Index of the node in the children list. Access bounds are not checked.
     * @return SparseOctreeNode&
     */
    SparseOctreeNode &getChild(size_t index);
    const SparseOctreeNode &getChild(size_t index) const;

    /**
     * @brief Tries to find a children by the slot it occupies.
     *
     * @param slot The slot of the children
     * @return SparseOctreeNode* if found, nullptr otherwise
     */
    SparseOctreeNode *findChild(uint8_t slot);
    const SparseOctreeNode *findChild(uint8_t slot) const;

    /**
     * @brief Returns the number of children on this node
     *
     * @return size_t
     */
    size_t getNumChildren() const;

    /**
     * @brief Return the children mask. The children mask is 1 byte flag indicating which child slots are present.
     *
     * @return uint8_t
     */
    uint8_t getChildrenMask() const;

    /**
     * @brief Adds a new children node. If the node already has a children on that slot, this function has no effect.
     *
     * @param slot The slot the node will occupy.
     */
    void addChildNode(uint8_t slot);

    /**
     * @brief Adds a new children leaf node (It appears on its children mask, but the node itself is not created). If
     * the node already has a children on that slot, this function has no effect.
     *
     * @param slot The slot the node will occupy.
     */
    void addChildLeaf(uint8_t slot);

    /**
     * @brief Returns the node slot mask. The slot mask is bit flag that indicates the position relative to the
     * parent node.
     *
     * @return uint8_t
     */
    uint8_t getSlotMask() const;

private:
    uint8_t _slotMask = 0;
    uint8_t _childMask;
    std::vector<SparseOctreeNode> _children;
};
}