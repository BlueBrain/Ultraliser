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
struct SparseOctreeNodeSlot
{
    static constexpr uint8_t root = 0;
    static constexpr const uint8_t backBottomLeft = 1;
    static constexpr const uint8_t backBottomRight = 2;
    static constexpr const uint8_t backTopLeft = 4;
    static constexpr const uint8_t backTopRight = 8;
    static constexpr const uint8_t frontBottomLeft = 16;
    static constexpr const uint8_t frontBottomRight = 32;
    static constexpr const uint8_t frontTopLeft = 64;
    static constexpr const uint8_t frontTopRight = 128;
};

class SparseOctreeNode
{
public:
    SparseOctreeNode(uint8_t slotMask);

    /**
     * @brief Transform node into leaf node if it has 8 children
     */
    void compact();

    /**
     * @brief Return a specific children of this node
     *
     * @param index Index of the node in the children list
     * @return SparseOctreeNode&
     */
    SparseOctreeNode &getChild(size_t index);
    const SparseOctreeNode &getChild(size_t index) const;

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
    void addChildNode(uint8_t slot);

    /**
     * @brief Returns the node slot mask. The slot mask is bit flag that indicates the position relative to the
     * parent node.
     *
     * @return uint8_t
     */
    uint8_t getSlotMask() const;

    /**
     * @brief Returns wether this node is a leaf node.
     *
     * @return true if the node has no children, false otherwise.
     */
    bool isLeaf() const;

private:
    friend class SparseOctree;

    uint8_t _slotMask;
    uint8_t _childMask;
    std::vector<SparseOctreeNode> _children;
};
}