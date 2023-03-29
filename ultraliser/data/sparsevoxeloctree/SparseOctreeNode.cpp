/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
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

#include "SparseOctreeNode.h"

#include <algorithm>
#include <cassert>

namespace
{
class ChildFinder
{
public:
    template<typename Container>
    static auto find(Container &&children, uint8_t childrenMask, uint8_t childSlot)
    {
        if (!(childrenMask & childSlot))
        {
            return static_cast<decltype(children.data())>(nullptr);
        }

        uint8_t i = 1;
        size_t index = 0;
        while (childSlot > i)
        {
            if (i & childrenMask)
            {
                ++index;
            }
            i <<= 1;
        }
        return &children[index];
    }
};

class CompacterVerifier
{
public:
    static bool isCompactable(const Ultraliser::SparseOctreeNode &node)
    {
        if (node.getNumChildren() == 8)
        {
            for (size_t i = 0; i < 8; ++i)
            {
                auto &child = node.getChild(i);
                if (child.getNumChildren() > 0 || child.getChildrenMask() != 0)
                {
                    return false;
                }
            }
            return true;
        }

        if (node.getNumChildren() == 0 && node.getChildrenMask() == 0xff)
        {
            return true;
        }

        return false;
    }
};

}

namespace Ultraliser
{
SparseOctreeNode::SparseOctreeNode(uint8_t slotMask)
    : _slotMask(slotMask)
    , _childMask(0)
{
}

void SparseOctreeNode::compact()
{
    if (!CompacterVerifier::isCompactable(*this))
    {
        return;
    }

    _children.clear();
    _childMask = 0xff;
}

SparseOctreeNode &SparseOctreeNode::getChild(size_t index)
{
    assert(index < _children.size());
    return _children[index];
}

const SparseOctreeNode &SparseOctreeNode::getChild(size_t index) const
{
    assert(index < _children.size());
    return _children[index];
}

SparseOctreeNode *SparseOctreeNode::findChild(uint8_t slot)
{
    return ChildFinder::find(_children, _childMask, slot);
}

const SparseOctreeNode *SparseOctreeNode::findChild(uint8_t slot) const
{
    return ChildFinder::find(_children, _childMask, slot);
}

size_t SparseOctreeNode::getNumChildren() const
{
    return _children.size();
}

uint8_t SparseOctreeNode::getChildrenMask() const
{
    return _childMask;
}

void SparseOctreeNode::addChildNode(uint8_t slot)
{
    if ((slot & _childMask))
    {
        return;
    }

    _childMask |= slot;
    _children.emplace_back(slot);

    // Sort children based on slot mask to keep access ordering correct
    std::sort(_children.begin(), _children.end(), [](auto &a, auto &b) { return a.getSlotMask() < b.getSlotMask(); });
}

void SparseOctreeNode::addChildLeaf(uint8_t slot)
{
    _childMask |= slot;
}

uint8_t SparseOctreeNode::getSlotMask() const
{
    return _slotMask;
}
}