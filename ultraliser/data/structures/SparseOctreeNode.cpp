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

#include "SparseOctreeNode.h"

#include <algorithm>

namespace sc
{       
SparseOctreeNode::SparseOctreeNode(const uint8_t value, 
                                    const uint8_t slotMask)
  : _nodeValue(value)
  , _slotMask(slotMask)
  , _childMask(0)
  , _weight(1)
{
}
    
/* 
  * Removes all children if they share the same value 
  */
void SparseOctreeNode::compact()
{
  // Do not compact if the node is not full
  if(_children.size() < 8)
    return;

  bool compact = true;
  for(const auto& child : _children)
  {
    // If a children has a different value (or has children itself)
    // do not compact
    if(child.getValue() != _nodeValue 
      || child.getChildrenMask() != 0)
    {
      compact = false;
      break;
    }
  }

  if(compact)
  {
    _children.clear();
    _childMask = 0;
  }
}

/* 
  * Return a list with the current node children 
  */
std::vector<SparseOctreeNode>& SparseOctreeNode::getChildren()
{
  return _children;
}

/* 
  * Return a list with the current node children 
  */
const std::vector<SparseOctreeNode> & SparseOctreeNode::getChildren() const
{
  return _children;
}

/* 
  * Return the current number of child nodes 
  */
size_t SparseOctreeNode::getNumChildren() const
{
  return _children.size();
}

/* 
  * Return the children mask (aggregated slot mask) of the child nodes
  */
uint8_t SparseOctreeNode::getChildrenMask() const
{
  return _childMask;
}

/* 
  * Adds a children node if its slot mask is not already 
  * pressent in the children mask
  */
void SparseOctreeNode::addChildNode(const SparseOctreeNode& node)
{
  // If the child is not present already
  if(!(node.getSlotMask() & _childMask))
  {
    _childMask |= node.getSlotMask();
    _children.push_back(node);

    // Sort children based on slot mask to keep access ordering
    // correct
    std::sort(_children.begin(), _children.end(), 
              [](SparseOctreeNode& a, SparseOctreeNode& b)
              {
                return a.getSlotMask() < b.getSlotMask();
              });
  }
}

/* 
  * Returns this node slot mask relative to its parent 
  */
uint8_t SparseOctreeNode::getSlotMask() const
{
  return _slotMask;
}

/*
  * Sets the node/voxel value
  */
void SparseOctreeNode::setValue(const uint8_t value)
{
  _nodeValue = value;
  _weight = 1;
}

void SparseOctreeNode::addValue(const uint8_t value)
{
  // Weighted current sample
  const float prevWeightedValue = static_cast<float>(_weight * _nodeValue);

  //  + 0.5 to round the number 
  const float newValue = (prevWeightedValue + static_cast<float>(value)) / static_cast<float>(_weight + 1) + 0.5;

  _nodeValue = static_cast<uint8_t>(newValue); 
  ++_weight;
}

uint8_t SparseOctreeNode::getValue() const
{
  return _nodeValue;
}
}