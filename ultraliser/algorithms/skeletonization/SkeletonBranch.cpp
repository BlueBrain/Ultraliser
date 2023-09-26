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

#define VALIDAITY_BIT_INDEX     0
#define VISITED_BIT_INDEX       1
#define TERMINAL_BIT_INDEX      2
#define ROOT_BIT_INDEX          3
#define DUPLICATE_BIT_INDEX     4


#include "SkeletonBranch.h"

namespace Ultraliser
{

SkeletonBranch::SkeletonBranch()
{
    _flags = new BitArray(8);
    _flags->clearAll();
}

SkeletonBranch::~SkeletonBranch()
{
    _flags->~BitArray();
}

void SkeletonBranch::printTree(size_t order)
{
    std::cout << std::string(order * 4, '-') << index << "\n";
    for (size_t i = 0; i < children.size(); ++i)
    {
        children[i]->printTree(order + 1);
    }
}

bool SkeletonBranch::hasTerminalNodes(const size_t& node1Index, const size_t& node2Index)
{
    const auto& frontNode = nodes.front();
    const auto& backNode = nodes.back();

    if (node1Index == backNode->index && node2Index == frontNode->index)
    {
        return true;
    }

    if (node1Index == frontNode->index && node2Index == backNode->index)
    {
        return true;
    }

    return false;
}

void SkeletonBranch::adjustDirection(const size_t& frontNodeIndex, const size_t& backNodeIndex)
{
    if (nodes.front()->index == backNodeIndex && nodes.back()->index == frontNodeIndex)
    {
        std::reverse(this->nodes.begin(), this->nodes.end());
    }
}

bool SkeletonBranch::isLoop() const
{
    return this->nodes.front()->index == this->nodes.back()->index;
}

void SkeletonBranch::setValid()
{
    _flags->setBit(VALIDAITY_BIT_INDEX);
}

void SkeletonBranch::setInvalid()
{
    _flags->clearBit(VALIDAITY_BIT_INDEX);
}

bool SkeletonBranch::isValid() const
{
    return _flags->bit(VALIDAITY_BIT_INDEX);
}

void SkeletonBranch::setVisited()
{
    _flags->setBit(VISITED_BIT_INDEX);
}

void SkeletonBranch::setUnvisited()
{
     _flags->clearBit(VISITED_BIT_INDEX);
}

bool SkeletonBranch::visited() const
{
    return _flags->bit(VISITED_BIT_INDEX);
}

void SkeletonBranch::setTerminal()
{
    _flags->setBit(TERMINAL_BIT_INDEX);
}

void SkeletonBranch::unsetTermainal()
{
    _flags->clearBit(TERMINAL_BIT_INDEX);
}

bool SkeletonBranch::isTerminal() const
{
    return _flags->bit(TERMINAL_BIT_INDEX);
}

void SkeletonBranch::setRoot()
{
    _flags->setBit(ROOT_BIT_INDEX);
}

void SkeletonBranch::unsetRoot()
{
    _flags->clearBit(ROOT_BIT_INDEX);
}

bool SkeletonBranch::isRoot() const
{
    return _flags->bit(ROOT_BIT_INDEX);
}

void SkeletonBranch::setDuplicate()
{
    _flags->setBit(DUPLICATE_BIT_INDEX);
}

void SkeletonBranch::unsetDuplicate()
{
    _flags->clearBit(DUPLICATE_BIT_INDEX);
}

bool SkeletonBranch::isDuplicate() const
{
    return _flags->bit(DUPLICATE_BIT_INDEX);
}

float SkeletonBranch::computeLength() const
{
    float length = 0;
    for (size_t i = 0; i < nodes.size() - 1; ++i)
    {
        auto p1 = nodes[i]->point;
        auto p2 = nodes[i + 1]->point;
        length += p1.distance(p2);
    }

    return length;
}

}
