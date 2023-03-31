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

#include "OctreeStatsPrinter.h"

#include <iostream>
#include <queue>
#include <vector>

namespace
{
class MemorySize
{
public:
    inline static std::vector<std::string_view> units = {"b", "Kb", "Mb", "Gb"};

    static std::string reduceUnit(size_t bytes)
    {
        size_t unitIndex = 0;
        auto pivot = static_cast<double>(bytes);

        constexpr auto invUnitStep = 1.0 / 1024.0;

        while (unitIndex + 1 < units.size() && pivot > 1024.0)
        {
            pivot *= invUnitStep;
            ++unitIndex;
        }

        return std::to_string(pivot) + " " + std::string(units[unitIndex]);
    }
};
}

namespace Ultraliser
{
void OctreeStatsPrinter::print(const SparseOctree &octree)
{
    size_t numIntermediateNodes = 0;
    size_t numLeafNodes = 0;

    auto queue = std::queue<const SparseOctreeNode *>();
    queue.push(&octree.getRoot());

    auto leafNodeHolders = std::vector<const SparseOctreeNode *>();

    while (!queue.empty())
    {
        auto next = queue.front();
        queue.pop();

        ++numIntermediateNodes;

        if (next->getNumChildren() == 0)
        {
            leafNodeHolders.push_back(next);
            continue;
        }

        for (size_t i = 0; i < next->getNumChildren(); ++i)
        {
            queue.push(&next->getChild(i));
        }
    }

    for (auto holder : leafNodeHolders)
    {
        auto mask = holder->getChildrenMask();
        for (size_t i = 0; i < 8; ++i)
        {
            if (mask & (1 << i))
            {
                ++numLeafNodes;
            }
        }
    }

    std::cout << "Octree stats:\n";
    std::cout << "\tMax depth: \t" << octree.getMaxDepth() << "\n";
    std::cout << "\tGrid resolution: \t" << (1u << octree.getMaxDepth()) << "\n";
    std::cout << "\tIntermediate nodes: \t" << numIntermediateNodes << "\n";
    std::cout << "\tLeaf nodes: \t" << numLeafNodes << std::endl;
    std::cout << "\tTotal memory: " << MemorySize::reduceUnit(numIntermediateNodes * sizeof(SparseOctreeNode)) << "\n";
    std::cout << std::endl;
}
}