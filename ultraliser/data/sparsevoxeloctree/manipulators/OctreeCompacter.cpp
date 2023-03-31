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

#include "OctreeCompacter.h"

#include <queue>
#include <stack>

namespace
{
class NodeCompacter
{
public:
    struct Compactable
    {
        Ultraliser::SparseOctreeNode *node;
        Compactable *parent = nullptr;
        bool compactable = true;
    };

    static void compact(Ultraliser::SparseOctree &octree)
    {
        auto compactCandidates = std::stack<Compactable>();

        auto searchQueue = std::queue<Compactable>();
        searchQueue.push({&octree.getRoot()});

        while (!searchQueue.empty())
        {
            auto compactable = searchQueue.front();
            searchQueue.pop();

            compactCandidates.push(compactable);
            auto parent = &compactCandidates.top();

            auto node = compactable.node;
            for (size_t i = 0; i < node->getNumChildren(); ++i)
            {
                searchQueue.push({&node->getChild(i), parent});
            }
        }

        while (!compactCandidates.empty())
        {
            auto candidate = compactCandidates.top();
            compactCandidates.pop();

            auto parent = candidate.parent;
            if (!candidate.compactable)
            {
                if (parent)
                {
                    parent->compactable = false;
                }
                continue;
            }

            auto node = candidate.node;
            if (!node->tryCompact())
            {
                parent->compactable = false;
            }
        }
    }
};
}

namespace Ultraliser
{
void OctreeCompacter::compact(SparseOctree &octree)
{
    NodeCompacter::compact(octree);
}
}