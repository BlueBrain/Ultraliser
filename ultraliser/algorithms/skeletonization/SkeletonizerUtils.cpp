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

#include "SkeletonizerUtils.h"

namespace Ultraliser
{

void removeEdgeNode(SkeletonNode* node, SkeletonNode* edgeNodeMarkedForRemoval)
{
    auto it = std::find(node->edgeNodes.begin(), node->edgeNodes.end(), edgeNodeMarkedForRemoval);

    if(it != node->edgeNodes.end())
        node->edgeNodes.erase(it);
    node->edgeNodes.shrink_to_fit();
}

void collapseTriangleIntoNode(SkeletonNodes& nodes,
                              SkeletonNode* n1, SkeletonNode* n2, SkeletonNode* n3)
{
    SkeletonNode* centerNode = new SkeletonNode();
    centerNode->point = (n1->point + n2->point + n3->point) / 3.f;
    centerNode->radius = (n1->radius + n2->radius + n3->radius) / 3.f;

    centerNode->index = nodes.back()->index + 1;
    centerNode->branching = true;
    centerNode->terminal = false;

    removeEdgeNode(n1, n2); removeEdgeNode(n1, n3);
    removeEdgeNode(n2, n1); removeEdgeNode(n2, n3);
    removeEdgeNode(n3, n1); removeEdgeNode(n3, n2);

    centerNode->edgeNodes.push_back(n1);
    centerNode->edgeNodes.push_back(n2);
    centerNode->edgeNodes.push_back(n3);

    n1->edgeNodes.push_back(centerNode);
    n2->edgeNodes.push_back(centerNode);
    n3->edgeNodes.push_back(centerNode);

    nodes.push_back(centerNode);
}

bool areConnected(const SkeletonNode* n1, const SkeletonNode* n2)
{
    for (size_t i = 0; i < n1->edgeNodes.size(); ++i)
    {
        if (n1->edgeNodes[i]->index == n2->index)
        {
            return true;
        }
    }

    return false;
}

bool isTriangleNode(const SkeletonNode* n, SkeletonNodes& connectedEdgeNodes)
{
    for (size_t i = 0; i < n->edgeNodes.size(); ++i)
    {
        for (size_t j = 0; j < n->edgeNodes.size(); ++j)
        {
            if (i == j) continue;
            if (n->edgeNodes[i]->index == n->index) continue;
            if (n->edgeNodes[j]->index == n->index) continue;

            if (areConnected(n->edgeNodes[i], n->edgeNodes[j]))
            {
                connectedEdgeNodes.push_back(n->edgeNodes[j]);
                connectedEdgeNodes.push_back(n->edgeNodes[i]);
                return true;
            }
        }
    }

    return false;
}

}
