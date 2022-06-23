/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marco Attene < IMATI-GE / CNR >
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
 *
 * The content of this file is based on MeshFix. The code has been modified under the terms of
 * the GNU General Public License as published by the Free Software Foundation either version 3 of
 * the License, or (at your option) any later version.
 * MeshFix has a dual license for free and commercial use. For further information, please refer
 * to the original repository at < https://github.com/MarcoAttene/MeshFix-V2.1>.
 **************************************************************************************************/

#pragma once

#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/primitives/AdvancedEdge.h>
#include <data/meshes/advanced/primitives/AdvancedPoint.h>
#include <data/meshes/advanced/primitives/AdvancedTriangle.h>

namespace Ultraliser
{

/**
 * @brief The AdvancedLoopSplit class
 * This LoopSplit is part of the AdvancedMesh.
 * AdvancedMesh is based on TMesh of MeshFix, which preserves the connectivity information.
 */
class AdvancedLoopSplit
{
public:

    /**
     * @brief edge
     */
    AdvancedEdge *edge;

    /**
     * @brief point
     */
    AdvancedPoint point;

    /**
     * @brief LoopSplit
     * Constructor
     *
     * @param inputEdge
     * @param md
     */
    AdvancedLoopSplit(AdvancedEdge *inputEdge, int md)
    {
        edge = inputEdge;

        // If the edge is not part of a triangle
        if (edge->t1 == nullptr || edge->t2 == nullptr || md)
        {
            point = ((*edge->v1) + (*edge->v2)) * 0.5f;
        }
        else
        {
            AdvancedVertex *ov1 = edge->t1->oppositeVertex(edge);
            AdvancedVertex *ov2 = edge->t2->oppositeVertex(edge);

            point = ((((*edge->v1) + (*edge->v2)) * 3.0f) + ((*ov1) + (*ov2))) / 8.0f;
        }
    }
};

}

