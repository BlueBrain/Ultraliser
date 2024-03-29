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

#include <data/meshes/advanced/AdvancedMesh.h>

namespace Ultraliser
{

#define DI_MAX_NUMBER_OF_CELLS      10000
#define DI_EPSILON_POINT            AdvancedPoint(1.0e-9, 1.0e-9, 1.0e-9)

/**
 * @brief The IntersectionCell class
 */
class IntersectionCell
{
public:

    /**
     * @brief pMin
     */
    AdvancedPoint pMin;

    /**
     * @brief pMax
     */
    AdvancedPoint pMax;

    /**
     * @brief triangles
     */
    List triangles;

    /**
     * @brief IntersectionCell
     */
    IntersectionCell() { }

    /**
     * @brief IntersectionCell
     * @param input
     * @param useAll
     */
    IntersectionCell(AdvancedMesh *input, bool useAll = true);

    /**
     * @brief isTriangleBoundingBoxInCell
     * @param triangle
     * @return
     */
    bool isTriangleBoundingBoxInCell(AdvancedTriangle *triangle) const;

    /**
     * @brief fork
     * @return
     */
    IntersectionCell *fork();

    /**
     * @brief selectIntersections
     * @param justProper
     */
    void selectIntersections(bool justProper = false);

    /**
     * @brief doesNotIntersectForSure
     * @return
     */
    bool doesNotIntersectForSure();
};

} 

