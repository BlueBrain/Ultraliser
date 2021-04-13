/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_INDEXED_TRIANGLE_HH
#define ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_INDEXED_TRIANGLE_HH

#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/primitives/AdvancedTriangle.h>

namespace Ultraliser
{

/**
 * @brief The AdvancedIndexedTriangle struct
 * This Triangle is part of the AdvancedMesh.
 * AdvancedMesh is based on TMesh of MeshFix, which preserves the connectivity information.
 */
struct AdvancedIndexedTriangle
{
public:

    /**
     * @brief AdvancedIndexedTriangle
     * Constructor
     */
    AdvancedIndexedTriangle()
        : _tri(nullptr)
        , _v1(0)
        , _v2(0)
        , _v3(0)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief AdvancedIndexedTriangle
     * Constructor
     *
     * @param input
     * Input triangle.
     * @param v1
     * First vertex
     * @param v2
     * Second vertex
     * @param v3
     * Third vertex
     */
    AdvancedIndexedTriangle(const AdvancedTriangle *input,
                    uint64_t v1,
                    uint64_t v2,
                    uint64_t v3)
        : _tri(input)
        , _v1(v1)
        , _v2(v2)
        , _v3(v3)
    {
        /// EMPTY CONSTRUCTOR
    }

    /**
     * @brief containsVertex
     * Checks if this AdvancedIndexedTriangle contains a vertex with a given vertex index or not.
     * @param vIndex
     * The index of the vertex we search.
     * @return
     * True or false
     */
    bool containsVertex(const uint64_t vIndex) const
    {
        return (_v1 == vIndex || _v2 == vIndex || _v3 == vIndex);
    }

public:

    /**
     * @brief _tri
     * Simply, an AdvancedTriangle
     */
    const AdvancedTriangle* _tri;

    /**
     * @brief _v1
     * First vertex
     */
    const uint64_t _v1;

    /**
     * @brief _v2
     * Second vertex
     */
    const uint64_t _v2;

    /**
     * @brief _v3
     * Third vertex
     */
    const uint64_t _v3;
};

}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_PRIMITIVES_ADVANCED_INDEXED_TRIANGLE_HH
