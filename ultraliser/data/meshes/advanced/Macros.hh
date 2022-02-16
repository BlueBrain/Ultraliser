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

#ifndef ULTRALISER_MESH_ADVANCED_MACROS_HH
#define ULTRALISER_MESH_ADVANCED_MACROS_HH

#define TO_TRIANGLE(POINTER)    static_cast< AdvancedTriangle* > (POINTER)

/// For each vertex
#define FOR_EACH_VERTEX(Tt, n)                                                                      \
    for (n = _vertices.head(), Tt = (n) ? (TO_VERTEX(n->data)) : nullptr;                           \
         n != nullptr; n = n->next(), Tt = (n) ? (TO_VERTEX(n->data)) : nullptr)

/// For each edge
#define FOR_EACH_EDGE(Tt, n)                                                                        \
    for (n = _edges.head(), Tt = (n) ? (TO_EDGE(n->data)) : nullptr;                                \
         n != nullptr; n = n->next(), Tt = (n) ? (TO_EDGE(n->data)) : nullptr)

/// For each triangle
#define FOR_EACH_TRIANGLE(Tt, n)                                                                    \
    for (n = _triangles.head(), Tt = (n) ? (TO_TRIANGLE(n->data)) : nullptr;                        \
         n != nullptr; n = n->next(), Tt = (n) ? (TO_TRIANGLE(n->data)) : nullptr)

#define TO_UCHAR(VALUE)                 static_cast< unsigned char >(VALUE)
#define MARK_VISIT(A)                   ((A)->mask |= (TO_UCHAR(1)))
#define IS_VISITED(A)                   ((A)->mask &  (TO_UCHAR(1)))
#define UNMARK_VISIT(A)                 ((A)->mask &= (~(TO_UCHAR(1))))

#define MARK_VISIT2(A)                  ((A)->mask |= (TO_UCHAR(2)))
#define IS_VISITED2(A)                  ((A)->mask &  (TO_UCHAR(2)))
#define UNMARK_VISIT2(A)                ((A)->mask &= (~(TO_UCHAR(2))))

#define MARK_BIT(A, B)                  ((A)->mask |= (TO_UCHAR((1 << B))))
#define IS_BIT(A, B)                    ((A)->mask &  (TO_UCHAR((1 << B))))
#define UNMARK_BIT(A, B)                ((A)->mask &= (~(TO_UCHAR((1 << B)))))

#define TAG_SHARPEDGE(A)                (MARK_BIT((A),7))
#define IS_SHARPEDGE(A)                 (IS_BIT((A),7))
#define UNTAG_SHARPEDGE(A)              (UNMARK_BIT((A),7))

#endif // ULTRALISER_MESH_ADVANCED_MACROS_HH
