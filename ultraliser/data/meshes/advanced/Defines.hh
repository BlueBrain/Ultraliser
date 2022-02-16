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

#ifndef ULTRALISER_DATA_MESHES_ADVANCED_DEFINES_HH
#define ULTRALISER_DATA_MESHES_ADVANCED_DEFINES_HH

#define FABS(a) (((a)<0)?(-(a)):(a))
#define LOG2(a) (log(a)/log(2))
#define PI2	(M_PI/2.0)
#ifndef MIN
#define MIN(a,b) (((a)<(b))?(a):(b))
#endif
#ifndef MAX
#define MAX(a,b)(((a)>(b))?(a):(b))
#endif

// Casting
#define NODE_TO_DOUBLE(NUMBER)         (static_cast< double >(NUMBER))
#define NODE_TO_FLOAT(NUMBER)          (static_cast< float >(NUMBER))
#define NODE_TO_INT(NUMBER)            (static_cast< int >(NUMBER))

// Getting vertices of a triangle
#define VERTEX_1(MESH)                  (NODE_TO_INT(((AdvancedTriangle *) MESH->data)->v1()->x))
#define VERTEX_2(MESH)                  (NODE_TO_INT(((AdvancedTriangle *) MESH->data)->v2()->x))
#define VERTEX_3(MESH)                  (NODE_TO_INT(((AdvancedTriangle *) MESH->data)->v3()->x))

#define DETERMINANT3X3(A11, A12, A13, A21, A22, A23, A31, A32, A33)                                 \
    ((A11) * ((A22) * (A33) - (A23) * (A32)) -                                                      \
     (A12) * ((A21) * (A33) - (A23) * (A31)) +                                                      \
     (A13) * ((A21) * (A32) - (A22) * (A31)))


#endif // ULTRALISER_DATA_MESHES_ADVANCED_DEFINES_HH
