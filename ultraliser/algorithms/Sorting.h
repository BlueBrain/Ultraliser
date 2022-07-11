/***************************************************************************************************
 * Copyright (c) 2013 - 2021
 * Consiglio Nazionale delle Ricerche, Sezione di Genova, IMATI-GE / CNR
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
 * This file has been adapted from TMesh under the terms of the GNU General Public License as
 * published by the Free Software Foundation of version 3. TMesh has a dual license for free and
 * commercial use. For further information, please contact Marco Attene.
 **************************************************************************************************/

#pragma once

namespace Ultraliser
{

/**
 * @brief jqSort
 * A generic QuickSort function.
 * The jqSort() function sorts an array with a given number of elements @numberElements.
 * The v argument points to the start of the array of elements casted to void *.
 * The contents of the array are sorted in ascending order according to a comparison function
 * pointed to by comp, which is called with two arguments that point to the objects being compared.
 * The comparison function must return an integer less than, equal to, or greater than zero if the
 * first argument is considered to be respectively less than, equal to, or greater than the second.
 * If two members compare as equal, their order in the sorted array is undefined.
 * See the manpage of the standard library qsort() function for further information.
 *
 * @param v
 * The v argument points to the start of the array of elements casted to void *.
 * @param numberElements
 * Number of elements in the array.
 */
extern void jqSort(void *v[], int numberElements,
                   int (*comparisonFunction)(const void *, const void *));

} 
