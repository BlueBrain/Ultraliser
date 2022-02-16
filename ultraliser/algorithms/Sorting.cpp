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
 * This file has been adapted from MeshFix under the terms of the GNU General Public License as
 * published by the Free Software Foundation of version 3. TMesh has a dual license for free and
 * commercial use. For further information, please contact Marco Attene.
 **************************************************************************************************/

#define USE_STD_SORT

#ifdef USE_STD_SORT
#include <vector>
#include <algorithm>
#endif

namespace Ultraliser
{

class ComparisonObject
{
    int(*compare)(const void *, const void *);

public:
    ComparisonObject(int(*c)(const void *, const void *))
    {
        compare = c;
    }

    bool operator()(void *a, void *b)
    {
        return (compare(a, b) < 0);
    }
};

void jqSort(void *v[], int numberElements, int(*comparisonFunction)(const void *, const void *))
{
    ComparisonObject a(comparisonFunction);
    std::sort(v, v + numberElements, a);
}

} 

