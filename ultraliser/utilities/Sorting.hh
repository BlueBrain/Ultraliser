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

#pragma once

#include <common/Headers.hh>

namespace Ultraliser
{

/**
 * @brief secondComparator
 *  Comparator function to sort pairs according to second value.
 * @param a
 * First element.
 * @param b
 * Second element.
 * @return
 *  The largest value.
 */
static bool secondComparator(std::pair< size_t, float >& a, std::pair< size_t, float >& b)
{
    return (a.second < b.second);
}

/**
 * @brief sortIndexRadiusMap
 * Function to sort the map according to value in a (key-value) pairs.
 * @param indexRadiusMap
 * Input map.
 */
static std::vector< std::pair< size_t , float > >
sortIndexRadiusMap(std::map< size_t, float >& indexRadiusMap)
{
    // Declare vector of pairs
    std::vector< std::pair< size_t , float > > pairsVector;

    // Copy key-value pair from Map to vector of pairs
    for (auto& it : indexRadiusMap)
    {
        pairsVector.push_back(std::make_pair(it.first, it.second));
    }

    // Sort using comparator function
    sort(pairsVector.begin(), pairsVector.end(), secondComparator);

    return pairsVector;
}

}
