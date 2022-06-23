/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#pragma once

#include <vector>
#include <common/Headers.hh>

namespace Ultraliser
{

/**
 * @brief The VasculatureH5Connectivity struct
 * Connectivity between two sections.
 */
struct VasculatureH5Connectivity
{
    /**
     * @brief parentSectionIndex
     * The index of the parent section.
     */
    int64_t parentSectionIndex;

    /**
     * @brief childSectionIndex
     * The index of the child section.
     */
    int64_t childSectionIndex;

};

/**
 * @brief VasculatureH5ConnectivityList
 * A list of connectivity entries in the morphology.
 */
typedef std::vector< VasculatureH5Connectivity > VasculatureH5ConnectivityList;

}
