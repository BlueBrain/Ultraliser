/***************************************************************************************************
 * Copyright (c) 2016 - 2022
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

#include <geometry/Geometry.h>
#include <data/morphologies/Sample.h>

namespace Ultraliser
{

/**
 * @brief The EndfootPatch class
 */
class EndfootPatch
{
public:

    /**
     * @brief EndfootPatch
     * Constructor
     * @param sample0
     * First vertex of the triangle, a point including a radius.
     * @param sample1
     * Second vertex of the triangle.
     * @param sample2
     * Third vertex of the triangle.
     */
    EndfootPatch(Sample* sample0, Sample* sample1, Sample* sample2)
      : sample0(sample0), sample1(sample1), sample2(sample2)
    {
        /// EMPTY CONSTRUCTOR
    };

public:

    /**
     * @brief sample0
     * First sample of the triangle.
     */
    Sample* sample0;

        /**
     * @brief sample1
     * Second sample of the triangle.
     */
    Sample* sample1;

    /**
     * @brief sample2
     * Third sample of the triangle.
     */
    Sample* sample2;
};

/**
 * @brief EndfeetPatches
 * A list of endfeet patches.
 */
typedef std::vector<EndfootPatch*> EndfootPatches;

/**
 * @brief EndfeetPatches
 */
typedef std::vector<EndfootPatches> EndfeetPatches;


}
