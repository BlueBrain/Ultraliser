/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
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
 **************************************************************************************************/

#pragma once

#include <common/Headers.hh>

namespace Ultraliser
{

/**
 * @brief The OccpuanyRange class
 */
struct OccpuancyRange
{
public:

    /**
     * @brief OccpuanyRange
     * @param lower
     * @param upper
     */
    OccpuancyRange(size_t lower, size_t upper)
    {
        this->lower = lower;
        this->upper = upper;
    }

    /**
     * @brief lower
     */
    size_t lower;

    /**
     * @brief upper
     */
    size_t upper;
};

/**
 * @brief OccpuancyRanges
 * A list of OccpuanyRange's along a single dimension, for example Y line
 */
typedef std::vector< OccpuancyRange > OccpuancyRanges;

/**
 * @brief SliceOccpuancyRanges
 * A list of OccpuanyRange's for a single slice (plane), for example ZY slice
 */
typedef std::vector< OccpuancyRanges > SliceOccpuancyRanges;

/**
 * @brief VolumeOccpuancyRanges
 * A list of OccpuanyRange's for the entire volume
 */
typedef std::vector < SliceOccpuancyRanges > VolumeOccpuancyRanges;

}

