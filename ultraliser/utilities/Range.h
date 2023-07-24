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

#include <common/Common.h>

namespace Ultraliser
{

// Forward declaration
class Range;

/**
 * @brief Ranges
 */
typedef std::vector< Range > Ranges;

/**
 * @brief The Range class
 */
class Range
{
public:

    /**
     * @brief Range
     * Constructor
     * @param minValue
     * The lowe limit or the minimum value of the range.
     * @param maxValue
     * The upper limit or the maximum value of the range.
     */
    Range(const int64_t& minValue, const int64_t& maxValue);

public:

    /**
     * @brief i1
     * The lower bound.
     */
    int64_t i1;

    /**
     * @brief i2
     * The upper bound.
     */
    int64_t i2;

public:

    /**
     * @brief printRange
     * Print the bounds of the range.
     */
    void printRange() const;

    /**
     * @brief decomposeToRanges
     * Decomposes the range into a list of ranges based on the given intervals.
     * @param intervals
     * Number of intervals with which the range will be decomposed.
     * @return
     * A list of range of type std::vector< Range >.
     */
    Ranges decomposeToRanges(const size_t& intervals) const;

    /**
     * @brief decomposeToBlocks
     * @param blockSize
     * @return
     */
    Ranges decomposeToBlocks(const size_t& blockSize);

public:

    /**
     * @brief decomposeToRanges
     * Decomposes the range into a list of ranges based on the given intervals.
     * @param minValue
     * The minimum value of the range.
     * @param maxValue
     * The maximum value of the range.
     * @param intervals
     * Number of intervals with which the range will be decomposed.
     * @return
     * A list of range of type std::vector< Range >.
     */
    static Ranges decomposeToRanges(const int64_t& minValue,
                                    const int64_t& maxValue,
                                    const size_t& intervals);

    /**
     * @brief addTwoSidedOverlaps
     * @param inputRanges
     * @param overlap
     * @return
     */
    static Ranges addTwoSidedOverlaps(const Ranges& inputRanges, const size_t& overlap);
};

}
