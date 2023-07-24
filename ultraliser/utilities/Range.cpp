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

#include "Range.h"

namespace Ultraliser
{

Range::Range(const int64_t& minValue, const int64_t& maxValue)
{
    i1 = minValue < maxValue ? minValue : maxValue;
    i2 = maxValue > minValue ? maxValue : minValue;
}

void Range::printRange() const
{
    std::cout << i1 << ", " << i2 << ", delta = " << i2 - i1 << std::endl;
}

Ranges Range::decomposeToRanges(const size_t& intervals) const
{
    // Returned list
    Ranges ranges;

    // In case the two values are the same, only a single range is returned
    if (i1 == i2)
    {
        ranges.push_back(Range(i1, i2));
        return ranges;
    }

    // Compute the delta value
    const size_t delta = std::ceil((i2 - i1) / intervals);

    // Compute the ranges
    int64_t r1, r2;
    for (size_t i = 0; i < intervals; ++i)
    {
        // Compute the lower limit
        if (i == 0) r1 = 0; else r1 = r2 + 1;

        // Compute the upper limit
        r2 = r1 + delta + 1;
        if (r2 > i2) r2 = i2;

        ranges.push_back(Range(r1, r2));
    }

    // Return the ranges
    return ranges;
}

Ranges Range::decomposeToBlocks(const size_t& blockSize)
{
    // Returned list
    Ranges ranges;

    // Ensure that the blockSize is smaller than the current range
    const size_t currentRange = i2 - i1;
    if (blockSize > currentRange)
    {
        ranges.push_back(Range(i1, i2));
        return ranges;
    }
    else
    {
        // Determine the number of blocks
        const size_t numberBlocks = static_cast< size_t >(
                    std::ceil(float(currentRange) / float(blockSize)));

        int64_t r1, r2;
        for (size_t i = 0; i < numberBlocks; ++i)
        {
            // Compute the lower limit
            if (i == 0) r1 = i1; else r1 = r2 + 1;

            // Compute the upper limit
            r2 = r1 + blockSize + 1;
            if (r2 > i2) r2 = i2;

            ranges.push_back(Range(r1, r2));
        }
    }

    return ranges;
}

Ranges Range::decomposeToRanges(const int64_t& minValue,
                                const int64_t& maxValue,
                                const size_t& intervals)
{
    // Returned list
    Ranges ranges;

    // In case the two values are the same, only a single range is returned
    if (minValue == maxValue)
    {
        ranges.push_back(Range(minValue, maxValue));
        return ranges;
    }

    // Verify the limits
    const int64_t i1 = minValue < maxValue ? minValue : maxValue;
    const int64_t i2 = maxValue > minValue ? maxValue : minValue;

    // Compute the delta value
    const size_t delta = std::ceil((i2 - i1) / intervals);

    // Compute the ranges
    int64_t r1, r2;
    for (size_t i = 0; i < intervals; ++i)
    {
        // Compute the lower limit
        if (i == 0) r1 = i1; else r1 = r2 + 1;

        // Compute the upper limit
        r2 = r1 + delta + 1;
        if (r2 > i2) r2 = i2;

        ranges.push_back(Range(r1, r2));
    }

    // Return the ranges
    return ranges;
}

Ranges Range::addTwoSidedOverlaps(const Ranges& inputRanges, const size_t& overlap)
{
    // Another list for the output ranges
    Ranges outputRanges;

    // First range
    outputRanges.push_back(Range(inputRanges.front().i1,
                                 inputRanges.front().i2 + overlap));

    // Intermediate ranges
    if (inputRanges.size() > 2)
    {
        for (size_t i = 1; i < inputRanges.size() - 1; ++i)
        {
            outputRanges.push_back(Range(inputRanges.at(i).i1 - overlap,
                                         inputRanges.at(i).i2 + overlap));
        }
    }

    // Last range
    outputRanges.push_back(Range(inputRanges.back().i1 - overlap,
                                 inputRanges.back().i2));

    // Return the result
    return outputRanges;
}


}
