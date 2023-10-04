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

#include "Skeletonizer.h"

namespace Ultraliser
{

void Skeletonizer::applyVolumeThinningToVolume(Volume* volume, const bool& displayProgress)
{
    std::unique_ptr< Thinning6Iterations > thinningKernel = std::make_unique<Thinning6Iterations>();

    if (displayProgress)
    {
        // Parameters to calculate the loop progress
        size_t initialNumberVoxelsToBeDeleted = 0;
        size_t loopCounter = 0;

        TIMER_SET;
        LOG_STATUS("Thinning Volume");
        LOOP_STARTS("Thinning Loop");
        LOOP_PROGRESS(0, 100);
        while(1)
        {
            size_t numberDeletedVoxels = volume->deleteCandidateVoxelsParallel(thinningKernel);

            // Updating the progess bar
           if (loopCounter == 0) initialNumberVoxelsToBeDeleted = numberDeletedVoxels;
           LOOP_PROGRESS(initialNumberVoxelsToBeDeleted - numberDeletedVoxels,
                         initialNumberVoxelsToBeDeleted);

           if (numberDeletedVoxels == 0)
               break;

           loopCounter++;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        // Parameters to calculate the loop progress
        size_t initialNumberVoxelsToBeDeleted = 0;
        size_t loopCounter = 0;
        while(1)
        {
            size_t numberDeletedVoxels = volume->deleteCandidateVoxelsParallel(thinningKernel);
           if (numberDeletedVoxels == 0)
               break;
        }
    }
}

void Skeletonizer::applyVolumeThinningWithDomainDecomposition()
{
    // Start the timer
    TIMER_SET;

    // Copy the input volume into a reference volume
    // The reference volume will be used to retrieve the new bricks, and the _volume will be
    // used to write the skeletonization result
    Volume* referenceVolume = new Volume(_volume->getWidth(),
                                         _volume->getHeight(),
                                         _volume->getDepth());

    const size_t subdivisions = 4;
    const size_t overlappingVoxels = 5;
    const size_t numberZeroVoxels = 2;

    Ranges xRanges = Range::decomposeToRanges(int64_t(0), _volume->getWidth() - 1, subdivisions);
    Ranges yRanges = Range::decomposeToRanges(int64_t(0), _volume->getHeight() - 1, subdivisions);
    Ranges zRanges = Range::decomposeToRanges(int64_t(0), _volume->getDepth() - 1, subdivisions);

    // Add the overlaps
    Ranges xRangesOverlapping = Range::addTwoSidedOverlaps(xRanges, overlappingVoxels);
    Ranges yRangesOverlapping = Range::addTwoSidedOverlaps(yRanges, overlappingVoxels);
    Ranges zRangesOverlapping = Range::addTwoSidedOverlaps(zRanges, overlappingVoxels);

    LOG_STATUS("Skeletonizing Volume Bricks");
    LOOP_STARTS("Skeletonization");
    int64_t progress = 0;
    for (size_t i = 0; i < xRanges.size(); ++i)
    {
        LOOP_PROGRESS(progress, xRanges.size());
        for (size_t j = 0; j < yRanges.size(); ++j)
        {
            // OMP_PARALLEL_FOR
            for (size_t k = 0; k < zRanges.size(); ++k)
            {
                // Extract the brick from the volume
                auto brick = _volume->extractBoundedBrickFromVolume(
                            xRangesOverlapping[i].i1, xRangesOverlapping[i].i2,
                            yRangesOverlapping[j].i1, yRangesOverlapping[j].i2,
                            zRangesOverlapping[k].i1, zRangesOverlapping[k].i2,
                            numberZeroVoxels, false);

                // Skeletonize the brick
                applyVolumeThinningToVolume(brick, false);

                size_t xOverlapping, yOverlapping, zOverlapping = 0;
                if (i > 0) xOverlapping = overlappingVoxels; else xOverlapping = 0;
                if (j > 0) yOverlapping = overlappingVoxels; else yOverlapping = 0;
                if (k > 0) zOverlapping = overlappingVoxels; else zOverlapping = 0;

                referenceVolume->insertOverlappingBoundedBrickToVolume(
                            brick,
                            xRanges[i].i1, xRanges[i].i2,
                            yRanges[j].i1, yRanges[j].i2,
                            zRanges[k].i1, zRanges[k].i2,
                            xOverlapping, yOverlapping, zOverlapping, numberZeroVoxels,
                            false);

                brick->~Volume();
            }
        }

        progress++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    _volume->insertBrickToVolume(referenceVolume,
                                 0, referenceVolume->getWidth() - 1,
                                 0, referenceVolume->getHeight() - 1,
                                 0, referenceVolume->getDepth() - 1);
    referenceVolume->~Volume();
}

void Skeletonizer::thinVolumeBlockByBlock(const size_t& blockSize,
                                          const size_t& numberOverlappingVoxels,
                                          const size_t& numberZeroVoxels)
{
    // Start the timer
    TIMER_SET;

    // Copy the input volume into a reference volume
    // The reference volume will be used to retrieve the new bricks, and the _volume will be
    // used to write the skeletonization result
    Volume* referenceVolume = new Volume(_volume->getWidth(),
                                         _volume->getHeight(),
                                         _volume->getDepth());

    // Initially, the range of the volume is identified
    Range xRange(0, _volume->getWidth() - 1);
    Range yRange(0, _volume->getHeight() - 1);
    Range zRange(0, _volume->getDepth() - 1);

    // Decompose the range into ranges, based on the blockSize
    Ranges xRanges = xRange.decomposeToBlocks(blockSize);
    Ranges yRanges = yRange.decomposeToBlocks(blockSize);
    Ranges zRanges = zRange.decomposeToBlocks(blockSize);

    for (const auto& range: xRanges)
        range.printRange();
    for (const auto& range: yRanges)
        range.printRange();
    for (const auto& range: zRanges)
        range.printRange();


    // Add the overlaps
    Ranges xRangesOverlapping = Range::addTwoSidedOverlaps(xRanges, numberOverlappingVoxels);
    Ranges yRangesOverlapping = Range::addTwoSidedOverlaps(yRanges, numberOverlappingVoxels);
    Ranges zRangesOverlapping = Range::addTwoSidedOverlaps(zRanges, numberOverlappingVoxels);

    LOG_STATUS("Skeletonizing Volume Bricks");
    LOOP_STARTS("Skeletonization");
    int64_t progress = 0;
    for (size_t i = 0; i < xRanges.size(); ++i)
    {
        LOOP_PROGRESS(progress, xRanges.size());
        for (size_t j = 0; j < yRanges.size(); ++j)
        {
            // OMP_PARALLEL_FOR
            for (size_t k = 0; k < zRanges.size(); ++k)
            {
                // Extract the brick from the volume
                auto brick = _volume->extractBoundedBrickFromVolume(
                            xRangesOverlapping[i].i1, xRangesOverlapping[i].i2,
                            yRangesOverlapping[j].i1, yRangesOverlapping[j].i2,
                            zRangesOverlapping[k].i1, zRangesOverlapping[k].i2,
                            numberZeroVoxels, false);

                // Skeletonize the brick
                applyVolumeThinningToVolume(brick, false);

                size_t xOverlapping, yOverlapping, zOverlapping = 0;
                if (i > 0) xOverlapping = numberZeroVoxels; else xOverlapping = 0;
                if (j > 0) yOverlapping = numberZeroVoxels; else yOverlapping = 0;
                if (k > 0) zOverlapping = numberZeroVoxels; else zOverlapping = 0;

                referenceVolume->insertOverlappingBoundedBrickToVolume(
                            brick,
                            xRanges[i].i1, xRanges[i].i2,
                            yRanges[j].i1, yRanges[j].i2,
                            zRanges[k].i1, zRanges[k].i2,
                            xOverlapping, yOverlapping, zOverlapping, numberZeroVoxels,
                            false);

                brick->~Volume();
            }
        }

        progress++;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    _volume->insertBrickToVolume(referenceVolume,
                                 0, referenceVolume->getWidth() - 1,
                                 0, referenceVolume->getHeight() - 1,
                                 0, referenceVolume->getDepth() - 1);
    referenceVolume->~Volume();
}

}
