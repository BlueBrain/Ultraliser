/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include "VolumeGrid.h"
#include "Projection.h"

namespace Ultraliser
{

void VolumeGrid::composeProjections(const std::string &prefix,
                                    const bool &projectXY,
                                    const bool &projectYZ,
                                    const bool &projectZX,
                                    const bool &colorCodedProjections,
                                    const bool &verbose)
{
    if (projectXY || projectYZ || projectZX)
    {
        TIMER_SET;

        LOG_TITLE("Projecting Volume");
        LOG_STATUS("Compositing Projection(s)");


        if (projectXY)
            _composeProjectionXY(prefix, colorCodedProjections, verbose);

        if (projectYZ)
            _composeProjectionYZ(prefix, colorCodedProjections, verbose);

        if (projectZX)
            _composeProjectionXZ(prefix, colorCodedProjections, verbose);

        // Statistics
        LOG_STATUS_IMPORTANT("Volume Projection Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void VolumeGrid::_composeProjectionXY(const std::string &prefix,
                                      const bool &colorCodedProjections,
                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getWidth();
    const size_t& projectionHeight = getHeight();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("XY Projection - Z Axis");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - i - 1) + getWidth() * j] +=
                                getValueF64(i, j, k);
                    }
                }
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - i - 1) + getWidth() * j] +=
                                getValueF64(i, j, k);
                    }
                }
            }
        }
    }

    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + XY_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

void VolumeGrid::_composeProjectionYZ(const std::string &prefix,
                                      const bool &colorCodedProjections,
                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getHeight();
    const size_t& projectionHeight = getDepth();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("YZ Projection - X Axis");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - j - 1) + getHeight() * k] +=
                                getValueF64(i, j, k);
                    }
                }
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - j - 1) + getHeight() * k] +=
                                getValueF64(i, j, k);
                    }
                }
            }
        }
    }

    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + YZ_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

void VolumeGrid::_composeProjectionXZ(const std::string &prefix,
                                      const bool &colorCodedProjections,
                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getWidth();
    const size_t& projectionHeight = getDepth();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("XZ Projection - Y Axis");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - i - 1) + getWidth() * k] +=
                                getValueF64(i, j, k);
                    }
                }
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for (int64_t i = 0; i < getWidth(); i++)
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[(projectionWidth - i - 1) + getWidth() * k] +=
                                getValueF64(i, j, k);
                    }
                }
            }
        }
    }


    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + XZ_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

void VolumeGrid::composeProjectionsFromFilledVoxels(const std::string &prefix,
                                                    const bool &projectXY,
                                                    const bool &projectYZ,
                                                    const bool &projectZX,
                                                    const bool &colorCodedProjections,
                                                    const bool &verbose)

{
    if (projectXY || projectYZ || projectZX)
    {
        TIMER_SET;

        if (verbose) LOG_TITLE("Projecting Volume *");
        if (verbose) LOG_STATUS("Compositing Projection(s)");

        // Ensure that the filled voxels data structure is created
        auto filledVoxels = getFilledVoxels();

        if (projectXY)
            _composeProjectionFromFilledVoxelsXY(prefix, colorCodedProjections, verbose);

        if (projectYZ)
            _composeProjectionFromFilledVoxelsYZ(prefix, colorCodedProjections, verbose);

        if (projectZX)
            _composeProjectionFromFilledVoxelsXZ(prefix, colorCodedProjections, verbose);

        // Statistics
        if (verbose) LOG_STATUS_IMPORTANT("Volume Projection Stats.");
        if (verbose) LOG_STATS(GET_TIME_SECONDS);
    }
}

void VolumeGrid::_composeProjectionFromFilledVoxelsXY(const std::string &prefix,
                                                      const bool &colorCodedProjections,
                                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getWidth();
    const size_t& projectionHeight = getHeight();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    // Get a reference to the FilledVoxels structure
    auto filledVoxels = getFilledVoxels(verbose);

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("XY Projection - Z Axis *");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for(size_t sliceIndex = 0; sliceIndex < getDepth(); ++sliceIndex)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.z == sliceIndex)
                {
                    projectionImage[(projectionWidth - voxel.x - 1) + getWidth() * voxel.y] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }

            // Update the progress bar
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for(size_t sliceIndex = 0; sliceIndex < getDepth(); ++sliceIndex)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.z == sliceIndex)
                {
                    projectionImage[(projectionWidth - voxel.x - 1) + getWidth() * voxel.y] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }
        }
    }

    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + XY_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

void VolumeGrid::_composeProjectionFromFilledVoxelsYZ(const std::string &prefix,
                                                      const bool &colorCodedProjections,
                                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getHeight();
    const size_t& projectionHeight = getDepth();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    // Get a reference to the FilledVoxels structure
    auto filledVoxels = getFilledVoxels();

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("YZ Projection - X Axis *");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for(size_t sliceIndex = 0; sliceIndex < getWidth(); ++sliceIndex)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.x == sliceIndex)
                {
                    projectionImage[(projectionWidth - voxel.y - 1) + getHeight() * voxel.z] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for(size_t sliceIndex = 0; sliceIndex < getWidth(); ++sliceIndex)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.x == sliceIndex)
                {
                    projectionImage[(projectionWidth - voxel.y - 1) + getHeight() * voxel.z] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }
        }
    }

    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + YZ_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

void VolumeGrid::_composeProjectionFromFilledVoxelsXZ(const std::string &prefix,
                                                      const bool &colorCodedProjections,
                                                      const bool &verbose)
{
    // Dimensions
    const size_t& projectionWidth = getWidth();
    const size_t& projectionHeight = getDepth();
    const size_t& projectionSize = projectionWidth * projectionHeight;

    // Create a projection array (float)
    double* projectionImage = new double[projectionSize]();

    // Get a reference to the FilledVoxels structure
    auto filledVoxels = getFilledVoxels();

    if (verbose)
    {
        // Starts the timer
        TIMER_SET;

        LOOP_STARTS("XZ Projection - Y Axis *");
        PROGRESS_SET;
        OMP_PARALLEL_FOR
        for(size_t i = 0; i < getHeight(); ++i)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.y == i)
                {
                    projectionImage[(projectionWidth - voxel.x - 1) + getWidth() * voxel.z] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }
            LOOP_PROGRESS(PROGRESS, getWidth());
            PROGRESS_UPDATE;
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        OMP_PARALLEL_FOR
        for(size_t i = 0; i < getHeight(); ++i)
        {
            for (size_t n = 0; n < filledVoxels.size(); n++)
            {
                const auto& voxel = filledVoxels[n];
                if (voxel.y == i)
                {
                    projectionImage[(projectionWidth - voxel.x - 1) + getWidth() * voxel.z] +=
                            getValueF64(voxel.x, voxel.y, voxel.z);
                }
            }
        }
    }

    // Write the projections
    writeProjections(prefix + PROJECTION_SUFFIX + XZ_SUFFIX,
                     projectionImage, projectionWidth, projectionHeight,
                     colorCodedProjections, verbose);
}

}
