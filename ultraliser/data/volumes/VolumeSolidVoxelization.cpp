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

#include "Volume.h"
#include <common/Common.h>
#include <data/volumes/utilities/VolumeType.hh>
#include <geometry/Intersection.h>
#include <geometry/Utilities.h>
#include <utilities/Utilities.h>
#include <math/Functions.h>
#include <data/images/TIFFImage.h>
#include <data/volumes/grids/VolumeGrid.h>
#include <data/volumes/voxels/DMCVoxel.h>
#include <data/volumes/grids/BitVolumeGrid.h>
#include <data/volumes/grids/UnsignedVolumeGrid.h>
#include <data/volumes/grids/Grids.h>
#include <data/meshes/simple/VolumeMesh.h>
#include <data/meshes/simple/MeshOperations.h>
#include <data/volumes/utilities/VolumeReader.h>
#include <data/meshes/simple/TriangleOperations.h>

namespace Ultraliser
{

void Volume::solidVoxelization(const SOLID_VOXELIZATION_AXIS& axis,
                               const bool& useAcceleratedStructure,
                               const bool& verbose)
{
    if (verbose) LOG_TITLE("Solid Voxelization");

    // The 2D flood filling is only supported for the solid voxelization
    if (verbose) LOG_STATUS("Flood-filling Volume");
    _floodFill2D(axis, useAcceleratedStructure, verbose);

    if (verbose) LOG_STATUS_IMPORTANT("Solid Voxelization Stats.");
    LOG_STATS(_solidVoxelizationTime);
}

void Volume::_floodFillX(VolumeGrid* grid,
                         const bool& useAcceleratedStructure,
                         const bool& verbose)
{
    // Maximum number of slices
    const int64_t& numberSlices = getWidth();

    if (verbose)
    {
        // Start the timer
        TIMER_SET;

        // Disable buffering
        setbuf(stdout, nullptr);

        PROGRESS_SET;
        if (useAcceleratedStructure)
        {
            LOOP_STARTS("Accelerated Slice Flood-filling (X-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_X(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        else
        {
            LOOP_STARTS("Slice Flood-filling (X-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_X(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        if (useAcceleratedStructure)
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_X(i);
            }
        }
        else
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_X(i);
            }
        }
    }
}

void Volume::_floodFillY(VolumeGrid* grid, const bool& useAcceleratedStructure, const bool &verbose)
{
    // Maximum number of slices
    const int64_t& numberSlices = getHeight();

    if (verbose)
    {
        // Start the timer
        TIMER_SET;

        // Disable buffering
        setbuf(stdout, nullptr);

        PROGRESS_SET;
        if (useAcceleratedStructure)
        {
            LOOP_STARTS("Accelerated Slice Flood-filling (Y-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_Y(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        else
        {
            LOOP_STARTS("Slice Flood-filling (X-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_Y(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        if (useAcceleratedStructure)
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_Y(i);
            }
        }
        else
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_Y(i);
            }
        }
    }
}

void Volume::_floodFillZ(VolumeGrid* grid, const bool& useAcceleratedStructure, const bool &verbose)
{
    // Maximum number of slices
    const int64_t& numberSlices = getDepth();

    if (verbose)
    {
        // Start the timer
        TIMER_SET;

        // Disable buffering
        setbuf(stdout, nullptr);

        PROGRESS_SET;
        if (useAcceleratedStructure)
        {
            LOOP_STARTS("Accelerated Slice Flood-filling (Z-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_Z(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        else
        {
            LOOP_STARTS("Slice Flood-filling (X-axis)");
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_Z(i);

                // Update the progress bar
                LOOP_PROGRESS(PROGRESS, numberSlices);
                PROGRESS_UPDATE;
            }
            LOOP_DONE;
        }
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        if (useAcceleratedStructure)
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSliceUsingFilledVoxels_Z(i);
            }
        }
        else
        {
            OMP_PARALLEL_FOR
            for (int64_t i = 0 ; i < numberSlices; ++i)
            {
                grid->floodFillSlice_Z(i);
            }
        }

    }
}


}
