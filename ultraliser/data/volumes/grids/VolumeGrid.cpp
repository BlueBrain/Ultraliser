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
#include <algorithms/Algorithms.h>

namespace Ultraliser
{

VolumeGrid::VolumeGrid(const size_t &width, const size_t &height, const size_t &depth)
    : _width(width)
    , _height(height)
    , _depth(depth)
    , _numberVoxels(width * height * depth)
{
    /// EMPTY CONSTRUCTOR
}

VolumeGrid::VolumeGrid(Ultraliser::VolumeGrid *grid)
{
    // Copy the dimensions
    grid->getDimensions(_width, _height, _depth);

    // Compute the number of voxels
    _numberVoxels = _width * _height * _depth;

    /// TODO: Copy the data of the grid
}

void VolumeGrid::getDimensions(size_t& width, size_t& height, size_t& depth) const
{
    width = _width;
    height = _height;
    depth = _depth;
}

size_t VolumeGrid::mapTo1DIndexWOBC(const size_t &x, const size_t &y, const size_t &z) const
{
    return (x + (_width * y) + (_width * _height * z));
}

size_t VolumeGrid::mapTo1DIndex(const size_t &x, const size_t &y, const size_t &z,
                              bool &outlier) const
{
    if(x >= getWidth()  || x < 0 || y >= getHeight() || y < 0 || z >= getDepth()  || z < 0)
    {
        outlier = true;
        return 0;
    }
    else
    {
        outlier = false;
        return I2UI64(x + (_width * y) + (_width * _height * z));
    }
}

void VolumeGrid::fillVoxel(const int64_t &x, const int64_t &y, const int64_t &z)
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return;
    else
        fillVoxel(index);
}

void VolumeGrid::clearVoxel(const int64_t &x, const int64_t &y, const int64_t &z)
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return;
    else
        clearVoxel(index);
}

bool VolumeGrid::isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return false;
    else
        return isFilled(index);
}

bool VolumeGrid::isEmpty(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return true;
    else
        return isEmpty(index);
}

size_t VolumeGrid::getDimension(const size_t &i) const
{
    switch (i)
    {
        case 0 : return getWidth();
        case 1 : return getHeight();
        case 2 : return getDepth();
        default: return 0;
    }
}

size_t VolumeGrid::computeNumberNonZeroVoxelsPerSlice(int64_t z) const
{
    size_t numberNonZeroVoxels = 0;

    for (int64_t i = 0; i < getWidth(); ++i)
    {
        for (int64_t j = 0; j < getHeight(); ++j)
        {
            bool outlier;
            size_t index = mapTo1DIndex(i, j, z, outlier);
            bool filled = isFilled(index);
            if (filled && !outlier)
            {
                numberNonZeroVoxels += 1;
            }
        }
    }

    return numberNonZeroVoxels;
}

void VolumeGrid::floodFillSlice_X(const int64_t &sliceIndex, const size_t &padding)
{
    /// Select an YZ slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getHeight() + 2 * padding;
    const size_t sliceHeight = getDepth() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getHeight();
    const size_t sliceHeightWithinVolume = getDepth();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Fill the slice with the surface voxels (pixels in the image)
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(sliceIndex, i, j, outlier);
            if (isFilled(index) && !outlier)
                slice.setPixelColor(i + padding, j + padding, GRAY);
#else
            if (isFilled(mapTo1DIndexWOBC(sliceIndex, i, j)))
                slice.setPixelColor(i + padding, j + padding, GRAY);
#endif
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(sliceIndex, i, j, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            // Construct the 1D index of the voxel
            size_t index = mapTo1DIndexWOBC(sliceIndex, i, j);

            // If the slice is BLACK at this pixel, clear the corresponding voxel, otherwise set it
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(index);
            else
                fillVoxel(index);
#endif
        }
    }
}

void VolumeGrid::floodFillSliceUsingFilledVoxels_X(const int64_t &sliceIndex,
                                                   const size_t &padding)
{
    /// Select an YZ slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getHeight() + 2 * padding;
    const size_t sliceHeight = getDepth() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getHeight();
    const size_t sliceHeightWithinVolume = getDepth();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Use the filled voxels acceleration structure to label the shaded voxels from the previous
    // surface voxelization operation
    for (int64_t i = 0; i < _filledVoxels.size(); ++i)
    {
        if (_filledVoxels[i].x == sliceIndex)
        {
            slice.setPixelColor(_filledVoxels[i].y + padding, _filledVoxels[i].z + padding, GRAY);
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(sliceIndex, i, j, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            // If the slice is BLACK at this pixel, clear the corresponding voxel, otherwise set it
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(mapTo1DIndexWOBC(sliceIndex, i, j));
            else
                fillVoxel(mapTo1DIndexWOBC(sliceIndex, i, j));
#endif
        }
    }
}

void VolumeGrid::floodFillSlice_Y(const int64_t &sliceIndex, const size_t &padding)
{
    /// Select an XZ slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getWidth() + 2 * padding;
    const size_t sliceHeight = getDepth() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getWidth();
    const size_t sliceHeightWithinVolume = getDepth();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Fill the slice with the surface voxels (pixels in the image)
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            const size_t index = mapTo1DIndex(i, sliceIndex, j, outlier);
            if (isFilled(index) && !outlier)
                slice.setPixelColor(i + padding, j + padding, GRAY);
#else
            if (isFilled(mapTo1DIndexWOBC(i, sliceIndex, j)))
                slice.setPixelColor(i + padding, j + padding, GRAY);
#endif
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            const size_t index = mapTo1DIndex(i, sliceIndex, j, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            // If the slice is BLACK at this pixel, clear the corresponding voxel, otherwise set it
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(mapTo1DIndexWOBC(i, sliceIndex, j));
            else
                fillVoxel(mapTo1DIndexWOBC(i, sliceIndex, j));
#endif
        }
    }
}

void VolumeGrid::floodFillSliceUsingFilledVoxels_Y(const int64_t &sliceIndex, const size_t &padding)
{
    /// Select an XZ slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getWidth() + 2 * padding;
    const size_t sliceHeight = getDepth() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getWidth();
    const size_t sliceHeightWithinVolume = getDepth();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Use the filled voxels acceleration structure to label the shaded voxels from the previous
    // surface voxelization operation
    for (int64_t i = 0; i < _filledVoxels.size(); ++i)
    {
        if (_filledVoxels[i].y == sliceIndex)
        {
            slice.setPixelColor(_filledVoxels[i].x + padding, _filledVoxels[i].z + padding, GRAY);
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            const size_t index = mapTo1DIndex(i, sliceIndex, j, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            // If the slice is BLACK at this pixel, clear the corresponding voxel, otherwise set it
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(mapTo1DIndexWOBC(i, sliceIndex, j));
            else
                fillVoxel(mapTo1DIndexWOBC(i, sliceIndex, j));
#endif
        }
    }
}

void VolumeGrid::floodFillSlice_Z(const int64_t &sliceIndex,
                                  const size_t &padding)
{
    /// Select an XY slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getWidth() + 2 * padding;
    const size_t sliceHeight = getHeight() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getWidth();
    const size_t sliceHeightWithinVolume = getHeight();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Fill the slice with the surface voxels (pixels in the image)
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(i, j, sliceIndex, outlier);
            if (isFilled(index) && !outlier)
                slice.setPixelColor(i + padding, j + padding, GRAY);
#else
            if (isFilled(mapTo1DIndexWOBC(i, j, sliceIndex)))
                slice.setPixelColor(i + padding, j + padding, GRAY);
#endif
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(i, j, sliceIndex, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(mapTo1DIndexWOBC(i, j, sliceIndex));
            else
                fillVoxel(mapTo1DIndexWOBC(i, j, sliceIndex));
#endif
        }
    }
}

void VolumeGrid::floodFillSliceUsingFilledVoxels_Z(const int64_t &sliceIndex,
                                                   const size_t &padding)
{
    /// Select an XY slice from the volume, based on the given @sliceIndex and flood-fill it
    // Compute the dimensions of the slice after taking the zero-padding into consideration
    const size_t sliceWidth = getWidth() + 2 * padding;
    const size_t sliceHeight = getHeight() + 2 * padding;

    // The dimensions of the original slice must be preserved, without the zero-padded pixels
    const size_t sliceWidthWithinVolume = getWidth();
    const size_t sliceHeightWithinVolume = getHeight();

    // Create a slice and set its initial colors by default to WHITE
    Image slice(sliceWidth, sliceHeight, WHITE);

    // Use the filled voxels acceleration structure to label the shaded voxels from the previous
    // surface voxelization operation
    for (int64_t i = 0; i < _filledVoxels.size(); ++i)
    {
        if (_filledVoxels[i].x == sliceIndex)
        {
            slice.setPixelColor(_filledVoxels[i].x + padding, _filledVoxels[i].y + padding, GRAY);
        }
    }

    // Apply the flood-filling operation
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Copy the content of the slice back to the volume after the flood-filling operation
    for (int64_t i = 0; i < sliceWidthWithinVolume; ++i)
    {
        for (int64_t j = 0; j < sliceHeightWithinVolume; ++j)
        {
#ifdef ULTRALISER_DEBUG
            bool outlier;
            size_t index = mapTo1DIndex(i, j, sliceIndex, outlier);
            if (slice.getPixelColor(i , j) == BLACK && !outlier)
                clearVoxel(index);
            else
                fillVoxel(index);
#else
            if (slice.getPixelColor(i , j) == BLACK)
                clearVoxel(mapTo1DIndexWOBC(i, j, sliceIndex));
            else
                fillVoxel(mapTo1DIndexWOBC(i, j, sliceIndex));
#endif
        }
    }
}


void VolumeGrid::floodFillSliceAlongAxis(const int64_t &sliceIndex,
                                         const AXIS &axis,
                                         const size_t &padding)
{      
    // Slice dimensions
    int64_t sliceWidth, sliceHeight, sliceSize;

    // Volume dimensions along the slice
    int64_t volumeWidth, volumeHeight;

    switch (axis)
    {
    // YZ axis
    case AXIS::X:
    {
        sliceWidth = getHeight() + 2 * padding;
        sliceHeight = getDepth() + 2 * padding;
        sliceSize = sliceWidth * sliceHeight;

        volumeWidth = getHeight();
        volumeHeight = getDepth();
        break;
    }

        // XZ axis
    case AXIS::Y:
    {
        sliceWidth = getWidth() + 2 * padding;
        sliceHeight = getDepth() + 2 * padding;
        sliceSize = sliceWidth * sliceHeight;

        volumeWidth = getWidth();
        volumeHeight = getDepth();
        break;
    }

        // XY axis
    case AXIS::Z:
    {
        // Dimensions
        sliceWidth = getWidth() + 2 * padding;
        sliceHeight = getHeight() + 2 * padding;
        sliceSize = sliceWidth * sliceHeight;

        volumeWidth = getWidth();
        volumeHeight = getHeight();
        break;
    }
    }

    // Create an X-slice
    Image slice(sliceWidth, sliceHeight);

    // Make it blank
    slice.fill(WHITE);




    // Assuming that there is a point cloud


//    std::vector< Vec3ui_16 > a;
//    for (size_t n = 0; n < a.size(); ++n)
//    {
//        // X
//        if (a[n].x() == sliceIndex)
//        {
//            slice.setPixelColor(a[n].x() + padding, a[n].y() + padding, GRAY);
//        }
//    }


    switch (axis)
    {
    case AXIS::X:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(sliceIndex, i, j, outlier);
                if (isFilled(index) && !outlier)
                    slice.setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice.setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    case AXIS::Y:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i, sliceIndex, j, outlier);
                if (isFilled(index) && !outlier)
                    slice.setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice.setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    case  AXIS::Z:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i, j, sliceIndex, outlier);
                if (isFilled(index) && !outlier)
                    slice.setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice.setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    }

    // Flood Filler
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(&slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Update the volume back
    switch (axis)
    {
    case AXIS::X:
    {
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(sliceIndex, i, j, outlier);
                if (slice.getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;

    case AXIS::Y:
    {
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i, sliceIndex, j, outlier);
                if (slice.getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;

    case  AXIS::Z:
    {
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i, j, sliceIndex, outlier);
                if (slice.getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;

    }
}


void VolumeGrid::floodFillSliceAlongAxisROI(const int64_t &sliceIndex,
                                            const AXIS &axis,
                                            const size_t& x1, const size_t x2,
                                            const size_t& y1, const size_t y2,
                                            const size_t& z1, const size_t z2,
                                            const size_t &padding)
{
    // Slice dimensions
    int64_t sliceWidth, sliceHeight, sliceSize;

    // The dimensions of the ROI within the volume itself
    size_t widthROI, heightROI;

    switch (axis)
    {
    // YZ axis
    case AXIS::X:
    {
        // Set the dimensions of the ROI (slice) within the volume itself
        widthROI = y2 - y1 + 1;
        heightROI = z2 - z1 + 1;
    } break;

    // XZ axis
    case AXIS::Y:
    {
        // Set the dimensions of the ROI (slice) within the volume itself
        widthROI = (x2 - x1 + 1);
        heightROI = (z2 - z1 + 1);
    } break;

    // XY axis
    case AXIS::Z:
    {
        // Set the dimensions of the ROI (slice) within the volume itself

        widthROI = (x2 - x1 + 1);
        heightROI = (y2 - y1 + 1);
    } break;
    }

    // Set the dimensions of the slice that will be used for the flood-filling, add the padding
    sliceWidth = widthROI + (2 * padding);
    sliceHeight = heightROI + (2 * padding);
    sliceSize = sliceWidth * sliceHeight;

    // Create an X-slice
    Image* slice = new Image(sliceWidth, sliceHeight);

    // Make it blank
    slice->fill(WHITE);

    switch (axis)
    {
    case AXIS::X:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(sliceIndex, i + y1, j + z1, outlier);
                if (isFilled(index) && !outlier)
                    slice->setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice->setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    case AXIS::Y:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i + x1, sliceIndex, j + z1, outlier);
                if (isFilled(index) && !outlier)
                    slice->setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice->setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    case  AXIS::Z:
    {
        // Fill the slice with the surface voxels (pixels in the image)
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i + x1, j + y1, sliceIndex, outlier);
                if (isFilled(index) && !outlier)
                    slice->setPixelColor(i + padding, j + padding, GRAY);
                else
                    slice->setPixelColor(i + padding, j + padding, WHITE);
            }
        }
    } break;

    }

    // Flood Filler
    PIXEL_COLOR newColor = WHITE;
    PIXEL_COLOR oldColor = BLACK;
    FloodFiller::fill(slice, sliceWidth, sliceHeight, 0, 0, newColor, oldColor);

    // Update the volume back
    switch (axis)
    {
    case AXIS::X:
    {
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(sliceIndex, i + y1, j + z1, outlier);
                if (slice->getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;

    case AXIS::Y:
    {
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i + x1, sliceIndex, j + z1, outlier);
                if (slice->getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;

    case  AXIS::Z:
    {
        for (int64_t i = 0; i < widthROI; ++i)
        {
            for (int64_t j = 0; j < heightROI; ++j)
            {
                bool outlier;
                size_t index = mapTo1DIndex(i + x1, j + y1, sliceIndex, outlier);
                if (slice->getPixelColor(i , j) == BLACK && !outlier)
                    clearVoxel(index);
                else
                    fillVoxel(index);
            }
        }
    } break;
    }

    delete slice;
}

void VolumeGrid::projectVolume(const std::string &prefix,
                               const bool &xyProjection,
                               const bool &xzProjection,
                               const bool &zyProjection,
                               const bool &colorCodedProjection)
{
    size_t numberProjections = 0;
    if (xyProjection ) numberProjections++;
    if (xzProjection) numberProjections++;
    if (zyProjection) numberProjections++;

    if (numberProjections > 1)
    {
        composeProjectionsFromFilledVoxels(prefix, xyProjection, zyProjection, xzProjection, colorCodedProjection, true);
    }
    else
    {
        composeProjections(prefix, xyProjection, zyProjection, xzProjection, colorCodedProjection, true);
    }

}

size_t VolumeGrid::computeNumberNonZeroVoxels() const
{
    // Starts the timer
    TIMER_SET;

    // The total number of of non zero voxels that will be computed for the entire slice
    size_t numberNonZeroVoxels = 0;

    // Number of non-zero voxels for every slice along the volume
    auto numberNonZeroVoxelsPerSlice = std::make_unique< uint64_t []>(I2UI64(getDepth()));

    // Initialize the arrays to zero
    OMP_PARALLEL_FOR
    for (int64_t k = 0; k < getDepth(); k++)
        numberNonZeroVoxelsPerSlice[k] = 0;

    LOOP_STARTS("Computing Filled Voxels")
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0; i < getDepth(); ++i)
    {
        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, getDepth());
        PROGRESS_UPDATE;

        numberNonZeroVoxelsPerSlice[i] = computeNumberNonZeroVoxelsPerSlice(i);
    }

    for (int64_t k = 0; k < getDepth(); k++)
    {
        numberNonZeroVoxels += numberNonZeroVoxelsPerSlice[k];
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Return the total number of non zero voxels in the entire volume
    return numberNonZeroVoxels;
}

void VolumeGrid::_buildVolumeOccupancy()
{
    TIMER_SET;
    _volumeOccupancy.resize(getWidth());

    PROGRESS_SET;
    LOOP_STARTS("Building Volume Occupancy Acceleration Structure");
    OMP_PARALLEL_FOR
    for (size_t i = 0; i < getWidth(); ++i)
    {
        auto& sliceOccupancy = _volumeOccupancy[i];
        sliceOccupancy.resize(getHeight());

        for (size_t j = 0; j < getHeight(); ++j)
        {
            for (size_t k = 0; k < getDepth(); ++k)
            {
                // If this voxel if solid, i.e. is filled
                if (isFilled(i, j, k))
                {
                    // Add the range to the list
                    sliceOccupancy[j].push_back(k);
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

VolumeOccpuancy VolumeGrid::getVolumeOccupancy()
{
    if (_volumeOccupancy.size() == 0)
        _buildVolumeOccupancy();

    return _volumeOccupancy;
}

void VolumeGrid::_queryFilledVoxels(const bool& verbose)
{
    // Start the timer
    TIMER_SET;

    // Process the volume in parallel
    // TODO: Make the query based on the dimension with the largest side
    std::vector< VoxelsXYZUI16 > filledVoxelsPerSlice;
    filledVoxelsPerSlice.resize(getWidth());

    if (verbose)
    {
        PROGRESS_SET;
        LOOP_STARTS("Querying Filled Voxels");
        OMP_PARALLEL_FOR
                for (size_t i = 0; i < getWidth(); ++i)
        {
            for(size_t j = 0; j < getHeight(); ++j)
            {
                for (size_t k = 0; k < getDepth(); ++k)
                {
                    if (isFilled(i, j, k))
                    {
                        filledVoxelsPerSlice[i].push_back(VoxelXYZUI16(i, j, k));
                    }
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
        for (size_t i = 0; i < getWidth(); ++i)
        {
            for(size_t j = 0; j < getHeight(); ++j)
            {
                for (size_t k = 0; k < getDepth(); ++k)
                {
                    if (isFilled(i, j, k))
                    {
                        filledVoxelsPerSlice[i].push_back(VoxelXYZUI16(i, j, k));
                    }
                }
            }
        }
        LOOP_DONE;
    }

    // Clear the list
    _filledVoxels.clear();
    _filledVoxels.shrink_to_fit();
    for (size_t i = 0; i < filledVoxelsPerSlice.size(); ++i)
    {
        _filledVoxels.insert(_filledVoxels.end(),
                             filledVoxelsPerSlice[i].begin(), filledVoxelsPerSlice[i].end());

        filledVoxelsPerSlice[i].clear();
        filledVoxelsPerSlice[i].shrink_to_fit();
    }
}

VoxelsXYZUI16 VolumeGrid::getFilledVoxels(const bool& verbose)
{
    if (_filledVoxels.size() == 0)
        _queryFilledVoxels(verbose);
    return _filledVoxels;
}

VOLUME_TYPE VolumeGrid::getType(const std::string &typeString)
{
    if (typeString == "bit")
    {
        return VOLUME_TYPE::BIT;
    }
    else if (typeString == "byte")
    {
        return VOLUME_TYPE::UI8;
    }
    else
    {
        LOG_WARNING("The volume type [ %s ] is not correct, using [bit]");
        return VOLUME_TYPE::BIT;
    }
}

VOLUME_TYPE VolumeGrid::getVolumeTypeFromHdrFile(const std::string& filePrefix)
{
    // Get the header file path from its prefix
    const std::string filePath = filePrefix + HEADER_EXTENSION;

    // Open the header file
    std::ifstream hdrFileStream(filePath.c_str());

    // Read the type
    std::string type;
    hdrFileStream >> type;

    // Close the stream
    hdrFileStream.close();

    if (type == FORMAT_BIT)
    {
        return VOLUME_TYPE::BIT;
    }
    else if (type == FORMAT_8UI)
    {
        return VOLUME_TYPE::UI8;
    }
    else if (type == FORMAT_16UI)
    {
        return VOLUME_TYPE::UI16;
    }
    else if (type == FORMAT_32UI)
    {
        return VOLUME_TYPE::UI32;
    }
    else if (type == FORMAT_64UI)
    {
        return VOLUME_TYPE::UI64;
    }
    else if (type == FORMAT_F32)
    {
        return VOLUME_TYPE::F32;
    }
    else if (type == FORMAT_F64)
    {
        return VOLUME_TYPE::F64;
    }
    else
    {
        LOG_ERROR("Volume type is NOT defined!");
        return VOLUME_TYPE::BIT;
    }
}

std::string VolumeGrid::getTypeString(const VOLUME_TYPE& type)
{
    if (type == VOLUME_TYPE::BIT)
    {
        return std::string("Bit");
    }
    else if (type == VOLUME_TYPE::UI8)
    {
        return std::string("Byte");
    }
    else
    {
        LOG_ERROR("The volume type is not correct");
        return std::string("");
    }
}

uint8_t VolumeGrid::getValueUI8(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueUI8(index);
}

uint16_t VolumeGrid::getValueUI16(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueUI16(index);
}

uint32_t VolumeGrid::getValueUI32(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueUI32(index);
}

uint64_t VolumeGrid::getValueUI64(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueUI64(index);
}

float VolumeGrid::getValueF32(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueF32(index);
}

double VolumeGrid::getValueF64(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    size_t index = mapTo1DIndex(x, y, z, outlier);
    if (outlier)
        return 0;
    else
        return getValueF64(index);
}

VolumeGrid::~VolumeGrid()
{
    /// EMPTY DESTRUCTOR
}

}
