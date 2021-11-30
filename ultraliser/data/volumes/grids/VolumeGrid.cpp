/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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
#include <algorithms/Algorithms.h>

namespace Ultraliser
{

VolumeGrid::VolumeGrid(const Vec3i_64& dimensions)
{
    // Update dimensions
    _dimensions = dimensions;

    // Update number of voxels
    _numberVoxels = _dimensions.v[0] * _dimensions.v[1] * _dimensions.v[2];
}

VolumeGrid::VolumeGrid(const int64_t &width,
                       const int64_t &height,
                       const int64_t &depth)
{
    // Update dimensions
    _dimensions.v[0] = width;
    _dimensions.v[1] = height;
    _dimensions.v[2] = depth;

    // Update number of voxels
    _numberVoxels = width * height * depth;
}

VolumeGrid::VolumeGrid(VolumeGrid *grid)
{
    // Update dimensions
    _dimensions = grid->getDimensions();

    // Update number of voxels
    _numberVoxels = _dimensions.v[0] * _dimensions.v[1] * _dimensions.v[2];
}

Vec3i_64 VolumeGrid::getDimensions() const
{
    return _dimensions;
}

int64_t VolumeGrid::getWidth() const
{
    return _dimensions.v[0];
}

int64_t VolumeGrid::getHeight() const
{
    return _dimensions.v[1];
}

int64_t VolumeGrid::getDepth() const
{
    return _dimensions.v[2];
}

int64_t VolumeGrid::getNumberVoxels() const
{
    return _numberVoxels;
}

uint64_t VolumeGrid::mapToIndex(const int64_t &x, const int64_t &y, const int64_t &z, bool &outlier) const
{
    if(x >= getWidth()  || x < 0 || y >= getHeight() || y < 0 || z >= getDepth()  || z < 0)
    {
        outlier = true;
        return 0;
    }
    else
    {
        outlier = false;
        return I2UI64(x + (_dimensions.v[0] * y) + (_dimensions.v[0] * _dimensions.v[1] * z));
    }
}

void VolumeGrid::fillVoxel(const int64_t &x, const int64_t &y, const int64_t &z)
{
    bool outlier;
    uint64_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return;
    else
        fillVoxel(index);
}

void VolumeGrid::clearVoxel(const int64_t &x, const int64_t &y, const int64_t &z)
{
    bool outlier;
    uint64_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return;
    else
        clearVoxel(index);
}

bool VolumeGrid::isFilled(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    uint64_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return false;
    else
        return isFilled(index);
}

bool VolumeGrid::isEmpty(const int64_t &x, const int64_t &y, const int64_t &z) const
{
    bool outlier;
    uint64_t index = mapToIndex(x, y, z, outlier);
    if (outlier)
        return true;
    else
        return isEmpty(index);
}

int64_t VolumeGrid::getDimension(const int32_t &i) const
{
    switch (i)
    {
        case 0 : return getWidth();
        case 1 : return getHeight();
        case 2 : return getDepth();
        default: return 0;
    }
}

uint64_t VolumeGrid::computeNumberNonZeroVoxelsPerSlice(int64_t z) const
{
    uint64_t numberNonZeroVoxels = 0;

    for (int64_t i = 0; i < getWidth(); ++i)
    {
        for (int64_t j = 0; j < getHeight(); ++j)
        {
            bool outlier;
            uint64_t index = mapToIndex(i, j, z, outlier);
            bool filled = isFilled(index);
            if (filled && !outlier)
            {
                numberNonZeroVoxels += 1;
            }
        }
    }

    return numberNonZeroVoxels;
}


void VolumeGrid::floodFillSliceAlongAxis(const int64_t &sliceIndex,
                                         const AXIS &axis,
                                         const uint64_t &padding)
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
    Image* slice = new Image(sliceWidth, sliceHeight);

    // Make it blank
    slice->fill(WHITE);

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
                uint64_t index = mapToIndex(sliceIndex, i, j, outlier);
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
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                uint64_t index = mapToIndex(i, sliceIndex, j, outlier);
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
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                uint64_t index = mapToIndex(i, j, sliceIndex, outlier);
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
        for (int64_t i = 0; i < volumeWidth; ++i)
        {
            for (int64_t j = 0; j < volumeHeight; ++j)
            {
                bool outlier;
                uint64_t index = mapToIndex(sliceIndex, i, j, outlier);
                if (slice->getPixelColor(i , j) == BLACK && !outlier)
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
                uint64_t index = mapToIndex(i, sliceIndex, j, outlier);
                if (slice->getPixelColor(i , j) == BLACK && !outlier)
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
                uint64_t index = mapToIndex(i, j, sliceIndex, outlier);
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

void VolumeGrid::writeProjection(const std::string &prefix,
                                 const PROJECTION &projection,
                                 const bool &projectColorCoded)
{
    // Starts the timer
    TIMER_SET;

    // Projection dimensions
    int64_t projectionWidth, projectionHeight, projectionSize;

    // Projection prefix
    std::stringstream prefixStream;

    // Projection string
    std::string projectionString;

    switch (projection)
    {
    case PROJECTION::XY_PROJECTION:
    {
        // Dimensions
        projectionWidth = getWidth();
        projectionHeight = getHeight();
        projectionSize = projectionWidth * projectionHeight;

        // Prefix
        prefixStream << prefix << PROJECTION_SUFFIX << XY_SUFFIX;

        // String
        projectionString = "XY Projection";
        break;
    }
    case PROJECTION::XZ_PROJECTION:
    {
        projectionWidth = getWidth();
        projectionHeight = getDepth();
        projectionSize = projectionWidth * projectionHeight;


        // Prefix
        prefixStream << prefix << PROJECTION_SUFFIX << XZ_SUFFIX;

        // String
        projectionString = "XZ Projection";
        break;
    }
    case PROJECTION::ZY_PROJECTION:
    {
        projectionWidth = getDepth();
        projectionHeight = getHeight();
        projectionSize = projectionWidth * projectionHeight;

        // Prefix
        prefixStream << prefix << PROJECTION_SUFFIX << YZ_SUFFIX;

        // String
        projectionString = "ZY Projection";
        break;
    }
    }

    // Create a projection array (float)
    std::vector< float > projectionImage(projectionSize);

    // Create normalized projection array (0 - 255)
    std::vector< uint8_t > normalizedProjectionImage(projectionSize);

    // Initialize the projections to zero to avoid garbage
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t index = 0; index < projectionSize; ++index)
    {
        projectionImage[index] = 0.f;
        normalizedProjectionImage[index] = 0;
    }

    LOOP_STARTS(projectionString.c_str());
    PROGRESS_RESET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0; i < getWidth(); i++)
    {
        switch (projection)
        {
        case PROJECTION::XY_PROJECTION:
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[i + getWidth() * j] += 1.0;
                    }
                }
            }
        } break;

        case PROJECTION::XZ_PROJECTION:
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[i + getWidth() * k] += 1.0;
                    }
                }
            }
        } break;

        case PROJECTION::ZY_PROJECTION:
        {
            for (int64_t j = 0; j < getHeight(); j++)
            {
                for (int64_t k = 0; k < getDepth(); k++)
                {
                    if (isFilled(i, j, k))
                    {
                        projectionImage[k + getDepth() * j] += 1.0;
                    }
                }
            }
        } break;

        }

        // Update the progress bar
        LOOP_PROGRESS(PROGRESS, getWidth());
        PROGRESS_UPDATE;
    }
    LOOP_DONE;

    // Get the maximum value
    float maxValue = 0.f;
    for (int64_t index = 0; index < projectionSize; ++index)
    {
        if (projectionImage[index] > maxValue)
            maxValue = projectionImage[index];
    }

    // Construct the normalized projection
    OMP_PARALLEL_FOR
    for (int64_t index = 0; index < projectionSize; ++index)
    {
        // Compute float pixel value
        float pixelValue = float(255.0f) * projectionImage[index] / float(maxValue);

        // Convert to uint8_t to be able to write it to the image
        normalizedProjectionImage[index] = F2UI8(pixelValue);
    }

    // Save the projection into a PPM image
    Utilities::savePPMLuminanceImage(prefixStream.str(), normalizedProjectionImage.data(),
                                     projectionWidth, projectionHeight);

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Save color coded projections with all possible color-maps
    if (projectColorCoded)
    {
        saveColorMappedProjectionWithAllColorMaps(prefixStream.str(), projectionImage.data(),
                                                  projectionWidth, projectionHeight, 0, maxValue);
    }
}
void VolumeGrid::projectVolume(const std::string &prefix,
                               const bool &xyProjection,
                               const bool &xzProjection,
                               const bool &zyProjection,
                               const bool &colorCodedProjection)
{
    if (xyProjection || xzProjection || zyProjection)
    {
        TIMER_SET;

        LOG_TITLE("Projecting Volume");
        LOG_STATUS("Compositing Projection(s)");

        if (xyProjection)
        {
            writeProjection(prefix, PROJECTION::XY_PROJECTION, colorCodedProjection);
        }

        if (xzProjection)
        {
            writeProjection(prefix, PROJECTION::XZ_PROJECTION, colorCodedProjection);
        }

        if (zyProjection)
        {
            writeProjection(prefix, PROJECTION::ZY_PROJECTION, colorCodedProjection);
        }

        // Statistics
        _projectionTime = GET_TIME_SECONDS;
        LOG_STATUS_IMPORTANT("Volume Projection Stats.");
        LOG_STATS(_projectionTime);
    }
}

uint64_t VolumeGrid::computeNumberNonZeroVoxels() const
{
    // Starts the timer
    TIMER_SET;

    // The total number of of non zero voxels that will be computed for the entire slice
    uint64_t numberNonZeroVoxels = 0;

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

VolumeGrid::~VolumeGrid()
{
    /// EMPTY
}

VolumeGrid::TYPE VolumeGrid::getType(const std::string &typeString)
{
    if (typeString == "bit")
    {
        return TYPE::BIT;
    }
    else if (typeString == "byte")
    {
        return TYPE::BYTE;
    }
    else if (typeString == "voxel")
    {
        return TYPE::VOXEL;
    }
    else
    {
        LOG_WARNING("The volume type [ %s ] is not correct, using [bit]");
        return TYPE::BIT;
    }
}

std::string VolumeGrid::getTypeString(const VolumeGrid::TYPE& type)
{
    if (type == TYPE::BIT)
    {
        return std::string("Bit");
    }
    else if (type == TYPE::BYTE)
    {
        return std::string("Byte");
    }
    else if (type == TYPE::VOXEL)
    {
        return std::string("Voxel");
    }
    else
    {
        LOG_ERROR("The volume type is not correct");
        return std::string("");
    }
}


}
