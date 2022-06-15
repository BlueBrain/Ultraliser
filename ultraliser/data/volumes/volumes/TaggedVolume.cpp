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

#include <data/volumes/volumes/TaggedVolume.h>
#include <data/volumes/grids/VolumeGrid.h>
#include <utilities/Image.h>
#include <data/volumes/grids/Projection.h>

namespace Ultraliser
{

TaggedVolume::TaggedVolume(const Vector3f& pMin,
                           const Vector3f& pMax,
                           const uint64_t &baseResolution,
                           const float &voxelPadding)
    : _pMin(pMin)
    , _pMax(pMax)
{
    // Compute the bounding box of the mesh
    Vector3f bbox = _pMax - _pMin;
    Vector3f delta;

    if (isZero(voxelPadding))
    {
        // If the resolution is greater than 64, use 2.5% padding
        if (baseResolution > 64)
            delta = bbox * 0.025f;

        // Otherwise, use 10%
        else
            delta = bbox * 0.25f;
    }

    // Use the specified voxel padding
    else
        delta = bbox * voxelPadding;

    // Increment the bounding box
    _pMax += delta;
    _pMin -= delta;

    // Compute the bounding box size
    Vector3f boundingBoxSize = (_pMax - _pMin);

    // Find the largest dimension of the mesh model to be able to create a
    // scaled grid.
    int32_t _largestDimensionIdx = getLargestDimension(boundingBoxSize);

    // Compute the voxel size
    const float _voxelSize = (boundingBoxSize[_largestDimensionIdx]) /
            I2F(baseResolution);

    _width = D2UI64(std::ceil(boundingBoxSize[0] / _voxelSize));
    _height = D2UI64(std::ceil(boundingBoxSize[1] / _voxelSize));
    _depth = D2UI64(std::ceil(boundingBoxSize[2] / _voxelSize));

    _numberVoxels = _width * _height * _depth;

    // Allocating the grid
    _allocateVolume();

    // Create the labeling color map
    _createLabelingColorMap();
}

TaggedVolume::TaggedVolume(const uint64_t width,
                           const uint64_t height,
                           const uint64_t depth,
                           Vector3f pMin, Vector3f pMax)
    : _width(width)
    , _height(height)
    , _depth(depth)
    , _pMin(pMin)
    , _pMax(pMax)
    , _numberVoxels(width * height * depth)
{
    _allocateVolume();
    _volumeAdditionTime = 0.f;

    // Create the labeling color map
    _createLabelingColorMap();
}

int32_t TaggedVolume::getLargestDimension(const Vector3f& dimensions)
{
    float value = dimensions[0];
    int index = 0;

    for (int32_t i = 1; i < DIMENSIONS; ++i)
    {
        if (value < dimensions[i])
        {
            index = i;
            value = dimensions[i];
        }
    }

    return index;
}

void TaggedVolume::updateColormap(const uint64_t&)
{

}

uint8_t TaggedVolume::getTag(const uint64_t &index) const
{
    return _data[index];
}

void TaggedVolume::setTag(const uint64_t &index, const uint8_t& tag)
{
    _data[index] = tag;
}

uint8_t TaggedVolume::getTag(const uint64_t &x,
                             const uint64_t &y,
                             const uint64_t &z) const
{
    return getTag(mapToIndex(x, y, z));
}

void TaggedVolume::setTag(const uint64_t &x,
                          const uint64_t &y,
                          const uint64_t &z,
                          const uint8_t tag)
{
    setTag(mapToIndex(x, y, z), tag);
}

uint64_t TaggedVolume::computeNumberNonZeroVoxelsPerSlice(uint64_t z) const
{
    uint64_t numberNonZeroVoxels = 0;

    for (uint64_t i = 0; i < getWidth(); i++)
    {
        for (uint64_t j = 0; j < getHeight(); j++)
        {
            if (_data[mapToIndex(i, j, z)] > 0)
            {
                numberNonZeroVoxels += 1;
            }
        }
    }

    return numberNonZeroVoxels;
}

uint64_t TaggedVolume::computeNumberNonZeroVoxels(void) const
{
    // Starts the timer
    TIMER_SET;

    uint64_t numberNonZeroVoxels = 0;
    uint64_t* numberNonZeroVoxelsPerSlice = new uint64_t[I2UI64(getDepth())];

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getDepth(); i++)
        numberNonZeroVoxelsPerSlice[i] = 0;

    LOOP_STARTS("Computing Filled Voxels")
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getDepth(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
        ++progress;
        if (omp_get_thread_num() == 0)
#else
        ++progress;
#endif
            LOOP_PROGRESS(progress, getDepth());

        numberNonZeroVoxelsPerSlice[i] =
                computeNumberNonZeroVoxelsPerSlice(i);
    }

    for (uint64_t i = 0; i < getDepth(); i++)
    {
        numberNonZeroVoxels += numberNonZeroVoxelsPerSlice[i];
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Delete the arrays
    delete [] numberNonZeroVoxelsPerSlice;

    return numberNonZeroVoxels;
}

float TaggedVolume::computeVolume3()
{
    // Compute the number of non zero voxels
    const uint64_t numberNonZeroVoxels = computeNumberNonZeroVoxels();

    // Get the voxel volume in units3
    // TODO: Update when you have the bounding box data
    const float voxelVolume = 1;

    // Return the result
    return voxelVolume * numberNonZeroVoxels;
}

uint64_t TaggedVolume::mapToIndex(const uint64_t &x,
                                  const uint64_t &y,
                                  const uint64_t &z) const
{
    if (x >= _width  || y >= _height ||z >= _depth)
    {
        return 0;
        // LOG_ERROR("Index out of bound [%d, %d, %d]", x, y, z);
    }

    return (x + (_width * y) + (_width * _height * z));
}

void TaggedVolume::addVolume(const Volume* volume, const uint8_t& index)
{
    // Start timer
    TIMER_SET;

    LOOP_STARTS("Adding Volume");
#ifdef ULTRALISER_USE_OPENMP
    uint64_t progress = 0;
    #pragma omp parallel for
#endif
    for (uint64_t k = 0; k < _depth; ++k)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
        ++progress;

        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _depth);
#else
        LOOP_PROGRESS(k, _depth);
#endif

        for (uint64_t i = 0; i < _width; i++)
        {
            for (uint64_t j = 0; j < _height; j++)
            {
                if (volume->isFilled(i, j, k))
                    _data[mapToIndex(i, j, k)] = index;
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    _volumeAdditionTime = GET_TIME_SECONDS;
}

void TaggedVolume::_writeHeader(const std::string &prefix) const
{
    std::string fileName = prefix + std::string(HEADER_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);
    header << _width << " " << _height << " " << _depth << std::endl;
    header.close();
}

void TaggedVolume::_allocateVolume(void)
{
    _data = new uint8_t[_numberVoxels];

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberVoxels; i++)
        _data[i] = 0;
}

void TaggedVolume::_createLabelingColorMap()
{
    // The color-map that is used for coloring the labeled dataset
    _labelingColorMap.resize(5);
    _labelingColorMap[0] = Vector4f(0.f, 0.f, 0.f, 1.f);
    _labelingColorMap[AXON_INDEX] = Vector4f(0.f, 0.f, 1.f, 1.0);
    _labelingColorMap[BASAL_DENDRITE_INDEX] = Vector4f(1.f, 0.f, 0.f, 1.f);
    _labelingColorMap[APICAL_DENDRITE_INDEX] = Vector4f(0.f, 1.f, 0.f, 1.f);
    _labelingColorMap[SOMA_INDEX] = Vector4f(1.f, 1.f, 0.f, 1.f);
}

void TaggedVolume::writeRAW(const std::string &prefix) const
{
    // Start timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Write the image file
    std::string fileName = prefix + std::string(RAW_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Byte)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel++)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);
        image << _data[voxel];
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void TaggedVolume::writeBIN(const std::string &prefix) const
{
    // Start timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Write the image file
    std::string fileName = prefix + std::string(BINARY_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out | std::ios::binary);

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Bit)");
    for (u_int64_t voxel = 0; voxel < _numberVoxels; voxel += 8)
    {
        uint8_t value = 0;
        for (uint64_t i = 0; i < 8; ++i)
        {
            if (_data[voxel + i])
            {
                value |= 1 << i;
            }
        }
        image << value;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void TaggedVolume::writeASCII(const std::string &prefix) const
{
    // Start timer
    TIMER_SET;

    // Write the header file
    _writeHeader(prefix);

    // Write the image file
    std::string fileName = prefix + std::string(ASCII_EXTENSION);
    std::fstream image;
    image.open(fileName.c_str(), std::ios::out);

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels (1 Bit)");
    for (uint64_t voxel = 0; voxel < _numberVoxels; voxel++)
    {
        image << uint32_t(_data[voxel]);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    image.close();
}

void TaggedVolume::writeNRRD(const std::string &prefix) const
{

    // Starts the timer
    TIMER_SET;

    // File name
    std::string fileName = prefix + std::string(NRRD_EXTENSION);
    FILE* fptr= fopen(fileName.c_str(), "w");

    // Header
    fprintf(fptr, "NRRD0001\n");
    fprintf(fptr, "content: \"Volume\"\n");
    fprintf(fptr, "type: unsigned char\n");
    fprintf(fptr, "dimension: 3\n");
    fprintf(fptr,"sizes: %" PRId64 " %" PRId64 " %" PRId64 "\n",
            getWidth(), getHeight(), getDepth());
    fprintf(fptr, "spacings: 1 1 1\n");
    fprintf(fptr, "encoding: raw\n");

    LOG_STATUS("Exporting Volume [ %s ]", fileName.c_str());

    LOOP_STARTS("Writing Voxels");
    for (int64_t voxel = 0; voxel < _numberVoxels; ++voxel)
    {
        LOOP_PROGRESS_FRACTION(voxel, _numberVoxels);

        fputc(_data[voxel], fptr);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Closing the file
    fclose(fptr);
}

void TaggedVolume::writeVolumes(const std::string &prefix,
                                const bool& binaryFormat,
                                const bool& rawFormat, const bool &nrrdFormat) const
{
    if (binaryFormat || rawFormat)
        LOG_TITLE("Writing Volumes");

    if (binaryFormat)
        writeBIN(prefix);

    if (rawFormat)
        writeRAW(prefix);

    if (nrrdFormat)
        writeNRRD(prefix);
}

void TaggedVolume::project(const std::string &prefix,
                           const bool &xy, const bool &zy, const bool &projectColorCoded) const
{
    if (xy || zy)
    {
        TIMER_SET;

        LOG_TITLE("Projecting Volume");
        LOG_STATUS("Compositing Projection(s)");

        if (xy)
            projectXY(prefix);

        if (zy)
            projectZY(prefix);

        LOG_STATUS_IMPORTANT("Volume Projection Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void TaggedVolume::projectXY(const std::string &prefix, const bool &projectColorCoded) const
{
    // Starts the timer
    TIMER_SET;

    uint8_t* normalizedProjection = new uint8_t[getWidth() * getHeight()];
    double* projection = new double[getWidth() * getHeight()];

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getWidth() * getHeight(); i++)
    {
        normalizedProjection[i] = 0;
        projection[i] = 0.0;
    }

    LOOP_STARTS("XY Projection");
#ifdef ULTRALISER_USE_OPENMP
    uint64_t progress = 0;
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
        ++progress;

        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, getWidth());
#else
        LOOP_PROGRESS(i, getWidth());
#endif

        for (uint64_t k = 0; k < getDepth(); k++)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                bool bit = _data[mapToIndex(i, j, k)] > 0;
                if (bit) projection[i + getWidth() * j] += 1.0;
            }
        }
    }
    LOOP_DONE;

    // Normalize the value
    double maxValue = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getWidth() * getHeight(); i++)
    {
        if (projection[i] > maxValue)
            maxValue = projection[i];
    }

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getWidth() * getHeight(); i++)
    {
        double pixelValue = 255.0 * projection[i] / maxValue;
        normalizedProjection[i] = uint8_t(pixelValue);
    }

    // Save the projection
    std::stringstream stream;
    stream << prefix << "_projection.xy";

    // Save PPM image
    Utilities::savePPMLuminanceImage(stream.str(), normalizedProjection,
                            UI2I64(getWidth()), UI2I64(getHeight()));

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Save color coded projections with all possible color-maps
    if (projectColorCoded)
    {
        saveColorMappedProjectionWithAllColorMaps(prefix, projection,
                                                  getWidth(), getHeight(), 0, maxValue);
    }

    // Release the array
    delete [] normalizedProjection;
    delete [] projection;
}

void TaggedVolume::projectZY(const std::string &prefix, const bool &projectColorCoded) const
{
    // Starts the timer
    TIMER_SET;

    uint8_t* normalizedProjection = new uint8_t[getDepth() * getHeight()];
    double* projection = new double[getDepth() * getHeight()];

#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getDepth() * getHeight(); i++)
    {
        normalizedProjection[i] = 0;
        projection[i] = 0.0;
    }

    LOOP_STARTS("ZY Projection");
#ifdef ULTRALISER_USE_OPENMP
    uint64_t progress = 0;
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
        ++progress;

        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, getWidth());
#else
        LOOP_PROGRESS(i, getWidth());
#endif

        for (uint64_t k = 0; k < getDepth(); k++)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                bool bit = _data[mapToIndex(i, j, k)] > 0;
                if (bit) projection[k + getDepth() * j] += 1.0;
            }
        }
    }
    LOOP_DONE;

    // Normalize the value
    double maxValue = 0.0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getDepth() * getHeight(); i++)
    {
        if (projection[i] > maxValue)
            maxValue = projection[i];
    }

#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for schedule(dynamic, 1)
#endif
    for (uint64_t i = 0; i < getDepth() * getHeight(); i++)
    {
        double pixelValue = 255.0 * projection[i] / maxValue;
        normalizedProjection[i] = uint8_t(pixelValue);
    }

    // Save the projection
    std::stringstream stream;
    stream << prefix << "_projection.zy";

    // Write a PPM image
    Utilities::savePPMLuminanceImage(stream.str(), normalizedProjection,
                            UI2I64(getDepth()), UI2I64(getHeight()));

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Save color coded projections with all possible color-maps
    if (projectColorCoded)
    {
        saveColorMappedProjectionWithAllColorMaps(prefix, projection,
                                                  getDepth(), getHeight(), 0, maxValue);
    }

    // Release the array
    delete [] normalizedProjection;
    delete [] projection;
}


void TaggedVolume::writeStackXY(const std::string &outputDirectory,
                                  const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "xy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: XY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t z = 0; z < getDepth(); ++z)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getDepth());

        // Create a slice
        auto slice = std::make_unique<Image>(getWidth(), getHeight());

        for (uint64_t i = 0; i < getWidth(); i++)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                uint64_t index = mapToIndex(i, j, z);

                if (_data[index] > 0)
                    slice->setPixelColor(i , j, WHITE);
                else
                    slice->setPixelColor(i , j, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << z;
        slice->writePPM(stream.str());
    }

    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);
}

void TaggedVolume::writeStackZY(const std::string &outputDirectory,
                                  const std::string &prefix) const
{
    // Starts the timer
    TIMER_SET;

    // Make a directory
    std::string directoryName = outputDirectory + "/" + "zy-" + prefix;
    mkdir(directoryName.c_str(), 0777);

    LOG_STATUS("Exporting Image Stack [ %s ]", directoryName.c_str());

    LOOP_STARTS("Writing Stack: ZY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        auto slice = std::make_unique<Image>(getDepth(), getHeight());

        for (uint64_t z = 0; z < getDepth(); z++)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                uint64_t index = mapToIndex(i, j, z);

                if (_data[index] > 0)
                    slice->setPixelColor(z , getHeight() - j - 1, WHITE);
                else
                    slice->setPixelColor(z , getHeight() - j - 1, BLACK);
            }
        }

        std::stringstream stream;
        stream << directoryName << "/" << i;
        slice->writePPM(stream.str());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);
}

void TaggedVolume::writeStacks(const std::string &outputDirectory,
                               const std::string &prefix,
                               const bool& xy,
                               const bool& zy) const
{
    if (xy || zy)
        LOG_TITLE("Writing Stacks");

    if (xy)
        writeStackXY(outputDirectory, prefix);

    if (zy)
        writeStackZY(outputDirectory, prefix);
}

void TaggedVolume::composeBrainbowXY(const std::string &prefix,
                                     const std::vector<Vector4f> colors,
                                     const float alpha) const
{
    // Start the timer
    TIMER_SET;
    Vector4f* projection = new Vector4f[getWidth() * getHeight()];

#pragma omp parallel for
    for (size_t i = 0; i < getWidth() * getHeight(); i++)
        projection[i] = Vector4f();

    LOOP_STARTS("Composing Tagged Image: XY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        for (uint64_t k = 0; k < getDepth(); ++k)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                const uint8_t voxelTag = _data[mapToIndex(i, j, k)];

                // Source voxel
                Vector4f src;
                if (voxelTag > 0)
                {
                    src = colors[voxelTag - 1];
                    Vector4f dst = projection[i + getWidth() * j];
                    projection[i + getWidth() * j] =
                            src * alpha + (dst * (1 - alpha));
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Save the projection
    std::stringstream stream;
    stream << prefix << ".xy_tagged";

    // Write a PPM image
    Utilities::savePPMColoredImage(stream.str(), projection,
                                   UI2I64(getWidth()), UI2I64(getHeight()));

    // Release the array
    delete [] projection;
}

void TaggedVolume::composeBrainbowZY(const std::string &prefix,
                                     const std::vector< Vector4f > colors,
                                     const float alpha) const
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Composing Tagged Image: ZY");
    Vector4f* projection = new Vector4f[getDepth() * getHeight()];

#pragma omp parallel for
    for (size_t i = 0; i < getDepth() * getHeight(); i++)
        projection[i] = Vector4f();

    LOOP_STARTS("Composing Labeled Image: ZY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        for (uint64_t k = 0; k < getDepth(); ++k)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                const uint8_t voxelTag = _data[mapToIndex(i, j, k)];

                // Source voxel
                Vector4f src;
                if (voxelTag > 0)
                {
                    src = colors[voxelTag - 1];
                    Vector4f dst = projection[k + getDepth() * j];
                    projection[k + getDepth() * j] =
                            src * alpha + (dst * (1 - alpha));
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Save the projection
    std::stringstream stream;
    stream << prefix << ".zy_tagged";

    Utilities::savePPMColoredImage(stream.str(), projection,
                                   UI2I64(getDepth()), UI2I64(getHeight()));

    // Release the array
    delete [] projection;
}

void TaggedVolume::composeBrainbow(const std::string &prefix,
                                   const std::vector< Vector4f > colors,
                                   const bool &xy, const bool &zy,
                                   const float alpha) const
{
    if (xy || zy)
    {
        TIMER_SET;

        LOG_TITLE("Projecting Brainbow Volume ");
        LOG_STATUS("Compositing Tagged Projection(s)");

        if (xy)
            composeBrainbowXY(prefix, colors, alpha);

        if (zy)
            composeBrainbowZY(prefix, colors, alpha);

        LOG_STATUS_IMPORTANT("Tagged Volume Projection Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}

void TaggedVolume::composeLabeledImageXY(const std::string &prefix,
                                         const float alpha) const
{
    // Start the timer
    TIMER_SET;

    Vector4f* projection = new Vector4f[getWidth() * getHeight()];
    float maxValue = 0.f;

    #pragma omp parallel for
    for (size_t i = 0; i < getWidth() * getHeight(); i++)
        projection[i] = Vector4f();


    LOOP_STARTS("Composing Labeled Image: XY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        for (uint64_t k = 0; k < getDepth(); ++k)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                const uint8_t voxelTag = _data[mapToIndex(i, j, k)];

                // Source voxel
                Vector4f src;
                if (voxelTag > 0)
                {
                    src = _labelingColorMap[voxelTag];
                    Vector4f dst = projection[i + getWidth() * j];
                    Vector4f value = src * alpha + (dst * (1 - alpha));
                    projection[i + getWidth() * j] = value;

                    // For normalization
                    if (value.abs() > maxValue)
                        maxValue = value.abs();
                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Normalization
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {
        for (uint64_t j = 0; j < getHeight(); j++)
        {
            projection[i + getWidth() * j] =
                    projection[i + getWidth() * j] / maxValue;
        }
    }

    // Save the projection
    std::stringstream stream;
    stream << prefix << ".xy_labeled";

    // Write a PPM image
    Utilities::savePPMColoredImage(stream.str(), projection,
                                   UI2I64(getWidth()), UI2I64(getHeight()));


    // Release the array
    delete [] projection;
}

void TaggedVolume::composeLabeledImageZY(const std::string &prefix,
                                         const float alpha) const
{
    // Start the timer
    TIMER_SET;

    Vector4f* projection = new Vector4f[getDepth() * getHeight()];
    float maxValue = 0.f;

    #pragma omp parallel for
    for (size_t i = 0; i < getDepth() * getHeight(); i++)
        projection[i] = Vector4f();

    LOOP_STARTS("Composing Labeled Image: ZY");
    int64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t i = 0; i < getWidth(); i++)
    {
#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getWidth());

        for (uint64_t k = 0; k < getDepth(); ++k)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                const uint8_t voxelTag = _data[mapToIndex(i, j, k)];

                // Source voxel
                Vector4f src;
                if (voxelTag > 0)
                {
                    src = _labelingColorMap[voxelTag];
                    Vector4f dst = projection[k + getDepth() * j];
                    Vector4f value = src * alpha + (dst * (1 - alpha));
                    projection[k + getDepth() * j] = value;

                    // For normalization
                    if (value.abs() > maxValue)
                        maxValue = value.abs();

                }
            }
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Normalization
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t k = 0; k < getDepth(); ++k)
    {
        for (uint64_t j = 0; j < getHeight(); j++)
        {
            projection[k + getDepth() * j] =
                    projection[k + getDepth() * j] / maxValue;
        }
    }

     // Save the projection
    std::stringstream stream;
    stream << prefix << ".zy_labeled";

    Utilities::savePPMColoredImage(stream.str(), projection,
                                   UI2I64(getDepth()), UI2I64(getHeight()));

    // Release the array
    delete [] projection;
}

void TaggedVolume::composeLabeledProjetions(const std::string &prefix,
                                            const bool &xy, const bool &zy,
                                            const float alpha) const
{
    if (xy || zy)
    {
        TIMER_SET;

        LOG_TITLE("Projecting Labeled Volume ");
        LOG_STATUS("Compositing Labeled Projection(s)");

        if (xy)
            composeLabeledImageXY(prefix, alpha);

        if (zy)
            composeLabeledImageZY(prefix, alpha);

        LOG_STATUS_IMPORTANT("Labeled Volume Projection Stats.");
        LOG_STATS(GET_TIME_SECONDS);
    }
}


Volume* TaggedVolume::getPlainVolume()
{
    Volume* plainVolume = new Volume(UI2I64(_width),
                                     UI2I64(_height),
                                     UI2I64(_depth),
                                     _pMin, _pMax);
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
#pragma omp parallel for
#endif
    for (uint64_t z = 0; z < getDepth(); ++z)
    {

#ifdef ULTRALISER_USE_OPENMP
        #pragma omp atomic
#endif
        ++progress;

#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
#endif
            LOOP_PROGRESS(progress, getDepth());

        for (uint64_t i = 0; i < getWidth(); i++)
        {
            for (uint64_t j = 0; j < getHeight(); j++)
            {
                uint64_t index = mapToIndex(i, j, z);
                if (_data[index] > 0)
                    plainVolume->fill(index);
            }
        }
    }

    return plainVolume;
}

void TaggedVolume::printVolumeStats(const std::string &reference,
                                    const std::string *prefix)
{
    LOG_TITLE("Volume Statistics");

    LOG_STATUS("Collecting Stats.");
    Vector3f bounds = _pMax - _pMin;
    float volumeSize = computeVolume3();

    // Write the statistics to a file
    if (prefix != nullptr)
    {
        // Create the file
        std::string fileName = *prefix + "-" + reference + VOLUME_INFO_EXTENSION;
        LOG_STATUS("Writing Info. [ %s ] \n", fileName.c_str());

        FILE* info = fopen(fileName.c_str(), "w");
        fprintf(info, "Stats. [ %s ] \n", reference.c_str());

        if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
        {
            fprintf(info, "\t* Bounding Box:         | [%f, %f, %f] \n",
                     F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
            fprintf(info, "\t* pMin:                 | [%f, %f, %f] \n",
                     F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
            fprintf(info, "\t* pMax:                 | [%f, %f, %f] \n",
                     F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
        }

        fprintf(info, "\t* Resolution            | [%d] x [%d] x [%d] \n",
                 I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
        fprintf(info, "\t* Number of Voxels      | %" PRIu64 " \n",
                getNumberVoxels());
        fprintf(info, "\t* Volume Format         | Byte \n");
        fprintf(info, "\t* Size in Memory        | %sBytes \n",
                 FORMAT(getNumberBytes()));
        fprintf(info, "\t* Volume                | %f³ \n",
                 F2D(volumeSize));

        // Close the file
        fclose(info);
    }

    LOG_STATUS("Volume [ %s ]", reference.c_str());

    if (bounds.x() > 0.f || bounds.y() > 0.f || bounds.z() > 0.f)
    {
        LOG_INFO("\t* Bounding Box:         | [%f, %f, %f]",
                 F2D(bounds.x()), F2D(bounds.y()), F2D(bounds.z()));
        LOG_INFO("\t* pMin:                 | [%f, %f, %f]",
                 F2D(_pMin.x()), F2D(_pMin.y()), F2D(_pMin.z()));
        LOG_INFO("\t* pMax:                 | [%f, %f, %f]",
                 F2D(_pMax.x()), F2D(_pMax.y()), F2D(_pMax.z()));
    }

    LOG_INFO("\t* Resolution            | [%d] x [%d] x [%d]",
             I2I32(getWidth()), I2I32(getHeight()), I2I32(getDepth()));
    LOG_INFO("\t* Number of Voxels      | %" PRIu64 "",
             getNumberVoxels());
    LOG_INFO("\t* Volume Format         | Byte");
    LOG_INFO("\t* Size in Memory        | %sBytes",
             FORMAT(getNumberBytes()));
    LOG_INFO("\t* Volume                | %f³",
             F2D(volumeSize));
}

TaggedVolume::~TaggedVolume()
{
    delete [] _data;
}

}
