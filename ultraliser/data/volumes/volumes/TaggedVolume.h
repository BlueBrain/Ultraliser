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

#ifndef ULTRALISER_TAGGED_VOLUME_H
#define ULTRALISER_TAGGED_VOLUME_H

#include <common/Common.h>
#include <data/volumes/volumes/Volume.h>

namespace Ultraliser
{

/**
 * @brief The TaggedVolume class
 */
class TaggedVolume
{
public:

    /**
     * @brief TaggedVolume
     * Constructs an 8-bit volume that can use 255 tags to label either the
     * different components of the neurons or different neurons loaded from a
     * circuit file with multiple tags.
     * @param pMin
     * pMin of the bounding box of the volume.
     * @param pMax
     * pMax of the bounding box of the volume.
     * @param baseResolution
     * The base resolution of the volume.
     * @param voxelPadding
     * A little value that is added to zero-pad the volume to avoid
     * intersections at the edge.
     */
    TaggedVolume(const Vector3f& pMin,
                 const Vector3f& pMax,
                 const uint64_t &baseResolution = 512,
                 const float &voxelPadding = 0.0);

    /**
     * @brief TaggedVolume
     * Constructs an 8-bit volume that can use 255 tags to label either the
     * different components of the neurons or different neurons loaded from a
     * circuit file with multiple tags.
     * @param width
     * Volume width in voxels.
     * @param height
     * Volume height in voxels.
     * @param depth
     * Volume depth in voxels.
     * @param pMin
     * pMin of the bounding box of the volume.
     * @param pMax
     * pMax of the bounding box of the volume.
     */
    TaggedVolume(const uint64_t width,
                 const uint64_t height,
                 const uint64_t depth,
                 Vector3f pMin, Vector3f pMax);
    ~TaggedVolume();

public:

    int32_t getLargestDimension(const Vector3f& dimensions);

    /**
     * @brief composeBrainbowXY
     * @param prefix
     * @param colors
     * @param alpha
     */
    void composeBrainbowXY(const std::string &prefix,
                           const std::vector< Vector4f > colors,
                           const float alpha = 0.1f) const;

    /**
     * @brief composeBrainbowZY
     * @param prefix
     * @param colors
     * @param alpha
     */
    void composeBrainbowZY(const std::string &prefix,
                           const std::vector< Vector4f > colors,
                           const float alpha = 0.1f) const;

    /**
     * @brief composeBrainbow
     * @param prefix
     * @param colors
     * @param alpha
     * @param xy
     * @param zy
     */
    void composeBrainbow(const std::string &prefix,
                         const std::vector< Vector4f > colors,
                         const bool& xy, const bool& zy,
                         const float alpha = 0.1f) const;

    void composeLabeledImageXY(const std::string &prefix,
                               const float alpha = 0.1f) const;

    /**
     * @brief composeLabeledImageZY
     * @param prefix
     * @param alpha
     */
    void composeLabeledImageZY(const std::string &prefix,
                               const float alpha) const;

    /**
     * @brief composeLabeledProjetions
     * @param prefix
     * @param colors
     * @param xy
     * @param zy
     * @param alpha
     */
    void composeLabeledProjetions(const std::string &prefix,
                                  const bool &xy, const bool &zy,
                                  const float alpha = 0.1f) const;

    /**
     * @brief printVolumeStats
     * @param reference
     * @param prefix
     */
    void printVolumeStats(const std::string &reference,
                          const std::string *prefix);

    /**
     * @brief updateColormap
     * @param length
     */
    void updateColormap(const uint64_t &length);

    /**
     * @brief getWidth
     * @return
     */
    uint64_t getWidth(void) const
    {
        return _width;
    }

    /**
     * @brief getHeight
     * @return
     */
    uint64_t getHeight(void) const
    {
        return _height;
    }

    /**
     * @brief getDepth
     * @return
     */
    uint64_t getDepth(void) const
    {
        return _depth;
    }

    /**
     * @brief getNumberVoxels
     * @return
     */
    uint64_t getNumberVoxels(void) const
    {
        return _width * _height * _depth;
    }

    /**
     * @brief computeNumberNonZeroVoxelsPerSlice
     * @param z
     * @return
     */
    uint64_t computeNumberNonZeroVoxelsPerSlice(uint64_t z) const;

    /**
     * @brief computeNumberNonZeroVoxels
     * @return
     */
    uint64_t computeNumberNonZeroVoxels(void) const;

    /**
     * @brief computeVolume3
     * @return
     */
    float computeVolume3();

    /**
     * @brief getNumberBytes
     * @return
     */
    uint64_t getNumberBytes(void) const
    {
        return getNumberVoxels();
    }

    /**
     * @brief getTag
     * @param index
     * @return
     */
    uint8_t getTag(const uint64_t &index) const;

    /**
     * @brief setTag
     * @param index
     * @param tag
     */
    void setTag(const uint64_t &index, const uint8_t &tag);

    /**
     * @brief getTag
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint8_t getTag(const uint64_t &x,
                   const uint64_t &y,
                   const uint64_t &z) const;

    /**
     * @brief setTag
     * @param x
     * @param y
     * @param z
     * @param tag
     */
    void setTag(const uint64_t &x,
                const uint64_t &y,
                const uint64_t &z,
                const uint8_t tag);

    /**
     * @brief mapToIndex
     * @param x
     * @param y
     * @param z
     * @return
     */
    uint64_t mapToIndex(const uint64_t &x,
                        const uint64_t &y,
                        const uint64_t &z) const;

    /**
     * @brief addVolume
     * @param volume
     */
    void addVolume(const Volume* volume, const uint8_t& index);

    /**
     * @brief writeRAW
     * @param prefix
     */
    void writeRAW(const std::string &prefix) const;

    /**
     * @brief writeVolumes
     * @param prefix
     * @param binaryFormat
     * @param rawFormat
     */
    void writeVolumes(const std::string &prefix,
                      const bool& binaryFormat = false,
                      const bool& rawFormat = false,
                      const bool& nrrdFormat = false) const;

    /**
     * @brief writeBIN
     * @param prefix
     */
    void writeBIN(const std::string &prefix) const;

    /**
     * @brief writeASCII
     * @param prefix
     */
    void writeASCII(const std::string &prefix) const;

    /**
     * @brief writeNRRD
     * @param prefix
     */
    void writeNRRD(const std::string &prefix) const;

    /**
     * @brief projectVolume
     * @param prefix
     * @param xy
     * @param zy
     * @param projectColorCoded
     */
    void project(const std::string &prefix,
                 const bool &xy = true,
                 const bool &zy = false,
                 const bool &projectColorCoded = false) const;

    /**
     * @brief projectXY
     * @param prefix
     */
    void projectXY(const std::string &prefix, const bool &projectColorCoded = false) const;

    /**
     * @brief projectZY
     * @param prefix
     */
    void projectZY(const std::string &prefix, const bool &projectColorCoded = false) const;

    /**
     * @brief writeStacks
     * @param outputDirectory
     * @param prefix
     * @param xy
     * @param zy
     */
    void writeStacks(const std::string &outputDirectory,
                     const std::string &prefix,
                     const bool& xy = false,
                     const bool& zy = false) const;

    /**
     * @brief writeStackXY
     * @param prefix
     */
    void writeStackXY(const std::string &outputDirectory,
                        const std::string &prefix) const;

    /**
     * @brief writeStackZY
     * @param outputDirectory
     * @param prefix
     */
    void writeStackZY(const std::string &outputDirectory,
                        const std::string &prefix) const;

    /**
     * @brief getVolumeAdditionTime
     */
    float getVolumeAdditionTime(void) const
    {
        return _volumeAdditionTime;
    }

    /**
     * @brief getPlainVolume
     * Returns a plain volume (1-bit per voxel) that can be used for meshing.
     * @return
     */
    Volume* getPlainVolume();

private:

    /**
     * @brief _createLabelingColorMap
     * Creates the labeling color map that is used to create a color-coded
     * projection image of the volume.
     */
    void _createLabelingColorMap();

    /**
     * @brief _writeHeader
     * @param prefix
     */
    void _writeHeader(const std::string &prefix) const;

    /**
     * @brief _allocateVolume
     */
    void _allocateVolume(void);

private:

    /**
     * @brief _width
     */
    uint64_t _width;

    /**
     * @brief _height
     */
    uint64_t _height;

    /**
     * @brief _depth
     */
    uint64_t _depth;

    /**
     * @brief _pMin
     */
    Vector3f _pMin;

    /**
     * @brief _pMax
     */
    Vector3f _pMax;

    /**
     * @brief _numberVoxels
     */
    uint64_t _numberVoxels;

    /**
     * @brief _data
     */
    uint8_t* _data;

    /**
     * @brief _volumeAdditionTime
     */
    double _volumeAdditionTime;

    /**
     * @brief _colormap
     */
    std::vector< Vector4f > _brainbowColorMap;

    /**
     * @brief _labelingColorMap
     */
    std::vector< Vector4f > _labelingColorMap;
};

}

#endif // ULTRALISER_TAGGED_VOLUME_H
