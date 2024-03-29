/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include <geometry/Geometry.h>

namespace Ultraliser
{

class Sample
{
public:

    /**
     * @brief Sample
     * Constructor
     * @param position
     * The Cartesian position of the sample.
     * @param radius
     * The radius of the sample.
     */
    Sample(const Vector3f &position, const float &radius, const size_t& index);

    /**
     * @brief Sample
     * @param sample
     */
    Sample(const Sample* sample);

public:

    /**
     * @brief getPosition
     * @return
     */
    Vector3f getPosition() const;

    /**
     * @brief setPosition
     * @param position
     */
    void setPosition(const Vector3f& position);

    /**
     * @brief getRadius
     * @return
     */
    float getRadius() const;

    /**
     * @brief setRadius
     * @param radius
     */
    void setRadius(const float& radius);

    /**
     * @brief getIndex
     * @return
     */
    size_t getIndex() const;

    /**
     * @brief setIndex
     * @param index
     */
    void setIndex(const size_t index);

    /**
     * @brief isLocatedInBoundingBox
     * Verifies if the sample is located within a given bounding box.
     * @param center
     * The center of the given bounding box.
     * @param width
     * The width of the given bounding box.
     * @param height
     * The height of the given bounding box.
     * @param depth
     * The depth of the given bounding box.
     * @return
     * True if the sample is located in the bounding box, and false otherwise.
     */
    bool isLocatedInBoundingBox(const Vector3f& center,
                                const float& width,
                                const float& height,
                                const float& depth) const;

private:

    /**
     * @brief _position
     * The Cartesian position of the sample.
     */
    Vector3f _position;

    /**
     * @brief _radius
     * The radius of the sample.
     */
    float _radius;

    /**
     * @brief _index
     * The unique index of the sample.
     */
    size_t _index;
};

/**
 * @brief Samples
 * A list of samples.
 */
typedef std::vector< Sample* > Samples;

/**
 * @brief Path
 * The path is composed of a list of samples.
 */
typedef Samples Path;

/**
 * @brief Paths
 * A path is simply a list of samples connecting starting point and an end point, for example, each
 * section in the morphology is a path with a list of connected samples. But this std::vector is
 * made to make it easy to collect all possible combinations of paths along the section starting
 * from any of its parents to any of its children.
 *
 * Each path in the list is a list of connected samples with type Samples or std::vector<Sample>.
 */
typedef std::vector< Path > Paths;
}
