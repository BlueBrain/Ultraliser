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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_SAMPLE_H
#define ULTRALISER_DATA_MORPHOLOGIES_SAMPLE_H

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
    Sample(const Vector3f &position, const float &radius, const uint64_t& index);

public:

    /**
     * @brief getPosition
     * @return
     */
    Vector3f getPosition() const;

    /**
     * @brief getRadius
     * @return
     */
    float getRadius() const;

    /**
     * @brief getIndex
     * @return
     */
    uint64_t getIndex() const;

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
    uint64_t _index;
};

/**
 * @brief Samples
 * A list of samples.
 */
typedef std::vector< Sample* > Samples;

/**
 * @brief Paths
 * A path is simply a list of samples connecting starting point and an end point, for example, each
 * section in the morphology is a path with a list of connected samples. But this std::vector is
 * made to make it easy to collect all possible combinations of paths along the section starting
 * from any of its parents to any of its children.
 *
 * Each path in the list is a list of connected samples with type Samples or std::vector<Sample>.
 */
typedef std::vector<Samples> Paths;

}
#endif // ULTRALISER_DATA_MORPHOLOGIES_SAMPLE_H
