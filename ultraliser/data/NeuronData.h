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

#pragma once

#include <string>
#include <math/Math.h>

namespace Ultraliser
{

/**
 * @brief The Neuron struct
 * A structure that contains all the data that describe a neuron structure
 * and position in a BBP circuit.
 */
struct Neuron
{
    /**
     * @brief gid
     */
    std::string gid;

    /**
     * @brief layer
     */
    std::string layer;

    /**
     * @brief column
     */
    std::string column;

    /**
     * @brief morphologyType
     */
    std::string morphologyType= "";

    /**
     * @brief morphologyLabel
     */
    std::string morphologyLabel = "";

    /**
     * @brief tag
     */
    uint8_t tag;

    /**
     * @brief somaPosition
     */
    Vector3f somaPosition;

    /**
     * @brief xOrientation
     */
    float xOrientation;

    /**
     * @brief yOrientation
     */
    float yOrientation;

    /**
     * @brief zOrientation
     */
    float zOrientation;

    /**
     * @brief somaMeanRadius
     */
    float somaMeanRadius;

    /**
     * @brief somaMinRadius
     */
    float somaMinRadius;

    /**
     * @brief somaMaxRadius
     */
    float somaMaxRadius;

    /**
     * @brief localToGlobalTransform
     */
    Matrix4f localToGlobalTransform;

    /**
     * @brief print
     */
    void print()
    {
        LOG_INFO("* Neuron %s", gid.c_str());
        LOG_INFO("\t Morphology Type: [ %s ]", morphologyType.c_str());
        LOG_INFO("\t Morphology Label: [ %s ]", morphologyLabel.c_str());
        LOG_INFO("\t Tag [%d]", tag);
        LOG_INFO("\t Layer [ %s ]", layer.c_str());
        LOG_INFO("\t Column [ %s ]", column.c_str());
    }
};

/**
 * @brief Neurons
 * A vector (or list) of neurons.
 */
typedef std::vector< Neuron > Neurons;

}
