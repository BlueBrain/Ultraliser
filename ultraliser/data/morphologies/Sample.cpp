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

#include "Sample.h"

namespace Ultraliser
{

Sample::Sample(const Vector3f &position, const float &radius, const size_t &index)
    : _position(position)
    , _radius(radius)
    , _index(index)
{
    /// EMPTY CONSTRUCTOR
}

Sample::Sample(const Sample* sample)
{
    _position = sample->getPosition();
    _radius = sample->getRadius();
    _index = sample->getIndex();
}

Ultraliser::Vector3f Sample::getPosition() const
{
    return _position;
}

float Sample::getRadius() const
{
    return _radius;
}

size_t Sample::getIndex() const
{
    return _index;
}

void Sample::setIndex(const size_t index)
{
    _index = index;
}

bool Sample::isLocatedInBoundingBox(const Vector3f& center,
                                    const float& width,
                                    const float& height,
                                    const float& depth) const
{
    const auto& xMax = center.x() + (width * 0.5);
    const auto& xMin = center.x() - (width * 0.5);

    const auto& yMax = center.y() + (height * 0.5);
    const auto& yMin = center.y() - (height * 0.5);

    const auto& zMax = center.z() + (depth * 0.5);
    const auto& zMin = center.z() - (depth * 0.5);

    // X verification
    if (_position.x() > xMax || _position.x() < xMin)
    {
        return false;
    }

    // Y verification
    if (_position.y() > yMax || _position.y() < yMin)
    {
        return false;
    }

    // Z verification
    if (_position.z() > zMax || _position.z() < zMin)
    {
        return false;
    }

    // All fine, return True
    return true;
}

}
