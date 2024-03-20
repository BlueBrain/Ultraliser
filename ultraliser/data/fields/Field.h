/***************************************************************************************************
 * Copyright (c) 2016 - 2024
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <math/Vector3f.h>
#include <vector>

namespace Ultraliser
{

/**
 * @brief The Field class
 */
class Field
{
public:

    /**
     * @brief Field
     * @param center
     * @param radius
     */
    Field(const Vector3f center, const float radius)
    {
        _center = center;
        _radius = radius;
    }

public:

    /**
     * @brief getCenter
     * @return
     */
    const Vector3f getCenter() const { return _center; }

    /**
     * @brief getRadius
     * @return
     */
    const float getRadius() const { return _radius; }

private:

    /**
     * @brief _center
     */
    Vector3f _center;

    /**
     * @brief _radius
     */
    float _radius;
};

/**
 * @brief FieldPtr
 */
typedef Field* FieldPtr;

/**
 * @brief Fields
 */
typedef std::vector< Field > Fields;

/**
 * @brief FieldPtrs
 */
typedef std::vector< FieldPtr > FieldPtrs;

}

