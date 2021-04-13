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

#ifndef ULTRALISER_MATH_EIGENVECTOR_HH
#define ULTRALISER_MATH_EIGENVECTOR_HH

namespace Ultraliser
{

/**
 * @brief The EigenVector
 */
struct EigenVector
{
    /**
     * @brief x1
     * The x-coordinate of first eigen vector
     */
    float x1;

    /**
     * @brief y1
     * The y-coordinate of first eigen vector.
     */
    float y1;

    /**
     * @brief z1
     * The z-coordinate of first eigen vector.
     */
    float z1;

    /**
     * @brief x2
     * The x-coordinate of second eigen vector.
     */
    float x2;
    /**
     * @brief y2
     * The y-coordinate of second eigen vector.
     */
    float y2;

    /**
     * @brief z2
     * The z-coordinate of second eigen vector.
     */
    float z2;

    /**
     * @brief x3
     * The x-coordinate of third eigen vector.
     */
    float x3;

    /**
     * @brief y3
     * The y-coordinate of third eigen vector.
     */
    float y3;

    /**
     * @brief z3
     * The z-coordinate of third eigen vector.
     */
    float z3;

    /**
     * @brief isValid
     * A bool to indicate if the eigen vector is valid or not.
     */
    bool isValid;
};

}

#endif // ULTRALISER_MATH_EIGENVECTOR_HH
