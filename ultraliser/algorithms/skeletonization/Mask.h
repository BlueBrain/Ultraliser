/***************************************************************************************************
 * Copyright (c) 2016 - 2023
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

#include <common/Common.h>

namespace Ultraliser
{

/**
 * @brief The Mask class
 *  Class  used  to  store a mask and its  rotations (90, 180 and 270 degrees)
 */
class Mask
{
public:

    /**
     * @brief The AXIS enum
     */
    enum AXIS { X, Y, Z };

    /**
     * @brief The DIRECTION enum
     */
    enum DIRECTION { U, D, E, W, N, S, UNSPECIFIED };

    Mask();
    ~Mask();

    void printMask();
    void printMasks();

    void generateRotations();					 // Generates  the rotated  masks
    void printDirection();

    void setDirection(char d);

    /**
     * @brief updateDirection
     * @param newDirection
     */
    void updateDirection(DIRECTION newDirection);

    // Tests if one of the four masks  (0, 90, 180 and 270) matchs with the neighborhood of the
    // voxel p  in the volume Vol
    bool matches(const int8_t *subVolume);


    // Generates the four masks from the mask "umask" in up direction
    void set_mask_from_u(int8_t umask[]);

private:

    int8_t *_mask0;
    int8_t *_mask90;
    int8_t *_mask180;
    int8_t *_mask270;

    char direction;

    DIRECTION _direction;

    // Static functions to match and rotate the vectors
    static bool matchf(int ***Vol, int *vec);
    bool matchf1d(int*vol, int *vec);
    bool matches(const int8_t *vol, int8_t *vec);



    static void rotate(int8_t vector[], int8_t *vector_rot, char axis);
};

}
