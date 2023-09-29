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

#include "Mask.h"
#include "Neighbors.hh"
#include "RotationMatrices.hh"
#include <common/Headers.hh>

#define NEIGHBOURING_VOXELS 26

namespace Ultraliser
{

Mask::Mask()
{
    // Allocate the masks
    _mask0      = new int8_t[ NEIGHBOURING_VOXELS ];
    _mask90     = new int8_t[ NEIGHBOURING_VOXELS ];
    _mask180    = new int8_t[ NEIGHBOURING_VOXELS ];
    _mask270    = new int8_t[ NEIGHBOURING_VOXELS ];

    // Initialize all the masks with -1
    for (size_t i = 0; i < NEIGHBOURING_VOXELS; ++i)
    {
        _mask0[i] = _mask90[i] = _mask180[i] = _mask270[i] = -1;
    }

    // Initially, unspecified direction
    _direction = UNSPECIFIED;
    direction = ' ';
}

Mask::~Mask()
{
    // Delete all the masks
    delete[] _mask0;
    delete[] _mask90;
    delete[] _mask180;
    delete[] _mask270;
}

void Mask::setDirection(char d)
{
    direction = d;
}

void Mask::updateDirection(DIRECTION newDirection)
{
    _direction = newDirection;
}

void Mask::printDirection()
{
    printf("%c\n",direction);
}

/// Uses the umask in  order  to directions. Indeed, to generate  the  "n" version of mask M2,
/// showed in he  figure 3 of the paper, we just to do a rotation  of M2  about  the axis  x.
/// To  perform this transformation, we use:
/// >>> rotate(umask,vector,'x'), which rotates the mask umask about "x" and stores the result
/// in "vector".
void Mask::set_mask_from_u(int8_t umask[26])
{
    if (direction == ' ') return;

    rotate(umask, _mask0,'x');

    if (direction == 'n') return;

    rotate(_mask0, _mask0,'x');

    if (direction == 'd') return;

    rotate(_mask0, _mask0,'x');

    if (direction == 's') return;

    rotate(_mask0, _mask0,'x');

    if (direction == 'u') return;

    rotate(_mask0, _mask0,'z');

    if (direction == 'w') return;

    rotate(_mask0, _mask0,'z');
    rotate(_mask0, _mask0,'z');

    if (direction == 'e') return;
}

/// Rotations for the  directions (u,d), (n,s) and  (w,e) must be  performed  about the same
void Mask::generateRotations()
{
    char axis;
    if (direction == ' ') return;

    if (direction == 'u' || direction == 'd') axis = 'y';
    if (direction == 'n' || direction == 's') axis = 'z';
    if (direction == 'w' || direction == 'e') axis = 'x';

    // 90 rotation
    rotate(_mask0, _mask90, axis);

    // 90 rotation
    rotate(_mask90, _mask180, axis);

    // 90 rotation
    rotate(_mask180, _mask270, axis);
}

void Mask::printMask()
{
    for (int i=0;i<26;i++) printf("%d ",_mask0[i]);
    printf("\n");
}

void Mask::printMasks() {
    printDirection();
    printf("0\t=\t");
    for (int i=0;i<26;i++) printf("%d ",_mask0[i]);
    printf("\n");
    printf("90\t=\t");
    for (int i=0;i<26;i++) printf("%d ",_mask90[i]);
    printf("\n");
    printf("180\t=\t");
    for (int i=0;i<26;i++) printf("%d ",_mask180[i]);
    printf("\n");
    printf("270\t=\t");
    for (int i=0;i<26;i++) printf("%d ",_mask270[i]);
    printf("\n");
}

void Mask::rotate(int8_t input[26], int8_t *output,char axis)
{
    // Get the rotation matrix MR based on the axis;
    int* MR;
    if (axis=='x')
        MR = ROTATION_MATRIX_X;

    if (axis=='y')
        MR = ROTATION_MATRIX_Y;

    if (axis=='z')
        MR = ROTATION_MATRIX_Z;

    // Auxiliary matrix to store the indices after the rotation
    int aux[3][3][3];

    // OMP_PARALLEL_FOR
    for (size_t i = 0; i < 26; ++i)
    {
        // Get the index after the rotation operation (3D rotation with 90 degrees)
        const int x = MR[0] * VDX[i] + MR[1] * VDY[i] + MR[2] * VDZ[i];
        const int y = MR[3] * VDX[i] + MR[4] * VDY[i] + MR[5] * VDZ[i];
        const int z = MR[6] * VDX[i] + MR[7] * VDY[i] + MR[8] * VDZ[i];

        // Update the auxiliary matrix
        aux[x+1][y+1][z+1] = input[i];
    }

    // Update the output vector
    // TODO: OMP_PARALLEL_FOR
    for (size_t i = 0; i < NEIGHBOURING_VOXELS; ++i)
        output[i] = aux[VDX[i] + 1][VDY[i] + 1][VDZ[i] + 1];
}

// There are four types of indexes in the masks: 0, 1, 2 and 3. We need worry
// about the indexes 0, 1 and 3 bacause 2 indicates "does not matter", ie the
// value of this  voxel  in the  volume can be any.  Indexes 0 and 1  must be
// exactly the same as in the  mask as in the volume. Index 3 is a particular
// case. As described in the paper, at last one voxel of the volume with this
// index must be 1.

bool Mask::matches(const int8_t *vol, int8_t *vec)
{
    int v, q3 = 0;
    int i;

    for (uint8_t i = 0; i < 26; ++i)
    {
        v = vol[i];

        if (vec[i] == 0 && v!= 0)
            return false;

        if (vec[i] ==1 && v!= 1)
            return false;

        if (vec[i] == 3)
            q3 = (v == 1 || q3 == 2) ? 2 : 1;
    }

    // If  q3  remains  equal  to 0,
    // there is no values 3  in  the
    // the mask. If q3 == 2 there is
    // value 3 in  the mask  and  at
    // last one  voxel in the volume
    // has value 1.

    if (q3 == 1)
        return false;
    return true;
}

// Verify if the rotations match.
bool Mask::matches(const int8_t* subVolume)
{
    if (matches(subVolume, _mask0))
        return true;

    if (matches(subVolume, _mask90))
        return true;

    if (matches(subVolume, _mask180))
        return true;

    if (matches(subVolume, _mask270))
        return true;

    return false;
}

}
