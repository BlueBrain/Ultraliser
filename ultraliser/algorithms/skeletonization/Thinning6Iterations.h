#pragma once

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include "DirectionalMask.h"

namespace Ultraliser
{


/**
 * @brief The Thinning6Iterations class
 * An implementation of the 3D thinning algorithm proposed by Kalman Palagyi  and  Attila Kuba.
 * This approach checks if the current voxel in the volume matches a set of six directional masks
 * depending on the given direction.
 * The indices of the faces are set as follows:
 *
 *     0  1  2
 *     3  4  5
 *     6  7  8
 *
 *     9 10 11
 *     12    13
 *     14 15 16
 *
 *     17 18 19
 *     20 21 22
 *     23 24 25
 *
 */
class Thinning6Iterations
{
public:

    /**
     * @brief Thinning6Iterations
     * Constructor.
     */
    Thinning6Iterations();
    ~Thinning6Iterations();

    bool matches(const size_t& direction, const int8_t* subVolume);

private:

    /**
     * @brief _directionalMasks
     * A set of six directional masks that are used to check if the voxel should be removed or not.
     */
    DirectionalMask *_directionalMasks[6];
};

}
