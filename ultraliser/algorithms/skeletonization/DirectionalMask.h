#pragma once 

#include "Mask.h"

namespace Ultraliser
{


/**
 * @brief The DirectionalMask class
 * Store M1, M2, M3, M4, M5 and M6 masks for a given direction.
 * This directional mask contains the base mask and and the four rotations.
 * This class then includes 6 x 4 = 24 masks for the same direction.
 * NOTE: The i-th element of the base mask corresponds to the mask Mi
 */
class DirectionalMask
{
public:

    /**
     * @brief DirectionalMask
     * @param direction
     */
    DirectionalMask(char direction);
    ~DirectionalMask();

    /**
     * @brief matches
     * Does the mask match an input sub-volume?
     * @param subVolume
     * @return
     */
    bool matches(const int8_t *subVolume);

private:

    /**
     * @brief _baseMasks
     */
    Mask *_baseMasks;
};

}
