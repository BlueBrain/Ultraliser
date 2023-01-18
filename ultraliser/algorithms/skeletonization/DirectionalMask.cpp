#include "DirectionalMask.h"

namespace Ultraliser
{


// Templates for 27-neighborhood and  the  six masks of the up direction.
int8_t ssMu[6][26] =
{{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3 },
 { 2, 2, 2, 0, 0, 0, 0, 0, 0, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 2, 2, 2, 0, 0, 2, 0, 0, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 0, 0, 1, 0, 0, 0, 0, 0, 0, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2 },
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 0, 0, 0, 3, 1, 3, 3, 0, 3, 0, 0, 0 },
 { 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 0, 2, 2, 1, 2, 0, 0, 1, 0, 0, 2 }};

DirectionalMask::DirectionalMask(char direction)
{
    _baseMasks = new Mask[6];

    for (int i = 6; i--;)
    {
        _baseMasks[i].setDirection(direction);
        _baseMasks[i].set_mask_from_u(ssMu[i]);
        _baseMasks[i].generateRotations();
    }
}

DirectionalMask::~DirectionalMask()
{
    delete[]
    _baseMasks;
}


bool DirectionalMask::matches(const int8_t *subVolume)
{
    for (int i = 6; i--;)
    {
        if (_baseMasks[i].matches(subVolume))
            return true;
    }
    return false;
}

}

