#include "Thinning6Iterations.h"

namespace Ultraliser
{


Thinning6Iterations::Thinning6Iterations()
{
    // Generates the 24 masks for each direction (6x24 = 144)
    _directionalMasks[0] = new DirectionalMask('u');
    _directionalMasks[1] = new DirectionalMask('d');
    _directionalMasks[2] = new DirectionalMask('n');
    _directionalMasks[3] = new DirectionalMask('s');
    _directionalMasks[4] = new DirectionalMask('e');
    _directionalMasks[5] = new DirectionalMask('w');
}

Thinning6Iterations::~Thinning6Iterations() {
    delete _directionalMasks[0];
    delete _directionalMasks[1];
    delete _directionalMasks[2];
    delete _directionalMasks[3];
    delete _directionalMasks[4];
    delete _directionalMasks[5];
}

bool Thinning6Iterations::matches(const size_t &direction, const int8_t *subVolume)
{
    if (_directionalMasks[direction]->matches(subVolume))
        return true;

    return false;
}

}
