#pragma once

#include <ctype.h>
namespace Ultraliser
{

// these elements to access the neighbours with the shifts
static const int VDX[26] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1,-1, 0, 1};
static const int VDY[26] = {-1,-1,-1,-1,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
static const int VDZ[26] = { 1, 1, 1, 0, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0,-1,-1,-1, 1, 1, 1, 0, 0, 0,-1,-1,-1};

}
