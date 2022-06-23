#pragma once

namespace Ultraliser
{

/**
 * @brief The VOLUME_TYPE enum
 */
enum VOLUME_TYPE
{
    // A bit volume, where each voxel is stored in a single bit (class 1)
    BIT,

    // Each voxel is stored in a single byte (8 bits) as uint8_t (class 2)
    UI8,

    // Each voxel is stored in a single word (16 bits) as uint16_t (class 2)
    UI16,

    // Each voxel is stored in a doubleword (32 bits) as uint32_t (class 2)
    UI32,

    // Each voxel is stored in a quad-word (64 bits) as uint64_t (class 2)
    UI64,

    // Each voxel is stored as a single-precision floating-point or float value (class 3)
    F32,

    // Each voxel is stored as a double-precision floating-point or double value (class 3)
    F64
};

}
