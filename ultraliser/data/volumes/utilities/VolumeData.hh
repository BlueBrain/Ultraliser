#pragma once

#include <common/Headers.hh>
#include <geometry/Geometry.h>
#include <data/volumes/utilities/VolumeType.hh>

#ifdef ULTRALISER_USE_NRRD
#include <nrrdloader/DataMangler.h>
#endif

namespace Ultraliser
{

/**
 * @brief The VolumeData struct
 */
struct VolumeData
{
    VolumeData()
    {
        width = 0;
        height = 0;
        depth = 0;
    }

    /**
     * @brief width
     * Volume width
     */
    size_t width;

    /**
     * @brief height
     * Volume Height
     */
    size_t height;

    /**
     * @brief depth
     * Volume depth
     */
    size_t depth;

    /**
     * @brief type
     * The type of the volume
     */
    VOLUME_TYPE type;
};

#ifdef ULTRALISER_USE_NRRD
/**
 * @brief The NRRDVolumeData struct
 */
struct NRRDVolumeData : VolumeData
{
    NRRDVolumeData()
    {
        center = Vector3f(0.f);
        scale = Vector3f(1.f);
    }

    /**
     * @brief center
     * Volume center
     */
    Vector3f center;

    /**
     * @brief scale
     * Volume scale
     */
    Vector3f scale;

    /**
     * @brief data
     */
    libNRRD::IDataMangler* data;
};
#endif

/**
 * @brief The VolumeData struct
 */
struct UltraliserVolumeData : VolumeData
{
    /**
     * @brief data
     * Volume data stored in a string to be parsed and converted to whatever data type later.
     */
    std::string data;
};

}
