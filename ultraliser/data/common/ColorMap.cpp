/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

#include "ColorMap.h"
#include "ColorMaps.hh"
#include <common/Common.h>
#include <common/Logging.h>

namespace Ultraliser
{
namespace ColorMap
{

Vector3f getRGBColorF(const float& value,
                      const float& minValue,
                      const float& maxValue,
                      const std::string &colorMap)
{
    // The color map length, it will be always 256
    const uint64_t colorMapLength = 256;

    // Get the index of the color map
    const uint64_t colorMapIndex = F2UI64((colorMapLength - 1) *
                                          (value - minValue) / (maxValue - minValue));

    // Return the corresponding color from the selected color map
    if (colorMap == "VIRIDIS")
        return VIRIDIS_CM[colorMapIndex];

    else if (colorMap == "PLASMA")
        return PLASMA_CM[colorMapIndex];

    else if (colorMap == "INFERNO")
        return INFERNO_CM[colorMapIndex];

    else if (colorMap == "MAGMA")
        return MAGMA_CM[colorMapIndex];

    else if (colorMap == "CIVIDIS")
        return CIVIDIS_CM[colorMapIndex];

    else if (colorMap == "GREYS")
        return GREYS_CM[colorMapIndex];

    else if (colorMap == "PURPLES")
        return PURPLES_CM[colorMapIndex];

    else if (colorMap == "BLUES")
        return BLUES_CM[colorMapIndex];

    else if (colorMap == "GREENS")
        return GREENS_CM[colorMapIndex];

    else if (colorMap == "ORANGES")
        return ORANGES_CM[colorMapIndex];

    else if (colorMap == "REDS")
        return REDS_CM[colorMapIndex];

    else if (colorMap == "YLORBR")
        return YLORBR_CM[colorMapIndex];

    else if (colorMap == "YLORRD")
        return YLORRD_CM[colorMapIndex];

    else if (colorMap == "ORRD")
        return ORRD_CM[colorMapIndex];

    else if (colorMap == "PURD")
        return PURD_CM[colorMapIndex];

    else if (colorMap == "RDPU")
        return RDPU_CM[colorMapIndex];

    else if (colorMap == "BUPU")
        return BUPU_CM[colorMapIndex];

    else if (colorMap == "GNBU")
        return GNBU_CM[colorMapIndex];

    else if (colorMap == "PUBU")
        return PUBU_CM[colorMapIndex];

    else if (colorMap == "YLGNBU")
        return YLGNBU_CM[colorMapIndex];

    else if (colorMap == "PUBUGN")
        return PUBUGN_CM[colorMapIndex];

    else if (colorMap == "BUGN")
        return BUGN_CM[colorMapIndex];

    else if (colorMap == "YLGN")
        return YLGN_CM[colorMapIndex];

    else if (colorMap == "BINARY")
        return BINARY_CM[colorMapIndex];

    else if (colorMap == "GIST_YARG")
        return GIST_YARG_CM[colorMapIndex];

    else if (colorMap == "GIST_GRAY")
        return GIST_GRAY_CM[colorMapIndex];

    else if (colorMap == "NORMAL_GRAY")
        return GRAY_CM[colorMapIndex];

    else if (colorMap == "BONE")
        return BONE_CM[colorMapIndex];

    else if (colorMap == "PINK")
        return PINK_CM[colorMapIndex];

    else if (colorMap == "SPRING")
        return SPRING_CM[colorMapIndex];

    else if (colorMap == "SUMMER")
        return SUMMER_CM[colorMapIndex];

    else if (colorMap == "AUTUMN")
        return AUTUMN_CM[colorMapIndex];

    else if (colorMap == "WINTER")
        return WINTER_CM[colorMapIndex];

    else if (colorMap == "COOL")
        return COOL_CM[colorMapIndex];

    else if (colorMap == "WISTIA")
        return WISTIA_CM[colorMapIndex];

    else if (colorMap == "HOT")
        return HOT_CM[colorMapIndex];

    else if (colorMap == "AFMHOT")
        return AFMHOT_CM[colorMapIndex];

    else if (colorMap == "GIST_HEAT")
        return GIST_HEAT_CM[colorMapIndex];

    else if (colorMap == "COPPER")
        return COPPER_CM[colorMapIndex];

    else if (colorMap == "PIYG")
        return PIYG_CM[colorMapIndex];

    else if (colorMap == "PRGN")
        return PRGN_CM[colorMapIndex];

    else if (colorMap == "BRBG")
        return BRBG_CM[colorMapIndex];

    else if (colorMap == "PUOR")
        return PUOR_CM[colorMapIndex];

    else if (colorMap == "RDGY")
        return RDGY_CM[colorMapIndex];

    else if (colorMap == "RDBU")
        return RDBU_CM[colorMapIndex];

    else if (colorMap == "RDYLBU")
        return RDYLBU_CM[colorMapIndex];

    else if (colorMap == "RDYLGN")
        return RDYLGN_CM[colorMapIndex];

    else if (colorMap == "SPECTRAL")
        return SPECTRAL_CM[colorMapIndex];

    else if (colorMap == "COOLWARM")
        return COOLWARM_CM[colorMapIndex];

    else if (colorMap == "BWR")
        return BWR_CM[colorMapIndex];

    else if (colorMap == "SEISMIC")
        return SEISMIC_CM[colorMapIndex];

    else if (colorMap == "TWILIGHT")
        return TWILIGHT_CM[colorMapIndex];

    else if (colorMap == "TWILIGHT_SHIFTED")
        return TWILIGHT_SHIFTED_CM[colorMapIndex];

    else if (colorMap == "HSV")
        return HSV_CM[colorMapIndex];

    else if (colorMap == "PASTEL1")
        return PASTEL1_CM[colorMapIndex];

    else if (colorMap == "PASTEL2")
        return PASTEL2_CM[colorMapIndex];

    else if (colorMap == "PAIRED")
        return PAIRED_CM[colorMapIndex];

    else if (colorMap == "ACCENT")
        return ACCENT_CM[colorMapIndex];

    else if (colorMap == "DARK2")
        return DARK2_CM[colorMapIndex];

    else if (colorMap == "SET1")
        return SET1_CM[colorMapIndex];

    else if (colorMap == "SET2")
        return SET2_CM[colorMapIndex];

    else if (colorMap == "SET3")
        return SET3_CM[colorMapIndex];

    else if (colorMap == "TAB10")
        return TAB10_CM[colorMapIndex];

    else if (colorMap == "TAB20")
        return TAB20_CM[colorMapIndex];

    else if (colorMap == "TAB20B")
        return TAB20B_CM[colorMapIndex];

    else if (colorMap == "TAB20C")
        return TAB20C_CM[colorMapIndex];

    else if (colorMap == "FLAG")
        return FLAG_CM[colorMapIndex];

    else if (colorMap == "PRISM")
        return PRISM_CM[colorMapIndex];

    else if (colorMap == "OCEAN")
        return OCEAN_CM[colorMapIndex];

    else if (colorMap == "GIST_EARTH")
        return GIST_EARTH_CM[colorMapIndex];

    else if (colorMap == "TERRAIN")
        return TERRAIN_CM[colorMapIndex];

    else if (colorMap == "GIST_STERN")
        return GIST_STERN_CM[colorMapIndex];

    else if (colorMap == "GNUPLOT")
        return GNUPLOT_CM[colorMapIndex];

    else if (colorMap == "GNUPLOT2")
        return GNUPLOT2_CM[colorMapIndex];

    else if (colorMap == "CMRMAP")
        return CMRMAP_CM[colorMapIndex];

    else if (colorMap == "CUBEHELIX")
        return CUBEHELIX_CM[colorMapIndex];

    else if (colorMap == "BRG")
        return BRG_CM[colorMapIndex];

    else if (colorMap == "GIST_RAINBOW")
        return GIST_RAINBOW_CM[colorMapIndex];

    else if (colorMap == "RAINBOW")
        return RAINBOW_CM[colorMapIndex];

    else if (colorMap == "JET")
        return JET_CM[colorMapIndex];

    else if (colorMap == "NIPY_SPECTRAL")
        return NIPY_SPECTRAL_CM[colorMapIndex];

    else if (colorMap == "GIST_NCAR")
        return GIST_NCAR_CM[colorMapIndex];

    else
        return HSV_CM[colorMapIndex];
}

Vec3i_32 getRGBColorI(const float& index,
                      const float& minValue,
                      const float& maxValue,
                      const std::string &colorMap)
{
    // Get the color map float value
    Vector3f floatRGB = getRGBColorF(index, minValue, maxValue, colorMap);

    // Convert it to integer
    Vec3i_32 intRGB;
    intRGB[0] = F2I32(floatRGB[0] * 256.f);
    intRGB[1] = F2I32(floatRGB[1] * 256.f);
    intRGB[2] = F2I32(floatRGB[2] * 256.f);

    // Return the RGB color map
    return intRGB;
}

}
}
