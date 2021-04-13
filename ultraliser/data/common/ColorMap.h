/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_DATA_COMMON_COLORMAP_H
#define ULTRALISER_DATA_COMMON_COLORMAP_H

#include <data/common/CommonData.h>
#include <data/common/ColorMaps.hh>

namespace Ultraliser
{
namespace ColorMap
{

/**
 * @brief getRGBColorF
 * Return an RGB color for a given
 * @param value
 * @param minValue
 * @param maxValue
 * @param colorMap
 * @return
 */
Vector3f getRGBColorF(const float &value,
                      const float &minValue,
                      const float &maxValue,
                      const std::string &colorMap="VIRIDIS");

/**
 * @brief getRGBColorI
 * @param value
 * @param minValue
 * @param maxValue
 * @param colorMap
 * @return
 */
Vec3i_32 getRGBColorI(const float& value,
                      const float& minValue,
                      const float& maxValue,
                      const std::string &colorMap="VIRIDIS");

/**
 * A list of all the supported color-maps in Ultraliser.
 * These color-maps were generated with a python-based utility using Matplotlib.
 */
const static std::vector<std::string> MAPS =
{
    "VIRIDIS",
    "PLASMA",
    "INFERNO",
    "MAGMA",
    "CIVIDIS",
    "GREYS",
    "PURPLES",
    "BLUES",
    "GREENS",
    "ORANGES",
    "REDS",
    "YLORBR",
    "YLORRD",
    "ORRD",
    "PURD",
    "RDPU",
    "BUPU",
    "GNBU",
    "PUBU",
    "YLGNBU",
    "PUBUGN",
    "BUGN",
    "YLGN",
    "BINARY",
    "GIST_YARG",
    "GIST_GRAY",
    "GRAY",
    "BONE",
    "PINK",
    "SPRING",
    "SUMMER",
    "AUTUMN",
    "WINTER",
    "COOL",
    "WISTIA",
    "HOT",
    "AFMHOT",
    "GIST_HEAT",
    "COPPER",
    "PIYG",
    "PRGN",
    "BRBG",
    "PUOR",
    "RDGY",
    "RDBU",
    "RDYLBU",
    "RDYLGN",
    "SPECTRAL",
    "COOLWARM",
    "BWR",
    "SEISMIC",
    "TWILIGHT",
    "TWILIGHT_SHIFTED",
    "HSV",
    "PASTEL1",
    "PASTEL2",
    "PAIRED",
    "ACCENT",
    "DARK2",
    "SET1",
    "SET2",
    "SET3",
    "TAB10",
    "TAB20",
    "TAB20B",
    "TAB20C",
    "FLAG",
    "PRISM",
    "OCEAN",
    "GIST_EARTH",
    "TERRAIN",
    "GIST_STERN",
    "GNUPLOT",
    "GNUPLOT2",
    "CMRMAP",
    "CUBEHELIX",
    "BRG",
    "GIST_RAINBOW",
    "RAINBOW",
    "JET",
    "NIPY_SPECTRAL",
    "GIST_NCAR"
};

}
}

#endif // ULTRALISER_DATA_COMMON_COLORMAP_H
