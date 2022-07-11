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


#include "LSFM.h"
#include <common/Common.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{
namespace PBRT
{

/**
 * @brief createConfigLSFMXY
 * @param configFile
 * @param params
 * @return
 */
std::string createConfigLSFMXY(const std::string configFile,
                               const PBRTParameters params)
{
    // Parse the input configuration into array
    std::vector< std::string > inputConfig =
            File::parsePBRTConfig(configFile);

    // The output configuration has the same size of the input config
    std::string outputConfig("");

    // Update the symbols
    for (size_t i = 0; i < inputConfig.size(); i++)
    {
        std::string line = inputConfig[i];

        // Update the stream and then replace a substring
        std::stringstream stream;

        /// Frustum
        stream << 0.0;
        String::replaceSubstring(line, "X_CAMERA", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.z();
        String::replaceSubstring(line, "Z_CAMERA", stream.str());
        CLEAR_STREAM;

        /// Camera
        stream << params.pMin.x();
        String::replaceSubstring(line, "X_MIN_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.x();
        String::replaceSubstring(line, "X_MAX_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.y();
        String::replaceSubstring(line, "Y_MIN_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.y();
        String::replaceSubstring(line, "Y_MAX_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.focalDistance;
        String::replaceSubstring(line, "FOCAL_DISTANCE", stream.str());
        CLEAR_STREAM;

        /// Film
        Resolution resolution = computeFilmResolution
                (PROJECTION::XY, params.resolution, params.pMin, params.pMax);

        stream << resolution.x;
        String::replaceSubstring(line, "X_RESOLUTION", stream.str());
        CLEAR_STREAM;

        stream << resolution.y;
        String::replaceSubstring(line, "Y_RESOLUTION", stream.str());
        CLEAR_STREAM;

        stream << resolution.x;
        String::replaceSubstring(line, "X0", stream.str());
        CLEAR_STREAM;

        stream << resolution.y;
        String::replaceSubstring(line, "Y0", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 2;
        String::replaceSubstring(line, "X1", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 2;
        String::replaceSubstring(line, "Y1", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 4;
        String::replaceSubstring(line, "X2", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 4;
        String::replaceSubstring(line, "Y2", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 8;
        String::replaceSubstring(line, "X3", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 8;
        String::replaceSubstring(line, "Y3", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 16;
        String::replaceSubstring(line, "X4", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 16;
        String::replaceSubstring(line, "Y4", stream.str());
        CLEAR_STREAM;

        stream << params.filmPrefix << ".exr";
        String::replaceSubstring(line, "IMAGE_NAME", stream.str());
        CLEAR_STREAM;

        /// Volume Integrator
        // Compute the step size
        stream << params.step ;
        String::replaceSubstring(line, "STEP_SIZE", stream.str());
        CLEAR_STREAM;

        /// LSFM Illumination
        stream << params.pMin.x() - 1.0;
        String::replaceSubstring(line, "X_LIGHT_OFFSET", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.z() - params.focalDistance;
        String::replaceSubstring(line, "Z_LIGHT_OFFSET", stream.str());
        CLEAR_STREAM;

        stream << 90.0;
        String::replaceSubstring(line, "Y_LIGHT_ROTATION", stream.str());
        CLEAR_STREAM;

        stream << (params.pMax.y() - params.pMin.y()) * 2;
        String::replaceSubstring(line, "Y_LIGHT_SIZE", stream.str());
        CLEAR_STREAM;

        /// Volume Grid
        stream << params.pMin.x();
        String::replaceSubstring(line, "P0_X", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.y();
        String::replaceSubstring(line, "P0_Y", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.z();
        String::replaceSubstring(line, "P0_Z", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.x();
        String::replaceSubstring(line, "P1_X", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.y();
        String::replaceSubstring(line, "P1_Y", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.z();
        String::replaceSubstring(line, "P1_Z", stream.str());
        CLEAR_STREAM;

        stream << params.volumePrefix;
        String::replaceSubstring(line, "VOLUME_PREFIX", stream.str());
        CLEAR_STREAM;

        outputConfig += line + "\n";
    }

    // Return the final configuration
    return outputConfig;
}

/**
 * @brief createConfigLSFMZY
 * @param configFile
 * @param params
 * @return
 */
std::string createConfigLSFMZY(const std::string configFile,
                               const PBRTParameters params)
{
    // Parse the input config into array
    std::vector< std::string > inputConfig =
            File::parsePBRTConfig(configFile);

    // The output config has the same size of the input config
    std::string outputConfig("");

    // Update the symbols
    for (size_t i = 0; i < inputConfig.size(); i++)
    {
        std::string line = inputConfig[i];

        // Update the stream and then replace a substring
        std::stringstream stream;

        /// Frustum
        stream << params.pMax.x();
        String::replaceSubstring(line, "X_CAMERA", stream.str());
        CLEAR_STREAM;

        stream << 0.0;
        String::replaceSubstring(line, "Z_CAMERA", stream.str());
        CLEAR_STREAM;

        /// Camera
        stream << params.pMin.z();
        String::replaceSubstring(line, "X_MIN_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.z();
        String::replaceSubstring(line, "X_MAX_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.y();
        String::replaceSubstring(line, "Y_MIN_SCREEN", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.y();
        String::replaceSubstring(line, "Y_MAX_SCREEN", stream.str());
        CLEAR_STREAM;
        stream << params.focalDistance;
        String::replaceSubstring(line, "FOCAL_DISTANCE", stream.str());
        CLEAR_STREAM;

        /// Film
        Resolution resolution = computeFilmResolution
                (PROJECTION::ZY, params.resolution, params.pMin, params.pMax);

        stream << resolution.x;
        String::replaceSubstring(line, "X_RESOLUTION", stream.str());
        CLEAR_STREAM;

        stream << resolution.y;
        String::replaceSubstring(line, "Y_RESOLUTION", stream.str());
        CLEAR_STREAM;

        stream << resolution.x;
        String::replaceSubstring(line, "X0", stream.str());
        CLEAR_STREAM;

        stream << resolution.y;
        String::replaceSubstring(line, "Y0", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 2;
        String::replaceSubstring(line, "X1", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 2;
        String::replaceSubstring(line, "Y1", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 4;
        String::replaceSubstring(line, "X2", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 4;
        String::replaceSubstring(line, "Y2", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 8;
        String::replaceSubstring(line, "X3", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 8;
        String::replaceSubstring(line, "Y3", stream.str());
        CLEAR_STREAM;

        stream << resolution.x * 16;
        String::replaceSubstring(line, "X4", stream.str());
        CLEAR_STREAM;

        stream << resolution.y * 16;
        String::replaceSubstring(line, "Y4", stream.str());
        CLEAR_STREAM;

        stream << params.filmPrefix << ".exr";
        String::replaceSubstring(line, "IMAGE_NAME", stream.str());
        CLEAR_STREAM;

        /// Volume Integrator
        stream << params.step;
        String::replaceSubstring(line, "STEP_SIZE", stream.str());
        CLEAR_STREAM;

        /// LSFM Illumination
        stream << params.pMax.x() - params.focalDistance;;
        String::replaceSubstring(line, "X_LIGHT_OFFSET", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.z() - 1.0;
        String::replaceSubstring(line, "Z_LIGHT_OFFSET", stream.str());
        CLEAR_STREAM;

        stream << 0.0;
        String::replaceSubstring(line, "Y_LIGHT_ROTATION", stream.str());
        CLEAR_STREAM;

        stream << (params.pMax.y() - params.pMin.y()) * 2;
        String::replaceSubstring(line, "Y_LIGHT_SIZE", stream.str());
        CLEAR_STREAM;

        /// Volume Grid
        stream << params.pMin.x();
        String::replaceSubstring(line, "P0_X", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.y();
        String::replaceSubstring(line, "P0_Y", stream.str());
        CLEAR_STREAM;

        stream << params.pMin.z();
        String::replaceSubstring(line, "P0_Z", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.x();
        String::replaceSubstring(line, "P1_X", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.y();
        String::replaceSubstring(line, "P1_Y", stream.str());
        CLEAR_STREAM;

        stream << params.pMax.z();
        String::replaceSubstring(line, "P1_Z", stream.str());
        CLEAR_STREAM;

        stream << params.volumePrefix;
        String::replaceSubstring(line, "VOLUME_PREFIX", stream.str());
        CLEAR_STREAM;

        outputConfig += line + "\n";
    }

    // Return the final configuration
    return outputConfig;
}

}
}
