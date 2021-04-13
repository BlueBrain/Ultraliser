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

#ifndef ULTRALISER_UTILITIES_IMAGE_H
#define ULTRALISER_UTILITIES_IMAGE_H

#include <common/Common.h>
#include <math/Math.h>

namespace Ultraliser
{
namespace Utilities
{

/**
 * @brief savePPMLuminanceImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void savePPMLuminanceImage(const std::string &imageName,
                  const uint8_t *imageData,
                  const int64_t &width,
                  const int64_t &height);

void savePPMColoredImage(const std::string &imageName,
                         const Vector3f* imageData,
                         const int64_t &width,
                         const int64_t &height);

/**
 * @brief savePPMColoredImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void savePPMColoredImage(const std::string &imageName,
                         const Vector4f* imageData,
                         const int64_t &width,
                         const int64_t &height);

void saveEXRImage(const std::string &imageName,
                  float* imageData,
                  float* imageAlpha,
                  const int64_t &width, const int64_t &height,
                  const int32_t &x, const int32_t &y,
                  const int32_t &xOffset, const int32_t &yOffset);

/**
 * @brief saveEXRLuminanceImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void saveEXRLuminanceImage(const std::string &imageName,
                           const float* imageData,
                           const int64_t &width,
                           const int64_t &height);

/**
 * @brief saveEXRLuminanceImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void saveEXRLuminanceImage(const std::string &imageName,
                           const uint8_t* imageData,
                           const int64_t &width,
                           const int64_t &height);

/**
 * @brief saveEXRRGBImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void saveEXRRGBImage(const std::string &imageName,
                     const uint8_t* imageData,
                     const int64_t &width,
                     const int64_t &height);

/**
 * @brief saveBrainbowImage
 * @param imageName
 * @param imageData
 * @param width
 * @param height
 */
void saveBrainbowImage(const std::string &imageName,
                       const Vector4f *imageData,
                       const int64_t &width,
                       const int64_t &height);

}
}

#endif // ULTRALISER_UTILITIES_IMAGE_H