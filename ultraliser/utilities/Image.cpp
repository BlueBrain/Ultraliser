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

#include <utilities/Image.h>
#include <common/Common.h>
#include <utilities/TypeConversion.h>

namespace Ultraliser
{
namespace Utilities
{

void savePPMLuminanceImage(const std::string &imageName,
                           const uint8_t *imageData,
                           const int64_t &width,
                           const int64_t &height)
{
    // Make the file
    std::stringstream imageStream;
    imageStream << imageName << PPM_EXTENSION;
    FILE* image = fopen(imageStream.str().c_str(), "wb");

    // Add the header
    fprintf(image, "P6\n%ld %ld\n255\n", width, height);

    size_t index = 0;
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            uint8_t color[3];
            size_t index1D = I2UI64(width * height) - index;
            color[0] = imageData[index1D]; // R
            color[1] = imageData[index1D]; // G
            color[2] = imageData[index1D]; // B
            fwrite(color, 1, 3, image);
            index++;
        }
    }

    fclose(image);
}

void savePPMLuminanceImage(const std::string &imageName,
                           const uint16_t *imageData,
                           const int64_t &width,
                           const int64_t &height)
{
    // Make the file
    std::stringstream imageStream;
    imageStream << imageName << PPM_EXTENSION;
    FILE* image = fopen(imageStream.str().c_str(), "wb");

    // Add the header
    fprintf(image, "P6\n%ld %ld\n65535\n", width, height);

    size_t index = 0;
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            uint16_t color[3];
            size_t index1D = I2UI64(width * height) - index;
            color[0] = imageData[index1D]; // R
            color[1] = imageData[index1D]; // G
            color[2] = imageData[index1D]; // B
            fwrite(color, 2, 3, image);
            index++;
        }
    }

    fclose(image);
}

void savePPMColoredImage(const std::string &imageName,
                  const Vector3f* imageData,
                  const int64_t &width,
                  const int64_t &height)
{
    // Make the file
    std::stringstream imageStream;
    imageStream << imageName << PPM_EXTENSION;
    FILE* image = fopen(imageStream.str().c_str(), "wb");

    // Add the header
    fprintf(image, "P6\n%ld %ld\n255\n", width, height);

    size_t index = 0;
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            Vector3f normalizeColor = imageData[I2UI64(width * height) - index] * 255.0;

            uint8_t color[3];
            color[0] = F2UI8(normalizeColor.x()); // R
            color[1] = F2UI8(normalizeColor.y()); // G
            color[2] = F2UI8(normalizeColor.z()); // B

            fwrite(color, 1, 3, image);
            index++;
        }
    }

    fclose(image);
}

void savePPMColoredImage(const std::string &imageName,
                  const Vector4f* imageData,
                  const int64_t &width,
                  const int64_t &height)
{
    // Make the file
    std::stringstream imageStream;
    imageStream << imageName << PPM_EXTENSION;
    FILE* image = fopen(imageStream.str().c_str(), "wb");

    // Add the header
    fprintf(image, "P6\n%ld %ld\n255\n", width, height);

    size_t index = 0;
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            Vector4f normalizeColor = imageData[I2UI64(width * height) - index] * 255.0;

            uint8_t color[3];
            color[0] = static_cast<uint8_t> (normalizeColor.x()); // R
            color[1] = static_cast<uint8_t> (normalizeColor.y()); // G
            color[2] = static_cast<uint8_t> (normalizeColor.z()); // B

            fwrite(color, 1, 3, image);
            index++;
        }
    }

    fclose(image);
}

void saveBrainbowImage(const std::string &imageName,
                       const Vector4f *imageData,
                       const int64_t &width,
                       const int64_t &height)
{

    Vector4f* rgbImage = new Vector4f[width * height];

    size_t index = 0;
    for (int64_t j = 0; j < height; ++j)
    {
        for (int64_t i = 0; i < width; ++i)
        {
            Vector4f color = imageData[I2UI64(width * height)] * 255.0;

            rgbImage[3 * index].x() = color.x();
            rgbImage[3 * index].y() = color.y();
            rgbImage[3 * index].z() = color.z();
            index++;
        }
    }

    // Write the RGB image
    std::stringstream imageStream;
    imageStream << imageName << "_brainbow" << EXR_EXTENSION;
    savePPMColoredImage(imageStream.str().c_str(), rgbImage, width, height);

    // Release temporary image memory
    delete[] rgbImage;
}

}
}
