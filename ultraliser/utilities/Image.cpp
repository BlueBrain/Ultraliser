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

#include <utilities/Image.h>
#include <common/Common.h>
#include <utilities/TypeConversion.h>

using namespace Imf;
using namespace Imath;

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

    uint64_t index = 0;
    for (int64_t i = 0; i < width; ++i)
    {
        for (int64_t j = 0; j < height; ++j)
        {
            uint8_t color[3];
            uint64_t index1D = I2UI64(width * height) - index;
            color[0] = imageData[index1D]; // R
            color[1] = imageData[index1D]; // G
            color[2] = imageData[index1D]; // B
            fwrite(color, 1, 3, image);
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

    uint64_t index = 0;
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

    uint64_t index = 0;
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

static void writeEXRDataToFile(const std::string &imageName,
                               float *imagePixels,
                               float* imageAlpha,
                               const int64_t &width, const int64_t &height,
                               const int32_t& x, const int32_t& y,
                               const int32_t& xOffset, const int32_t& yOffset)
{
    // EXR RGBA array
    Rgba* rgbaArray = new Rgba[width * height];

    // Fill th e array
    for (int64_t i = 0; i < width * height; ++i)
        rgbaArray[i] = Rgba(imagePixels[3 * i],
                            imagePixels[3 * i + 1],
                            imagePixels[3 * i + 2],
                            imageAlpha ? imageAlpha[i]: 1.f);

    // Minor
    Box2i displayWindow(V2i(0, 0), V2i(x - 1, y - 1));
    Box2i dataWindow(V2i(xOffset, yOffset),
                     V2i(I2I32(xOffset + width) - 1,
                         I2I32(yOffset + height) - 1));

    try
    {
        RgbaOutputFile file(imageName.c_str(), displayWindow, dataWindow,
                            WRITE_RGBA);
        file.setFrameBuffer(rgbaArray - xOffset - yOffset * width, 1,
                            I2UI64(width));
        file.writePixels(I2I32(height));
    }
    catch(const std::exception&)
    {
        LOG_WARNING("Unable to write an EXR image file");
    }

    delete[] rgbaArray;
}

void saveEXRImage(const std::string &imageName,
                  float* imageData,
                  float* imageAlpha,
                  const int64_t &width, const int64_t &height,
                  const int32_t& x, const int32_t& y,
                  const int32_t& xOffset, const int32_t& yOffset)
{
    if (imageName.size() >= 5)
    {
        int64_t suffixOffset = I2I64(imageName.size() - 4);

        if (!strcmp(imageName.c_str() + suffixOffset, ".exr") ||
            !strcmp(imageName.c_str() + suffixOffset, ".EXR"))
            writeEXRDataToFile(imageName, imageData, imageAlpha, width, height,
                               x, y, xOffset, yOffset);
    }
}


void saveEXRLuminanceImage(const std::string &imageName,
                           const uint8_t* inputImage,
                           const int64_t &width,
                           const int64_t &height)
{
    // Convert the gray scale float image into an RGB one.
    const uint64_t numberPixels = I2UI64(width * height);
    float *rgbImage = new float[3 * numberPixels];

    for (int64_t i = 0; i < width; i++)
    {
        for (int64_t j = 0; j < height; j++)
        {
            uint64_t luminanceIndex = I2UI64(i + (j * width));
            uint64_t rgbIndex = I2UI64((width - 1 - i) + (height - 1 - j) * width);

            rgbImage[3 * rgbIndex + 0] = I2F(inputImage[luminanceIndex]);
            rgbImage[3 * rgbIndex + 1] = I2F(inputImage[luminanceIndex]);
            rgbImage[3 * rgbIndex + 2] = I2F(inputImage[luminanceIndex]);
        }
    }

    // Create the file
    std::stringstream stream;
    stream << imageName << EXR_EXTENSION;

    // Write RGB image
    saveEXRImage(stream.str().c_str(), rgbImage, nullptr, width, height,
                 I2I32(width), I2I32(height), 0, 0);

    // Release temporary image memory
    delete[] rgbImage;
}

void saveEXRRGBImage(const std::string fileName, const uint8_t* inputImage,
                     const int64_t &width, const int64_t &height)
{
    // Convert the gray scale float image into an RGB one.
    const int64_t numberPixels = width * height;
    float *rgbImage = new float[3 * numberPixels];

    int offset = 0;
    for (int64_t i = 0; i < width; i++)
    {
        for (int64_t j = 0; j < height; j++)
        {
            uint64_t luminanceIndex = I2UI64(i + (j  * width));
            uint64_t rgbIndex = I2UI64(i + (height - 1 - j) * width);

            rgbImage[3 * rgbIndex + 0] =
                    I2F(inputImage[3 * luminanceIndex + 0]);
            rgbImage[3 * rgbIndex + 1] =
                    I2F(inputImage[3 * luminanceIndex + 1]);
            rgbImage[3 * rgbIndex + 2] =
                    I2F(inputImage[3 * luminanceIndex + 2]);

            ++offset;
        }
    }

    // Write RGB image
    std::stringstream imageStream;
    imageStream << fileName << EXR_EXTENSION;
    saveEXRImage(imageStream.str().c_str(), rgbImage, nullptr, width, height,
                 I2I32(width), I2I32(height), 0, 0);

    // Release temporary image memory
    delete[] rgbImage;
}

void saveBrainbowImage(const std::string &imageName,
                       const Vector4f *imageData,
                       const int64_t &width,
                       const int64_t &height)
{

    float* rgbImage = new float[width * height * 3];

    uint64_t index = 0;
    for (int64_t j = 0; j < height; ++j)
    {
        for (int64_t i = 0; i < width; ++i)
        {
            Vector4f color = imageData[I2UI64(width * height) - index] * 255.0;

            rgbImage[3 * index + 0] = color.x();
            rgbImage[3 * index + 1] = color.y();
            rgbImage[3 * index + 2] = color.z();
            index++;
        }
    }

    // Write the RGB image
    std::stringstream imageStream;
    imageStream << imageName << "_brainbow" << EXR_EXTENSION;
    saveEXRImage(imageStream.str().c_str(), rgbImage, nullptr, width, height,
                 I2I32(width), I2I32(height), 0, 0);

    // Release temporary image memory
    delete[] rgbImage;
}




}
}
