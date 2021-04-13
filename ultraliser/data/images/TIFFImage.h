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

#ifndef ULTRALISER_DATA_IMAGE_TIFF_IMAGE_H
#define ULTRALISER_DATA_IMAGE_TIFF_IMAGE_H

#ifdef ULTRALISER_USE_TIFF
#include <tiffio.h>
#endif
#include <common/Common.h>

/**
 * @brief The TiffImage class
 */
class TiffImage
{
public:

    /**
     * @brief TiffImage
     */
    TiffImage();

    /**
     * @brief isPixelFilled
     * @param x
     * @param y
     * @return
     */
    bool isPixelFilled(int32_t x, int32_t y);

    /**
     * @brief setimageFile
     * @param imageFile
     */
    void setimageFile(const std::string &_imageFile);

    /**
     * @brief getimageFile
     * @return
     */
    std::string getimageFile();

    /**
     * @brief setCompressionOutput
     * @param compressionLevel
     */
    void setCompressionOutput(const uint32_t &compressionLevel);

    /**
     * @brief getCompressionOuput
     * @return
     */
    uint32_t getCompressionOuput();

    /**
     * @brief readImage
     */
    void readImage();

    /**
     * @brief writeImage
     * @param outFile
     * @return
     */
    bool writeImage(const std::string &outFile);

    /**
     * @brief getComplementaryColour
     * Provide some basic Functions that all image formats may find useful.
     * @param rgb
     * @return
     */
    int getComplementaryColour(const std::vector<int32_t> &rgb);

    /**
     * @brief rgbToHex
     * Converts rgb(a) values to hex values
     * @param rgb
     * @return
     */
    std::string rgbToHex(const std::vector<int32_t> &rgb);

    /**
     * @brief hextToRGB
     * Converts hex to rgb(a) values
     * @param hexval
     * @return
     */
    std::vector<int32_t> hextToRGB(const std::string &hexVal);


    /**
     * @brief subtractHex
     * Subtracs two hex values. in this case it subtractes the given hex value
     * from white (FFFFFF) which results in the complementary color.
     * @param hexval
     * @return
     */
    std::string subtractHex(const std::string &hexValue);

    /**
     * @brief hexToDec
     * Transform hex to decimal
     * @param hexVal
     * @return
     */
    int hexToDec(const std::string &hexVal);

private:

    /**
     * @brief _compressionOutput
     */
    uint32_t _compressionOutput;

    /**
     * @brief _imageBuffer
     */
    std::vector< uint32 > _imageBuffer;

    /**
     * @brief _imageWidth
     */
    int32_t _imageWidth;

    /**
     * @brief _imageHeight
     */
    int32_t _imageHeight;

    /**
     * @brief _imageOrientation
     */
    uint32 _imageOrientation;

    /**
     * @brief imageFile
     */
    std::string _imageFile;

    /**
     * @brief _hexadecimal
     */
    static std::string _hexadecimal;
};


#endif // ULTRALISER_DATA_IMAGE_TIFF_IMAGE_H
