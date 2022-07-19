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

#include "TIFFImage.h"
#include <utilities/Utilities.h>

std::string TiffImage::_hexadecimal = "0123456789ABCDEF";

TiffImage::TiffImage() : _compressionOutput(5)
{
    /// EMPTY CONSTRUCTOR
}

void TiffImage::setimageFile(const std::string &imageFile)
{
    this->_imageFile = imageFile;
}

std::string TiffImage::getimageFile()
{
    return this->_imageFile;
}

void TiffImage::setCompressionOutput(const uint32_t &compressionLevel)
{
    // 1        No Compression
    // 5        LZW
    // 7        JPEG
    // 32909    DEFLATE
    std::vector<int> compLevels = {1, 5, 7, 32909};

    if (find(compLevels.begin(),
             compLevels.end(), compressionLevel) != compLevels.end())
    {
        this->_compressionOutput = compressionLevel;
    }
    else
    {
        LOG_ERROR( "Supported values for compression are 5(LZW), 7(JPEG), 32909(DEFLATE)!");
    }
}

uint32_t TiffImage::getCompressionOuput()
{
    return this->_compressionOutput;
}

uint32_t rgbaFloatToInt(const uint32_t& r, const uint32_t& g, const uint32_t& b, const uint32_t& a)
{

    return a << 24 | b << 16 | g << 8 | r;
}


bool TiffImage::isPixelFilled(const uint32_t &x, const uint32_t &y)
{
#ifdef ULTRALISER_USE_TIFF
    size_t index = x + (_imageWidth * y);
    uint32_t value = _imageBuffer[index];

    uint32_t r = TIFFGetR(value);
    uint32_t g = TIFFGetG(value);
    uint32_t b = TIFFGetB(value);
    uint32_t a = TIFFGetA(value);

    if (r > 0 || g > 0 || b > 0)
        return true;
    return false;
#else
    return false;
#endif
}

void TiffImage::readImage()
{
#ifdef ULTRALISER_USE_TIFF
    // Open the TIFF image
    TIFF *tiffImage = TIFFOpen(this->getimageFile().c_str(), "r");

    // If opened
    if (tiffImage)
    {
        // Read the data
        int32_t numberPixels;
        uint32_t* raster;

        TIFFGetField(tiffImage, TIFFTAG_IMAGEWIDTH, &_imageWidth);
        TIFFGetField(tiffImage, TIFFTAG_IMAGELENGTH, &_imageHeight);
        TIFFGetField(tiffImage, TIFFTAG_ORIENTATION, &_imageOrientation);

        numberPixels = _imageWidth * _imageHeight;
        raster = static_cast< uint32_t* > (_TIFFmalloc(numberPixels * sizeof (uint32_t)));

        if (raster != nullptr)
        {
            if (TIFFReadRGBAImage(tiffImage, (_imageWidth), (_imageHeight), raster, 0))
            {
                // Process raster data
                // Write data in vector member so we can do with it what we want
                uint32_t *rasterPrint = raster;
                for ( int32_t n = 0; n < numberPixels; ++n)
                {
                    uint32_t r, g, b, a;
                    r = TIFFGetR(*rasterPrint);
                    g = TIFFGetG(*rasterPrint);
                    b = TIFFGetB(*rasterPrint);
                    a = TIFFGetA(*rasterPrint);

                    uint32_t value;
                    if (r > 0 || g > 0 || b > 0)
                    {
                        value = rgbaFloatToInt(255,255,255,255);
                    }
                    else
                    {
                        value = rgbaFloatToInt(r,g,b,a);
                    }

                    _imageBuffer.push_back(value);
                    *rasterPrint++;
                }
            }
            _TIFFfree(raster);

            // TIFFRedRGBAImage is starting in the lower left corner, so we
            // got to swap our vector. this means we hve to swap first row with
            // last and so on.
            uint32_t upBufPos, downBufPos;
            for (uint32_t i = 0 ; i < I2UI32(_imageHeight) / 2; ++i)
            {
                for (uint32_t j = 0 ; j < _imageWidth; ++j)
                {
                    upBufPos = i * _imageWidth + j;
                    if (i * j == 0)
                    {
                        upBufPos = i+j;
                    }

                    downBufPos = (((this->_imageHeight) - i - 1) *
                                   (this->_imageWidth)) + j;
                    std::swap(this->_imageBuffer[upBufPos],
                              this->_imageBuffer[downBufPos]);
                }
            }
        }
        TIFFClose(tiffImage);
    }
#endif
}

bool TiffImage::writeImage(const std::string &outFile)
{
#ifdef ULTRALISER_USE_TIFF
    bool retBool = false;
    TIFF *outputImage;

    // Open the TIFF file
    if ((outputImage = TIFFOpen(outFile.c_str(), "w")) == nullptr)
    {
        LOG_ERROR("Unable to write tif file [ %s ]", outFile.c_str());
    }

    // We need to set some values for basic tags before we can add any data
    TIFFSetField(outputImage, TIFFTAG_IMAGEWIDTH, this->_imageWidth);
    TIFFSetField(outputImage, TIFFTAG_IMAGELENGTH, this->_imageHeight);

    // TODO: We should read these values from the input picture and use
    // them accordingly
    TIFFSetField(outputImage, TIFFTAG_BITSPERSAMPLE, 8);
    TIFFSetField(outputImage, TIFFTAG_SAMPLESPERPIXEL, 4);
    TIFFSetField(outputImage, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);

    TIFFSetField(outputImage, TIFFTAG_COMPRESSION, this->getCompressionOuput());
    TIFFSetField(outputImage, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField(outputImage, TIFFTAG_ORIENTATION, this->_imageOrientation);

    // Write the information to the file
    if (TIFFWriteEncodedStrip(outputImage, 0,
                              &_imageBuffer[0], _imageWidth*_imageHeight * 4) > -1)
    {
        // If bigger than -1 writing was successful
        retBool= true;
    }


    // Close the file
    TIFFClose(outputImage);
    return retBool;
#else
    return false;
#endif
}

int TiffImage::getComplementaryColour(const std::vector< int32_t > &rgb)
{
    // RGB to HEX
    std::string hexValue = this->rgbToHex(rgb);

    // Get hex value of complementary color by subtracting from white
    std::string subtractedHex = this->subtractHex(hexValue);

    //Transform back to RGB
    // complementaryRGB = this->hextToRGB(subtractedHex);
    return this->hexToDec(subtractedHex);
}


std::string TiffImage::rgbToHex(const std::vector< int32_t > &rgb)
{
    // Rgb to hex
    std::string hexvalue = "";
    std::stringstream stringStream;
    for (auto const &i : rgb)
    {
        // Use std::hex and streams to transform an int to hexvalue
        stringStream << std::hex << i;
        hexvalue.append(stringStream.str());
        if (stringStream.str().compare("0") == 0)
        {
            hexvalue.append("0");
        }

        stringStream.str("");
    }

    // hexValue to upper
    transform(hexvalue.begin(), hexvalue.end(),hexvalue.begin(), ::toupper);
    return hexvalue;
}
std::vector<int32_t> TiffImage::hextToRGB(const std::string &hexVal)
{
    std::vector<int> complementaryRGB;
    std::stringstream stringStream;

    int j = 0;
    for (size_t i = 0 ; i < hexVal.length(); i = i + 2)
    {
        std::string rgbString;
        rgbString += hexVal[i];
        rgbString += hexVal[i + 1];

        stringStream << std::hex << rgbString;
        stringStream >> j;
        complementaryRGB.push_back(j);
    }
    return complementaryRGB;
}

std::string TiffImage::subtractHex(const std::string &hexValue)
{
    std::string hexWhite = "FFFFFF";
    std::string subtractedHex;
    std::stringstream stringStream;

    for (size_t i = hexValue.length() - 1; i >= 0 ; --i)
    {
        // Get hex value when we subtract two other hex values
        size_t posHexWhite = this->_hexadecimal.find(hexWhite[I2UI64(i)]);
        size_t posHexval = this->_hexadecimal.find(hexValue[I2UI64(i)]);
        size_t newPos = posHexWhite - posHexval;

        // Use a std::stringstream otherwise there are strange side effects
        stringStream << _hexadecimal.at(newPos);

        // We start to subtract at the rightmost position but the hex has to
        // be written from left to right
        subtractedHex.insert(0, stringStream.str());

        // Empty std::stringstream
        stringStream.str("");
    }

    return subtractedHex;
}

int TiffImage::hexToDec(const std::string &hexval)
{
    int retDec = 0;
    std::stringstream stringStream;
    stringStream << hexval;
    stringStream >> std::hex >> retDec;
    return retDec;
}
