/* Copyright (c) 2020, EPFL/Blue Brain Project
 * All rights reserved. Do not distribute without permission.
 * Responsible Author: Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
 *
 * This file is part of SimCrusher
 * <LINK>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
 * details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include "Texture3D.h"

#include <stdexcept>

namespace svorender
{
Texture3D::Texture3D(uint32_t dimX, uint32_t dimY, uint32_t dimZ, const std::vector<uint8_t> &volumeData)
    : _gridX(dimX)
    , _gridY(dimY)
    , _gridZ(dimZ)
    , _oglTexture(GL_INVALID_VALUE)
{
    if (volumeData.size() < _gridX * _gridY * _gridZ)
    {
        throw std::runtime_error("Texture3D: Given volume data is smaller than grid measures");
    }

    glGenTextures(1, &_oglTexture);
    if (_oglTexture == GL_INVALID_VALUE)
    {
        throw std::runtime_error("Texture3D: Could not create texture");
    }

    bind();

    glTexStorage3D(GL_TEXTURE_3D, 1, GL_R8, _gridX, _gridY, _gridZ);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    glTexSubImage3D(GL_TEXTURE_3D, 0, 0, 0, 0, _gridX, _gridY, _gridZ, GL_RED, GL_UNSIGNED_BYTE, volumeData.data());

    // glGenerateMipmap(GL_TEXTURE_3D);

    // Generate bindless handle
    _oglTextureHandle = glGetTextureHandleARB(_oglTexture);

    if (_oglTextureHandle == GL_INVALID_VALUE)
        throw std::runtime_error("Texture3D: Could not get bindless handle for texture");
}

void Texture3D::bind() const
{
    glBindTexture(GL_TEXTURE_3D, _oglTexture);
}

void Texture3D::unbind() const
{
    glBindTexture(GL_TEXTURE_3D, 0);
}

void Texture3D::makeHandleResident() const
{
    glMakeTextureHandleResidentARB(_oglTextureHandle);
}

void Texture3D::makeHandleNonResident() const
{
    glMakeTextureHandleNonResidentARB(_oglTextureHandle);
}

GLuint64 Texture3D::getBindlessHandle() const
{
    return _oglTextureHandle;
}

void Texture3D::destroy()
{
    makeHandleNonResident();
    glDeleteTextures(1, &_oglTexture);
}
}