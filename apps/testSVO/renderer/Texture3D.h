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

#pragma once

#include <glad/glad.h>

#include <glm/glm.hpp>

#include <cstdint>
#include <vector>

namespace svorender
{
class Texture3D
{
public:
    Texture3D(uint32_t dimX, uint32_t dimY, uint32_t dimZ, const std::vector<uint8_t> &volumeData);

    void bind() const;
    void unbind() const;

    void makeHandleResident() const;
    void makeHandleNonResident() const;

    GLuint64 getBindlessHandle() const;

    void destroy();

private:
    uint32_t _gridX;
    uint32_t _gridY;
    uint32_t _gridZ;

    GLuint _oglTexture;
    GLuint64 _oglTextureHandle;
};
}