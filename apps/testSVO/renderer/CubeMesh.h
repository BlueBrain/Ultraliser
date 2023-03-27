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

#include <glm/glm.hpp>

#include <glad/glad.h>

#include <vector>

namespace svorender
{
class CubeMesh
{
public:
    CubeMesh(const glm::vec3 &min = glm::vec3(0.f), const glm::vec3 &max = glm::vec3(1.f));

    void bind() const;
    void unbind() const;

    void destroy();

private:
    void _enableAttributes();

private:
    GLuint _gpuMeshVAO;
    GLuint _buffers[2];
};
}