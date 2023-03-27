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

#include "CubeMesh.h"

#include <array>
#include <cstring>
#include <stdexcept>

namespace
{
class VertexBufferFactory
{
public:
    static std::vector<float> generate(const glm::vec3 &min, const glm::vec3 &max)
    {
        const std::array<glm::vec3, 8> v = {
            glm::vec3(min.x, min.y, max.z),
            glm::vec3(max.x, min.y, max.z),
            min,
            glm::vec3(max.x, min.y, min.z),
            glm::vec3(min.x, max.y, max.z),
            max,
            glm::vec3(min.x, max.y, min.z),
            glm::vec3(max.x, max.y, min.z)};

        auto vertexBuffer = std::vector<float>(24);
        for (int i = 0; i < 8; ++i)
        {
            memcpy(&vertexBuffer[0] + 3 * i, &(v[i][0]), sizeof(float) * 3);
        }

        return vertexBuffer;
    }
};

class IndexBufferFactory
{
public:
    static std::vector<uint32_t> generate()
    {
        return std::vector<uint32_t>{0, 2, 1, 1, 2, 3, 5, 7, 4, 4, 7, 6, 3, 5, 1, 5, 3, 7,
                                     4, 6, 0, 0, 6, 2, 4, 0, 5, 5, 0, 1, 2, 6, 3, 3, 6, 7};
    }
};
}

namespace svorender
{
CubeMesh::CubeMesh(const glm::vec3 &min, const glm::vec3 &max)
    : _gpuMeshVAO(GL_INVALID_VALUE)
    , _buffers{GL_INVALID_VALUE, GL_INVALID_VALUE}
{
    auto vertices = VertexBufferFactory::generate(min, max);
    auto faces = IndexBufferFactory::generate();

    glGenVertexArrays(1, &_gpuMeshVAO);

    if (_gpuMeshVAO == GL_INVALID_VALUE)
    {
        throw std::runtime_error("CubeMesh: Error generating VAO");
    }

    bind();

    glGenBuffers(2, _buffers);

    if (_buffers[0] == GL_INVALID_VALUE || _buffers[1] == GL_INVALID_VALUE)
    {
        throw std::runtime_error("CubeMesh: Error generating VBOs");
    }

    glBindBuffer(GL_ARRAY_BUFFER, _buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), nullptr, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
    _enableAttributes();

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(uint32_t) * faces.size(), nullptr, GL_STATIC_DRAW);
    glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, sizeof(uint32_t) * faces.size(), faces.data());
}

void CubeMesh::bind() const
{
    glBindVertexArray(_gpuMeshVAO);
}

void CubeMesh::unbind() const
{
    glBindVertexArray(0);
}

void CubeMesh::destroy()
{
    glDeleteBuffers(2, _buffers);
    glDeleteVertexArrays(1, &_gpuMeshVAO);
}

void CubeMesh::_enableAttributes()
{
    static constexpr GLsizei vertexSize = 3 * sizeof(float);

    GLsizei stride = vertexSize;
    GLuint offset = 0;

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, static_cast<char *>(0) + (offset));

    glEnableVertexAttribArray(0);
    offset += vertexSize;
}
}