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

#include "Common.h"

#include <cstring>
#include <stdexcept>

namespace scr
{
  CubeMesh::CubeMesh(const glm::vec3& min, const glm::vec3& max)
   : _min(min)
   , _max(max)
   , _gpuMeshVAO(GL_INVALID_VALUE)
   , _buffers{GL_INVALID_VALUE, GL_INVALID_VALUE}
  {
    std::vector<float> vertices;
    std::vector<uint32_t> faces;
    generate(vertices, faces);

    glGenVertexArrays(1, &_gpuMeshVAO);

    if(_gpuMeshVAO == GL_INVALID_VALUE)
      throw std::runtime_error("CubeMesh: Error generating VAO");

    bind();

    glGenBuffers(2, _buffers);

    if(_buffers[0] == GL_INVALID_VALUE || _buffers[1] == GL_INVALID_VALUE)
      throw std::runtime_error("CubeMesh: Error generating VBOs");
    
    glBindBuffer(GL_ARRAY_BUFFER, _buffers[0]);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * vertices.size(), nullptr, GL_STATIC_DRAW);
    glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * vertices.size(), vertices.data());
    enableAttributes();

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

  void CubeMesh::generate(std::vector<float>& vertexBuffer,
                          std::vector<uint32_t>& faceBuffer)
  {
    vertexBuffer.resize(24);

    // Vertices
    const glm::vec3 v[8] = {
     {_min.x, _min.y, _max.z},
     {_max.x, _min.y, _max.z},
     _min,
     {_max.x, _min.y, _min.z},
     {_min.x, _max.y, _max.z},
     _max,
     {_min.x, _max.y, _min.z},
     {_max.x, _max.y, _min.z}
    };

    for(int i = 0; i < 8; i++)
      memcpy(&vertexBuffer[0] + 3 * i, &(v[i][0]), sizeof(float) * 3);


    // Faces
    faceBuffer.reserve(36);
    
    faceBuffer.push_back(0); faceBuffer.push_back(2); faceBuffer.push_back(1);
    faceBuffer.push_back(1); faceBuffer.push_back(2); faceBuffer.push_back(3);
    faceBuffer.push_back(5); faceBuffer.push_back(7); faceBuffer.push_back(4);
    faceBuffer.push_back(4); faceBuffer.push_back(7); faceBuffer.push_back(6);
    faceBuffer.push_back(3); faceBuffer.push_back(5); faceBuffer.push_back(1);
    faceBuffer.push_back(5); faceBuffer.push_back(3); faceBuffer.push_back(7);
    faceBuffer.push_back(4); faceBuffer.push_back(6); faceBuffer.push_back(0);
    faceBuffer.push_back(0); faceBuffer.push_back(6); faceBuffer.push_back(2);
    faceBuffer.push_back(4); faceBuffer.push_back(0); faceBuffer.push_back(5);
    faceBuffer.push_back(5); faceBuffer.push_back(0); faceBuffer.push_back(1);
    faceBuffer.push_back(2); faceBuffer.push_back(6); faceBuffer.push_back(3);
    faceBuffer.push_back(3); faceBuffer.push_back(6); faceBuffer.push_back(7);
  }

  void CubeMesh::enableAttributes()
  {
    static constexpr GLsizei vertexSize = 3 * sizeof(float);

    GLsizei stride = vertexSize;
    GLuint offset = 0;

    glVertexAttribPointer(0, 
                          3,
                          GL_FLOAT,
                          GL_FALSE,
                          stride,
                          BUFFER_OFFSET(offset));

    glEnableVertexAttribArray(0);
    offset += vertexSize;
  }
}