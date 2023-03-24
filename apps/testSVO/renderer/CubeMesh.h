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

#ifndef SIMCRUSHERRENDERER_CUBEMESH_H
#define SIMCRUSHERRENDERER_CUBEMESH_H

#include <glm/glm.hpp>

#include <glad/glad.h>

#include <vector>

namespace scr
{
  class CubeMesh
  {
    public:
      CubeMesh(const glm::vec3& min, const glm::vec3& max);

      void bind() const;
      void unbind() const;

      void destroy();

    private:
      void generate(std::vector<float>& vertexBuffer,
                    std::vector<uint32_t>& faceBuffer);
      void enableAttributes();
            
    private:
      const glm::vec3 _min;
      const glm::vec3 _max;

      GLuint _gpuMeshVAO;
      GLuint _buffers[2];
  };
}

#endif