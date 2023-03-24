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

#ifndef SIMCRUSHERRENDERER_SHADER_H
#define SIMCRUSHERRENDERER_SHADER_H

#include <glad/glad.h>

#include <glm/glm.hpp>

#include <string>
#include <cstdint>
#include <functional>
#include <unordered_map>

namespace scr
{
  struct ShaderConfig
  {
    std::string vertexShader;
    std::string geometryShader;
    std::string fragmentShader;
  };

  class Shader
  {
    public:
      Shader(const ShaderConfig& sc);

      void bind() const;
      void unbind() const;

      void set(const std::string& name, const int value) const;
      void set(const std::string& name, const float value) const;
      void set(const std::string& name, const glm::vec3& value) const;
      void set(const std::string& name, const glm::vec4& value) const;
      void set(const std::string& name, const glm::mat3& value) const;
      void set(const std::string& name, const glm::mat4& value) const;
      void set(const std::string& name, const GLuint64 handle) const;

      void destroy();

    private:
      void loadShaderCode(const std::string& filePath, std::string& buffer);
      GLuint compileShader(const std::string& filePath, GLenum shaderType);
      void compileAndLink(const ShaderConfig& sc);
      void parseInputs();
      bool getUniformId(const std::string& name, GLuint& buf) const;

    private:
      GLuint _oglShader;
      GLuint _oglvShader;
      GLuint _oglgShader;
      GLuint _oglfShader;

      std::unordered_map<std::string, GLuint> _attribs;
      std::unordered_map<std::string, GLuint> _uniforms;
      
  };
}

#endif