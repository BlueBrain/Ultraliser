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

#include "Shader.h"

#include <iostream>
#include <fstream>
#include <vector>

namespace scr
{
  Shader::Shader(const ShaderConfig& sc)
   : _oglShader(GL_INVALID_VALUE)
   , _oglvShader(GL_INVALID_VALUE)
   , _oglgShader(GL_INVALID_VALUE)
   , _oglfShader(GL_INVALID_VALUE)
  {
    compileAndLink(sc);
  }

  void Shader::bind() const
  {
    glUseProgram(_oglShader);
  }

  void Shader::unbind() const
  {
    glUseProgram(0);
  }

  void Shader::set(const std::string& name, const int value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniform1i(id, value);
  }

  void Shader::set(const std::string& name, const float value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniform1f(id, value);
  }

  void Shader::set(const std::string& name, const glm::vec3& value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniform3fv(id, 1, &value[0]);
  }

  void Shader::set(const std::string& name, const glm::vec4& value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniform4fv(id, 1, &value[0]);
  }

  void Shader::set(const std::string& name, const glm::mat3& value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniformMatrix3fv(id, 1, GL_FALSE, &(value[0][0]));
  }

  void Shader::set(const std::string& name, const glm::mat4& value) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniformMatrix4fv(id, 1, GL_FALSE, &(value[0][0]));
  }

  void Shader::set(const std::string& name, const GLuint64 handle) const
  {
    GLuint id;
    if(getUniformId(name, id))
      glUniformHandleui64ARB(id, handle);
  }

  void Shader::loadShaderCode(const std::string& filePath, std::string& buffer)
  {
    std::ifstream file;
    file.open(filePath, std::ios::in);
    if (!file)
    {
      std::string message = "Could not find shader file " + filePath;
      throw std::runtime_error(message.c_str());
    }

    char buf[0xfff];
    while (file.good())
    {
      file.getline(buf, 0xfff);
      buffer += std::string(buf) + "\n";  
    }

    file.close();
  }

  void Shader::destroy()
  {
    glDetachShader(_oglShader, _oglvShader);
    glDeleteShader(_oglvShader);
    glDetachShader(_oglShader, _oglfShader);
    glDeleteShader(_oglfShader);

    if(_oglgShader != GL_INVALID_VALUE)
    {
      glDetachShader(_oglShader, _oglgShader);
      glDeleteShader(_oglgShader);
    }

    glDeleteProgram(_oglShader);
  }

  GLuint Shader::compileShader(const std::string& filePath, GLenum shaderType)
  {
    std::string code;
    loadShaderCode(filePath, code);

    uint shaderId = glCreateShader(shaderType);

    const GLchar * codeCStr = static_cast<const GLchar*>(code.c_str());
    const GLint codeLen = static_cast<GLint>(code.length());
    glShaderSource(shaderId, 1, &codeCStr, &codeLen);

    glCompileShader(shaderId);

    GLint compiled;
    glGetShaderiv(shaderId, GL_COMPILE_STATUS, &compiled);
    if (!compiled)
    {
      GLint logLen;
      glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &logLen);
      std::vector<char> logString(logLen);
      glGetShaderInfoLog(shaderId, logLen, NULL, &logString[0]);
      glDeleteShader(shaderId);
      throw std::runtime_error(&logString[0]);
    }

    return shaderId;
  }

  void Shader::compileAndLink(const ShaderConfig& sc)
  {
    _oglShader = glCreateProgram();
    if(_oglShader == GL_INVALID_VALUE)
      throw std::runtime_error("Shader: Could not create program");
    
    _oglvShader = compileShader(sc.vertexShader, GL_VERTEX_SHADER);
    if(_oglvShader == GL_INVALID_VALUE)
      throw std::runtime_error("Shader: Could not create vertex shader");
    glAttachShader(_oglShader, _oglvShader);

    _oglfShader = compileShader(sc.fragmentShader, GL_FRAGMENT_SHADER);
    if(_oglfShader == GL_INVALID_VALUE)
      throw std::runtime_error("Shader: Could not create fragment shader");
    glAttachShader(_oglShader, _oglfShader);

    if(!sc.geometryShader.empty())
    {
      _oglgShader = compileShader(sc.geometryShader, GL_GEOMETRY_SHADER);
      if(_oglgShader == GL_INVALID_VALUE)
        throw std::runtime_error("Shader: Could not create geometry shader");

      glAttachShader(_oglShader, _oglgShader);
    }

    glLinkProgram(_oglShader);

    int linked;
    glGetProgramiv(_oglShader, GL_LINK_STATUS, &linked);
    if (!linked)
    {
      std::cerr << "Error linking shader program for:" << std::endl;
      std::cerr << "\tVertex shader: " << sc.vertexShader << std::endl;
      std::cerr << "\tGeometry shader: " << sc.geometryShader << std::endl;
      std::cerr << "\tFragment shader: " << sc.fragmentShader << std::endl;
      GLint logLen;
      glGetProgramiv(_oglShader, GL_INFO_LOG_LENGTH, &logLen);
      std::vector<char> logString(logLen);
      glGetProgramInfoLog(_oglShader, logLen, NULL, &logString[0]);
      throw std::runtime_error(&logString[0]);
    }

    parseInputs();
  }

  void Shader::parseInputs()
  {
    GLint i = 0;
    char nameBuffer[0xff];
    GLsizei nameLen = 0;
    GLint size = 0;
    GLenum typeInt;
    GLint active = 0;

    glGetProgramiv(_oglShader, GL_ACTIVE_UNIFORMS, &active);
    while(i < active)
    {
      glGetActiveUniform(_oglShader, i, 0xff, &nameLen, &size, &typeInt, nameBuffer);
      const std::string nameStr(nameBuffer);
      const GLuint id = glGetUniformLocation(_oglShader, nameBuffer);
      _uniforms[nameStr] = id;
      i++;
    }

    i = 0;
    glGetProgramiv(_oglShader, GL_ACTIVE_ATTRIBUTES, &active);
    while (i < active)
    {
      glGetActiveAttrib(_oglShader, i, 0xff, &nameLen, &size, &typeInt, nameBuffer);
      const std::string nameStr(nameBuffer);
      const GLuint id = glGetAttribLocation(_oglShader, nameBuffer);
      _attribs[nameStr] = id;
      i++;
    }
  }

  inline bool Shader::getUniformId(const std::string& name, GLuint& buf) const
  {
    auto it = _uniforms.find(name);
    if(it != _uniforms.end())
    {
      buf = it->second;
      return true;
    }
    return false;
  }
}