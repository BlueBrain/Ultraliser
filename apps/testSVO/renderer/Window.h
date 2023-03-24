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

#ifndef SIMCRUSHERRENDERER_WINDOW_H
#define SIMCRUSHERRENDERER_WINDOW_H

#include <glad/glad.h>

#include <GLFW/glfw3.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "Common.h"
#include "Camera.h"
#include "Model.h"

namespace scr
{
  class Window
  {
    public:
      Window(const uint32_t width, 
              const uint32_t height, 
              const std::string& title = "SimCrusher Volume Visual Debugger");

      void onWindowResize(const uint32_t width, const uint32_t height);
      void onKeyboardPress(SCR_UNUSED const char key, 
                           SCR_UNUSED const int scanMode, 
                           SCR_UNUSED const int action, 
                           SCR_UNUSED const int mods);
      void onMousePress(const int button, const int action, const int mods);
      void onMouseMove(const int xpos, const int ypos);
      void onMouseScroll(SCR_UNUSED const float xoffset, 
                         SCR_UNUSED const float yoffset);

      void renderLoop();

      uint32_t width() const;
      uint32_t height() const;

      void setCamera(const float nearPlane, 
                     const float farPlane, 
                     const float fov);
      void addModel(CubeMesh* mesh, Material* material);

      Camera& getCamera();
      const Camera& getCamera() const;
    
    private:
      void initContext(const std::string& windowTitle);
      void initOpenGL();

    private:
      uint32_t _winWidth;
      uint32_t _winHeight;

      GLFWwindow * _window;

      std::unique_ptr<Camera> _camera;
      std::vector<std::unique_ptr<Model>> _models;

      int _mouseButtonPressed;
      int _lastMouseX;
      int _lastMouseY;
      
  };
}

#endif