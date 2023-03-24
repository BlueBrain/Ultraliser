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

#include "Window.h"

#include "WindowInput.h"

namespace scr
{

  Window::Window(const uint32_t width, 
                 const uint32_t height, 
                 const std::string& title)
   : _winWidth(width)
   , _winHeight(height)
   , _window(nullptr)
   , _camera(std::unique_ptr<Camera>(new Camera(1.f, 500.f, 30.f)))
   , _mouseButtonPressed(-1)
   , _lastMouseX(-1)
   , _lastMouseY(-1)
  {
    initContext(title);
    initOpenGL();
    _camera->onScreenResize(width, height);
  }

  void Window::onWindowResize(const uint32_t width, const uint32_t height)
  {
    _winWidth = width;
    _winHeight = height;

    _camera->onScreenResize(width, height);

    glViewport(0, 0, width, height);
  }

  void Window::onKeyboardPress(SCR_UNUSED const char key, 
                               SCR_UNUSED const int scanMode, 
                               SCR_UNUSED const int action, 
                               SCR_UNUSED const int mods)
  {
 
  }

  void Window::onMousePress(const int button, const int action, const int mods)
  {
    (void)mods;
    switch (action)
    {
      case GLFW_PRESS:
        _mouseButtonPressed = button;
        break;
      case GLFW_RELEASE:
        _mouseButtonPressed = (button == _mouseButtonPressed)? 
                              -1 : _mouseButtonPressed;
        _lastMouseX = -1;
        _lastMouseY = -1;
        break;
    }
  }

  void Window::onMouseMove(const int xpos, const int ypos)
  {
    if(_mouseButtonPressed != 0)
      return;

    if(_lastMouseX == -1 || _lastMouseY == -1)
    {
      _lastMouseX = xpos;
      _lastMouseY = ypos;
      return;
    }

    const int deltaX = xpos - _lastMouseX;
    const int deltaY = ypos - _lastMouseY;

    const float deltaRotX= -static_cast<float>(deltaY) * CAMERA_MOV_SPEED * 0.01  ;
    const float deltaRotY= -static_cast<float>(deltaX) * CAMERA_MOV_SPEED * 0.01;

    const glm::vec3 deltaRot (deltaRotX, deltaRotY, 0.f);
    _camera->updateRotation(deltaRot);

    _lastMouseX = xpos;
    _lastMouseY = ypos;
  }

  void Window::onMouseScroll(SCR_UNUSED const float xoffset, 
                             SCR_UNUSED const float yoffset)
  {

      (void)xoffset;
      _camera->onZoom(static_cast<float>(yoffset));
  }

  void Window::renderLoop()
  {
    while(!glfwWindowShouldClose(_window))
    {
      glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
      
      _camera->update();

      for(const auto& model : _models)
      {
        model->bind();
        model->renderModel(*_camera.get());
      }

      glfwPollEvents();
      glfwSwapBuffers(_window);
    }
  }

  uint32_t Window::width() const
  {
    return _winWidth;
  }

  uint32_t Window::height() const
  {
    return _winHeight;
  }

  void Window::setCamera(const float nearPlane, 
                         const float farPlane, 
                         const float fov)
  {
    _camera.reset(new Camera(nearPlane, farPlane, fov));
  }

  void Window::addModel(CubeMesh* mesh, Material* material)
  {
    _models.push_back(std::unique_ptr<Model>(new Model(mesh, material)));
  }

  Camera& Window::getCamera()
  {
    return *_camera.get();
  }

  const Camera& Window::getCamera() const
  {
    return *_camera.get();
  }

  void Window::initContext(const std::string& windowTitle)
  {
    if(!glfwInit())
      throw std::runtime_error("Window: Could not initialize GLFW");
    
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);

    _window = glfwCreateWindow(_winWidth, _winHeight, windowTitle.c_str(), NULL, NULL);
    if (_window == nullptr)
    {
        glfwTerminate();
        throw std::runtime_error("Window: Could not create GLFW window");
    }

    glfwMakeContextCurrent(_window);

    glfwSetKeyCallback(_window, keyboardCallbackFunc);
    glfwSetCursorPosCallback(_window, mouseMovementCallbackFunc);
    glfwSetMouseButtonCallback(_window, mouseInputCallbackFunc);
    glfwSetScrollCallback(_window, scrollCallbackFunc);
    glfwSetFramebufferSizeCallback(_window, resizeCallbackFunc);
    glfwSetErrorCallback(errorCallback);

    glfwSetWindowUserPointer(_window, this);
    glfwSetWindowPos(_window, 0, 0);

    glfwSwapInterval(1);

    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
    {
        glfwTerminate();
        throw std::runtime_error("Window: Could not load OpenGL");
    }
  }

  void Window::initOpenGL()
  {
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glFrontFace(GL_CCW);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glEnable(GL_CULL_FACE);
    glViewport(0, 0, _winWidth, _winHeight);
  }
}