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

#include <iostream>

namespace
{
struct CameraConstants
{
    static constexpr float movementSpeed = 0.5f;
};

svorender::Window *getWindowHandle(GLFWwindow *window)
{
    auto win = reinterpret_cast<svorender::Window *>(glfwGetWindowUserPointer(window));
    if (!win)
    {
        throw std::runtime_error("WindowInput: GLFWwindow User pointer is null");
    }

    return win;
}

void keyboardCallbackFunc(GLFWwindow *window, int key, int scanCode, int action, int mods)
{
    auto win = getWindowHandle(window);
    win->onKeyboardPress(key, scanCode, action, mods);
}

void mouseInputCallbackFunc(GLFWwindow *window, int button, int action, int mods)
{
    auto win = getWindowHandle(window);
    win->onMousePress(button, action, mods);
}

void mouseMovementCallbackFunc(GLFWwindow *window, double xpos, double ypos)
{
    auto win = getWindowHandle(window);
    win->onMouseMove(static_cast<int>(xpos), static_cast<int>(ypos));
}

void scrollCallbackFunc(GLFWwindow *window, double xOffset, double yOffset)
{
    auto win = getWindowHandle(window);
    win->onMouseScroll(xOffset, yOffset);
}

void resizeCallbackFunc(GLFWwindow *window, int width, int height)
{
    auto win = getWindowHandle(window);
    win->onWindowResize(width, height);
}

void errorCallback(int errorCode, const char *message)
{
    std::cerr << "GLFW Error received!" << std::endl;
    std::cerr << "\tCode " << errorCode << std::endl;
    std::cerr << "\tMessage: " << message << std::endl;
}
}

namespace svorender
{

Window::Window(const uint32_t width, const uint32_t height, const std::string &title)
    : _winWidth(width)
    , _winHeight(height)
    , _window(nullptr)
    , _camera(1.f, 10000.f, 30.f)
    , _mouseButtonPressed(-1)
    , _lastMouseX(-1)
    , _lastMouseY(-1)
{
    _initContext(title);
    _initOpenGL();
    _camera.onScreenResize(width, height);
}

void Window::onWindowResize(uint32_t width, uint32_t height)
{
    _winWidth = width;
    _winHeight = height;
    _camera.onScreenResize(width, height);
    glViewport(0, 0, width, height);
}

void Window::onKeyboardPress(char key, int scanMode, int action, int mods)
{
    (void)key;
    (void)scanMode;
    (void)action;
    (void)mods;
}

void Window::onMousePress(int button, int action, int mods)
{
    (void)mods;
    switch (action)
    {
    case GLFW_PRESS:
        _mouseButtonPressed = button;
        break;
    case GLFW_RELEASE:
        _mouseButtonPressed = (button == _mouseButtonPressed) ? -1 : _mouseButtonPressed;
        _lastMouseX = -1;
        _lastMouseY = -1;
        break;
    }
}

void Window::onMouseMove(int xpos, int ypos)
{
    if (_mouseButtonPressed != 0)
        return;

    if (_lastMouseX == -1 || _lastMouseY == -1)
    {
        _lastMouseX = xpos;
        _lastMouseY = ypos;
        return;
    }

    auto deltaX = xpos - _lastMouseX;
    auto deltaY = ypos - _lastMouseY;
    auto deltaRotX = -static_cast<float>(deltaY) * CameraConstants::movementSpeed * 0.01;
    auto deltaRotY = -static_cast<float>(deltaX) * CameraConstants::movementSpeed * 0.01;
    auto deltaRot = glm::vec3(deltaRotX, deltaRotY, 0.f);
    _camera.updateRotation(deltaRot);

    _lastMouseX = xpos;
    _lastMouseY = ypos;
}

void Window::onMouseScroll(float xoffset, float yoffset)
{
    (void)xoffset;
    _camera.onZoom(static_cast<float>(yoffset));
}

void Window::renderLoop()
{
    while (!glfwWindowShouldClose(_window))
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        _camera.update();

        for (auto &model : _models)
        {
            model->render(_camera);
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

void Window::setCamera(const float nearPlane, const float farPlane, const float fov)
{
    _camera = Camera(nearPlane, farPlane, fov);
}

void Window::addModel(std::unique_ptr<Model> model)
{
    _models.push_back(std::move(model));
}

Camera &Window::getCamera()
{
    return _camera;
}

const Camera &Window::getCamera() const
{
    return _camera;
}

void Window::_initContext(const std::string &windowTitle)
{
    if (!glfwInit())
    {
        throw std::runtime_error("Window: Could not initialize GLFW");
    }

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

    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        glfwTerminate();
        throw std::runtime_error("Window: Could not load OpenGL");
    }
}

void Window::_initOpenGL()
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