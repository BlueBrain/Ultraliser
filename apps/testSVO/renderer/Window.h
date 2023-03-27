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

#include <GLFW/glfw3.h>

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "Camera.h"
#include "Model.h"

namespace svorender
{
class Window
{
public:
    Window(const uint32_t width, const uint32_t height, const std::string &title = "SVO Visualizer");

    void onWindowResize(uint32_t width, uint32_t height);
    void onKeyboardPress(char key, int scanMode, int action, int mods);
    void onMousePress(int button, int action, int mods);
    void onMouseMove(int xpos, int ypos);
    void onMouseScroll(float xoffset, float yoffset);

    void renderLoop();

    uint32_t width() const;
    uint32_t height() const;

    void setCamera(const float nearPlane, const float farPlane, const float fov);
    void addModel(std::unique_ptr<Model> model);

    Camera &getCamera();
    const Camera &getCamera() const;

private:
    void _initContext(const std::string &windowTitle);
    void _initOpenGL();

private:
    uint32_t _winWidth;
    uint32_t _winHeight;

    GLFWwindow *_window;

    Camera _camera;
    std::vector<std::unique_ptr<Model>> _models;

    int _mouseButtonPressed;
    int _lastMouseX;
    int _lastMouseY;
};
}