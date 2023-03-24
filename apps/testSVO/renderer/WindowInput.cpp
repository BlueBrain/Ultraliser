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

#include "WindowInput.h"

#include "Window.h"

#include <iostream>
#include <stdexcept>

namespace scr
{
  inline Window* getWindowHandle(GLFWwindow* window)
  {
    Window* win = reinterpret_cast<Window*>(glfwGetWindowUserPointer(window));
    if(!win)
      throw std::runtime_error("WindowInput: GLFWwindow User pointer is null");

    return win;
  }

  void keyboardCallbackFunc(SCR_UNUSED GLFWwindow * window, 
                            SCR_UNUSED int key, 
                            SCR_UNUSED int scanCode, 
                            SCR_UNUSED int action,
                            SCR_UNUSED int mods)
  {
     Window* win = getWindowHandle(window);
     win->onKeyboardPress(key, scanCode, action, mods);
  }

  void mouseInputCallbackFunc(SCR_UNUSED GLFWwindow * window, 
                              SCR_UNUSED int button, 
                              SCR_UNUSED int action,
                              SCR_UNUSED int mods)
  {
     Window* win = getWindowHandle(window);
     win->onMousePress(button, action, mods);
  } 

  void mouseMovementCallbackFunc(SCR_UNUSED GLFWwindow * window, 
                                 SCR_UNUSED double xpos, 
                                 SCR_UNUSED double ypos)
  {
     Window* win = getWindowHandle(window);
     win->onMouseMove(static_cast<int>(xpos), static_cast<int>(ypos));
  }

  void scrollCallbackFunc(SCR_UNUSED GLFWwindow * window,
                          SCR_UNUSED double xOffset,
                          SCR_UNUSED double yOffset)
  {
     Window* win = getWindowHandle(window);
     win->onMouseScroll(xOffset, yOffset);
  }

  void resizeCallbackFunc(SCR_UNUSED GLFWwindow * window, 
                          SCR_UNUSED int width, 
                          SCR_UNUSED int height)
  {
    Window* win = getWindowHandle(window);
    win->onWindowResize(width, height);
  }

  void errorCallback(SCR_UNUSED int errorCode, 
                     SCR_UNUSED const char * message)
  {
    std::cerr << "GLFW Error received!" << std::endl;
    std::cerr << "\tCode " << errorCode << std::endl;
    std::cerr << "\tMessage: " << message << std::endl;
  }
}