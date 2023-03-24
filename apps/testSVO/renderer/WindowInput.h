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

#ifndef SIMCRUSHERRENDERER_WINDOWINPUT_H
#define SIMCRUSHERRENDERER_WINDOWINPUT_H

#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "Common.h"

namespace scr
{
  void keyboardCallbackFunc(SCR_UNUSED GLFWwindow * window_, 
                            SCR_UNUSED int key_, 
                            SCR_UNUSED int scanCode_, 
                            SCR_UNUSED int action_,
                            SCR_UNUSED int mods_);

  void mouseInputCallbackFunc(SCR_UNUSED GLFWwindow * window_, 
                              SCR_UNUSED int button_, 
                              SCR_UNUSED int action_,
                              SCR_UNUSED int mods_);

  void mouseMovementCallbackFunc(SCR_UNUSED GLFWwindow * window_, 
                                 SCR_UNUSED double xpos_, 
                                 SCR_UNUSED double ypos_);

  void scrollCallbackFunc(SCR_UNUSED GLFWwindow * window_,
                          SCR_UNUSED double xOffset_,
                          SCR_UNUSED double yOffset_);

  void resizeCallbackFunc(SCR_UNUSED GLFWwindow * window_, 
                          SCR_UNUSED int width_, 
                          SCR_UNUSED int height_);

  void errorCallback(SCR_UNUSED int errorCode_, 
                     SCR_UNUSED const char * message_);
}

#endif