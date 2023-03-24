# Copyright (c) 2020, EPFL/Blue Brain Project
# All rights reserved. Do not distribute without permission.
# Responsible Author: Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
#
# This file is part of SimCrusher
# <LINK>
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License version 3.0 as published
# by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this library; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
find_package(glfw3 REQUIRED)
find_package(glm REQUIRED)
add_subdirectory(glad)

set(VOXELVIEWER_SOURCES
  testSVO.cpp
  renderer/Camera.cpp
  renderer/CubeMesh.cpp
  renderer/Material.cpp
  renderer/Model.cpp
  renderer/Shader.cpp
  renderer/Texture3D.cpp
  renderer/Window.cpp
  renderer/WindowInput.cpp
  renderer/materials/BoundRenderMaterial.cpp
  renderer/materials/VolumeMaterial.cpp
)

set(VOXELVIEWER_HEADERS
  renderer/Camera.h
  renderer/Common.h
  renderer/CubeMesh.h
  renderer/Material.h
  renderer/Model.h
  renderer/Shader.h
  renderer/Texture3D.h
  renderer/Window.h
  renderer/WindowInput.h
  renderer/materials/BoundRenderMaterial.h
  renderer/materials/VolumeMaterial.h
)

include_directories(${CMAKE_BINARY_DIR}/include)

add_executable(voxelviewer ${VOXELVIEWER_SOURCES} ${VOXELVIEWER_HEADERS})

target_link_libraries(voxelviewer ultraliser)
target_link_libraries(voxelviewer ${CMAKE_DL_LIBS} glad glfw)

target_include_directories(voxelviewer PRIVATE ${GLM_INCLUDE_DIRS})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/shaders DESTINATION ${CMAKE_BINARY_DIR})

add_custom_command(TARGET voxelviewer POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:voxelviewer>
  ${CMAKE_BINARY_DIR}/bin)