####################################################################################################
# Copyright (c) 2016 - 2023
# Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
#
# Author(s)
#       Marwan Abdellah <marwan.abdellah@epfl.ch>
#
# For complete list of authors, please see AUTHORS.md
#
# This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
#
# This library is free software; you can redistribute it and/or modify it under the terms of the
# GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
#
# This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# See the GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License along with this library;
# if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301 USA.
# You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
####################################################################################################

# Application name
set(MESH_APP_NAME ultraSkeletonizeVasculatureMesh)
set(VOLUME_APP_NAME ultraSkeletonizeVasculatureVolume)

# Add the current directory to include directories
include_directories(${CMAKE_BINARY_DIR}/include)

# Source files
file(GLOB ${MESH_APP_NAME}_SOURCE_FILES
    "SkeletonizeVasculatureMesh.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp")

file(GLOB ${VOLUME_APP_NAME}_SOURCE_FILES
    "SkeletonizeVasculatureVolume.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp")

# Header files
file(GLOB ${MESH_APP_NAME}_HEADER_FILES "*.h" "*.hh")
file(GLOB ${VOLUME_APP_NAME}_HEADER_FILES "*.h" "*.hh")

# Add the current directory to include directories
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../common/")

# Executable(s)
add_executable(${MESH_APP_NAME}
    ${${MESH_APP_NAME}_SOURCE_FILES}
    ${${MESH_APP_NAME}_HEADER_FILES}
)
add_executable(${VOLUME_APP_NAME}
    ${${VOLUME_APP_NAME}_SOURCE_FILES}
    ${${VOLUME_APP_NAME}_HEADER_FILES}
)

# Link against ultraliser
target_link_libraries(${MESH_APP_NAME} ultraliser)
target_link_libraries(${VOLUME_APP_NAME} ultraliser)

# Copy the executable after the build
add_custom_command(TARGET ${MESH_APP_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${MESH_APP_NAME}>
                   ${CMAKE_BINARY_DIR}/bin)
add_custom_command(TARGET ${VOLUME_APP_NAME} POST_BUILD
                  COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${VOLUME_APP_NAME}>
                  ${CMAKE_BINARY_DIR}/bin)

# Install the executable in the bin directory
install(TARGETS ${MESH_APP_NAME} DESTINATION bin)
install(TARGETS ${VOLUME_APP_NAME} DESTINATION bin)
