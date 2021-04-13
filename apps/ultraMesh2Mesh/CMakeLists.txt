####################################################################################################
# Copyright (c) 2016 - 2021
# Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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
set(APP_NAME ultraMesh2Mesh)

# Add the current directory to include directories
include_directories(${CMAKE_BINARY_DIR}/include)

# Source files
file(GLOB ${APP_NAME}_SOURCE_FILES "*.cpp")

# Header files
file(GLOB ${APP_NAME}_HEADER_FILES "*.h" "*.hh")

# Executable
add_executable(${APP_NAME}
    ${${APP_NAME}_SOURCE_FILES}
    ${${APP_NAME}_SOURCE_FILES}
)

# Link against ultraliser
target_link_libraries(${APP_NAME} ultraliser)

# Copy the executable after the build
add_custom_command(TARGET ${APP_NAME} POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${APP_NAME}>
                   ${CMAKE_BINARY_DIR}/bin)

# Install the executable in the bin directory
install(TARGETS ${APP_NAME} DESTINATION bin)
