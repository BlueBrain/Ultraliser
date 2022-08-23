####################################################################################################
# Copyright (c) 2016 - 2022
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

# Applications
set(ANALYSIS_APPS "")
list(APPEND ANALYSIS_APPS AnalyseVessMorpho)

foreach(APP_NAME ${ANALYSIS_APPS})

    # Add the current directory to include directories
    include_directories(${CMAKE_BINARY_DIR}/include)

    # Add the current directory to include directories
    include_directories("${CMAKE_CURRENT_SOURCE_DIR}/../common/")

    # Source files
    file(GLOB COMMON_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/../common/*.cpp")

    # Executables
    add_executable("ultra${APP_NAME}" "${APP_NAME}.cpp" ${COMMON_SOURCE_FILES})

    # Link against ultraliser
    target_link_libraries("ultra${APP_NAME}" ultraliser)

    # Copy the executables after the build
    add_custom_command(TARGET "ultra${APP_NAME}" POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:ultra${APP_NAME}>
                       ${CMAKE_BINARY_DIR}/bin)

    # Install the executable in the bin directory
    install(TARGETS "ultra${APP_NAME}" DESTINATION bin)

endforeach()
