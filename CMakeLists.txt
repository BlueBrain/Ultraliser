####################################################################################################
# Copyright (c) 2016 - 2021
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

# CMake 3.5 is required
cmake_minimum_required(VERSION 3.5)

# Ultraliser
project (Ultraliser)

#Add CMake customized modules
set(CMAKE_MODULE_PATH  ${CMAKE_SOURCE_DIR}/cmake)
include(OpenMP)
include(TIFF)
include(HDF5)
include(Eigen3)
include(GLM)
include(ZLIB)
include(BZip2)

#include(Doxygen)

include_directories("/usr/include")

# C++ 14
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Installation directories
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# C++ Flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 \
    -Wno-variadic-macros \
    -Wno-int-to-pointer-cast \
    -Wno-deprecated \
    -Wno-unreachable-code \
    -Wno-old-style-cast \
    -Wno-error=format-security \
    -Wno-float-equal")

set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_RELEASE")

OPTION(ENABLE_PROGRESS_BAR "Disable the progress bar to gain performance" ON)
if(ENABLE_PROGRESS_BAR)
  add_definitions(-DENABLE_PROGRESS_BAR)
endif()

# NRRD library, from Brayns
add_subdirectory(externals)

# Ultraliser library
add_subdirectory(ultraliser)

# Ultraliser applications
add_subdirectory(apps)

# Adding documentation
# add_subdirectory(doc)
