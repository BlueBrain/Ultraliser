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

# Ultraliser library
set(ULTRALISER_LIBRARY "ultraliser")

# Ultraliser files
set(ULTRALISER_SOURCE_FILES "")
set(ULTRALISER_HEADER_FILES "")

# A list of the directories that will be scanned
set(DIRECTORIES "algorithms"
                "algorithms/floodfill"
                "algorithms/geometry"
                "algorithms/mcs"
                "algorithms/simulation"
                "algorithms/skeletonization"
                "algorithms/skeletonization/thinning"
                "algorithms/utilities"
                "arguments"
                "common"
                "data"
                "data/common"
                "data/structures"
                "data/images"
                "data/meshes"
                "data/meshes/advanced"
                "data/meshes/advanced/math"
                "data/meshes/advanced/primitives"
                "data/meshes/advanced/io"
                "data/meshes/simple"
                "data/meshes/simple/primitives"
                "data/meshes/simple/io"
                "data/morphologies"
                "data/morphologies/h5"
                "data/morphologies/swc"
                "data/morphologies/vmv"
                "data/volumes"
                "data/volumes/grids"
                "data/volumes/utilities"
                "data/volumes/voxels"
                "geometry"
                "math"
                "microscopes"
                "utilities"
)

# Removing the include directory from the build one
file(REMOVE ${CMAKE_BINARY_DIR}/include)

# Add all the files in the project
foreach(DIRECTORY in ${DIRECTORIES})
    # Source files
    file(GLOB SOURCES
            "${DIRECTORY}/*.c" "${DIRECTORY}/*.cpp" "${DIRECTORY}/*.hpp")
    list(APPEND ULTRALISER_SOURCE_FILES ${SOURCES})

    # Header files
    file(GLOB HEADERS
            "${DIRECTORY}/*.h" "${DIRECTORY}/*.hh")
    list(APPEND ULTRALISER_HEADER_FILES ${HEADERS})

    # Copy the header files to the build directory for the apps
    file(COPY ${DIRECTORY}
            DESTINATION ${CMAKE_BINARY_DIR}/include
            FILES_MATCHING PATTERN "*.h" PATTERN "*.hh")
endforeach()

# Copy the interface
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Ultraliser.h
          ${CMAKE_BINARY_DIR}/include/Ultraliser.h COPYONLY)

# Add the current directory to include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Create a shared library
add_library(ultraliser ${ULTRALISER_SOURCE_FILES} ${ULTRALISER_HEADER_FILES})

# Link against nrrdlib, if compiled
if(${ULTRALISER_USE_NRRD})
    target_link_libraries(ultraliser nrrdlib)
endif(${ULTRALISER_USE_NRRD})

# Copy the artifacts for the apps
set_property(TARGET ultraliser PROPERTY IMPORTED_LOCATION
    ${CMAKE_BINARY_DIRECTORY}/lib${ULTRALISER_LIBRARY}.so)

# Install the files
install(TARGETS ${ULTRALISER_LIBRARY} DESTINATION lib)
install(FILES ${ULTRALISER_HEADER_FILES} DESTINATION include)
