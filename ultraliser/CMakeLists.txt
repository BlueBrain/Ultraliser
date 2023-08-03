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

# Removing the include directory from the build one
file(REMOVE_RECURSE ${CMAKE_BINARY_DIR}/include)

# Get all the cource files
FILE(GLOB_RECURSE ULTRALISER_SOURCE_FILES *.cpp *.c)
FILE(GLOB_RECURSE ULTRALISER_HEADER_FILES_FP *.h *.hh *.hpp)

#file(GLOB V_GLOB LIST_DIRECTORIES true "*")
#foreach(item ${V_GLOB})
#    message(${item})
#endforeach()

# Get the correct paths of the header files
foreach(HEADER_FILE ${ULTRALISER_HEADER_FILES_FP})
    string(REPLACE ${CMAKE_CURRENT_SOURCE_DIR}/ "" RELATIVE_HEADER ${HEADER_FILE})
    list(APPEND ULTRALISER_HEADER_FILES ${RELATIVE_HEADER})
endforeach()

# Copy the header files to the actual destinations in the build directory
foreach(item ${ULTRALISER_HEADER_FILES})
    get_filename_component(FILE_NAME ${item} NAME)
    get_filename_component(DIRECTORY_NAME ${item} DIRECTORY)
    file(COPY ${item} DESTINATION ${CMAKE_BINARY_DIR}/include/${DIRECTORY_NAME})
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
