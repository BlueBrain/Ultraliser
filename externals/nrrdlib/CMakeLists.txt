####################################################################################################
# Copyright (c) 2015 - 2022
# Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
#
# Author(s)
#       Marwan Abdellah <marwan.abdellah@epfl.ch>
#       Nadir Roman Guerrero <nadir.romanguerrero@epfl.ch>
#
# For complete list of authors, please see AUTHORS.md
#
# This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
#
# This project is extended from Brayns with the same license of Ultraliser.
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

# nrrdloader
set(NRRD_LIB "nrrdlib")

# Collect the header and source files
file(GLOB_RECURSE NRRD_READER_HEADERS "*.h")
file(GLOB_RECURSE NRRD_READER_SOURCES "*.cpp")

# Add the current directory to include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Create the library
add_library(${NRRD_LIB} ${NRRD_READER_HEADERS} ${NRRD_READER_SOURCES})

# Copy the artifacts for the apps
set_property(TARGET ${NRRD_LIB} PROPERTY IMPORTED_LOCATION
    ${CMAKE_BINARY_DIRECTORY}/lib${NRRD_LIB}.so)
    
# Link the libraries
target_link_libraries(${NRRD_LIB} PRIVATE ${ZLIB_LIBRARIES} ${BZIP2_LIBRARIES})
target_include_directories(${NRRD_LIB} PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>"
    "$<INSTALL_INTERFACE:include>")
    
target_include_directories(${NRRD_LIB} SYSTEM PUBLIC ${GLM_INCLUDE_DIRS})

target_compile_options(${NRRD_LIB} PRIVATE -Wall -Wextra -pedantic -Werror)

# Install the files
install(TARGETS ${NRRD_LIB} DESTINATION lib)

install(DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/nrrdloader
    DESTINATION include
    FILES_MATCHING PATTERN "*.h")


