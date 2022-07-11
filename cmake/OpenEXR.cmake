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

# Zlib
include(FindZLIB)

if(ZLIB_FOUND)
    # Library paths
    set(LIBRARY_PATHS
            /usr/lib
            /usr/local/lib
            /sw/lib
            /opt/local/lib
            $ENV{PROGRAM_FILES}/OpenEXR/lib/static)

    # Headers
    find_path(OPENEXR_INCLUDE_PATH ImfRgbaFile.h
        PATH_SUFFIXES OpenEXR
            /usr/include
            /usr/local/include
            /sw/include
            /opt/local/include)

    # Include directories
    include_directories(${OPENEXR_INCLUDE_PATH})
    include_directories(${OPENEXR_INCLUDE_PATH}/OpenEXR)
    include_directories(/usr/include)
    include_directories(/usr/local/include)
    include_directories(/opt/local/include)

    # Half
    find_library(OPENEXR_HALF_LIBRARY
        NAMES Half
        PATHS ${LIBRARY_PATHS})
    link_libraries(${OPENEXR_HALF_LIBRARY})

    # Iex
    find_library(OPENEXR_IEX_LIBRARY
        NAMES Iex
        PATHS ${LIBRARY_PATHS})
    link_libraries(${OPENEXR_IEX_LIBRARY})

    # Imath
    find_library(OPENEXR_IMATH_LIBRARY
        NAMES Imath
        PATHS ${LIBRARY_PATHS})
    link_libraries(${OPENEXR_IMATH_LIBRARY})

    # IlmImf
    find_library(OPENEXR_ILMIMF_LIBRARY
        NAMES IlmImf
        PATHS ${LIBRARY_PATHS})
    link_libraries(${OPENEXR_ILMIMF_LIBRARY})

    # IlmThread
    find_library(OPENEXR_ILMTHREAD_LIBRARY
        NAMES IlmThread
        PATHS ${LIBRARY_PATHS})
    link_libraries(${OPENEXR_ILMTHREAD_LIBRARY})
endif(ZLIB_FOUND)

if(OPENEXR_INCLUDE_PATH AND OPENEXR_IMATH_LIBRARY AND OPENEXR_ILMIMF_LIBRARY
        AND OPENEXR_IEX_LIBRARY AND OPENEXR_HALF_LIBRARY)

    set(OPENEXR_FOUND TRUE)
    set(OPENEXR_INCLUDE_PATHS ${OPENEXR_INCLUDE_PATH} CACHE STRING
        "The include paths needed to use OpenEXR")
    include_directories(${OPENEXR_INCLUDE_PATHS})

    # For BBP5 compilation
    include_directories(${ILMBASE_PACKAGE_PREFIX}/include/OpenEXR)
    link_directories(${ILMBASE_PACKAGE_PREFIX}/lib)

    # Libraries
    set(OPENEXR_LIBRARIES
            ${OPENEXR_IMATH_LIBRARY}
            ${OPENEXR_ILMIMF_LIBRARY}
            ${OPENEXR_IEX_LIBRARY}
            ${OPENEXR_HALF_LIBRARY}
            ${OPENEXR_ILMTHREAD_LIBRARY}
            ${ZLIB_LIBRARY} CACHE STRING "The libraries needed to use OpenEXR")
    link_libraries(${OPENEXR_LIBRARIES})

endif(OPENEXR_INCLUDE_PATH AND OPENEXR_IMATH_LIBRARY AND OPENEXR_ILMIMF_LIBRARY
    AND OPENEXR_IEX_LIBRARY AND OPENEXR_HALF_LIBRARY)

if(OPENEXR_FOUND)
    message(STATUS "Found OpenEXR: ${OPENEXR_ILMIMF_LIBRARY}")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DULTRALISER_USE_OPENEXR")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_OPENEXR")
else(OPENEXR_FOUND)
    message(STATUS "OpenEXR NOT Found")
endif(OPENEXR_FOUND)

# Advanced variables
mark_as_advanced(
    OPENEXR_INCLUDE_PATHS
    OPENEXR_LIBRARIES
    OPENEXR_ILMIMF_LIBRARY
    OPENEXR_IMATH_LIBRARY
    OPENEXR_IEX_LIBRARY
    OPENEXR_HALF_LIBRARY
)
