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

# Locate header file
set(fmt_HEADER_SEARCH_DIRS
    "/usr/include"
    "/usr/local/include"
)

find_path(fmt_INCLUDE_DIR "fmt/format.h"
    PATHS ${fmt_HEADER_SEARCH_DIRS}
)

if(NOT ${fmt_INCLUDE_DIR})
    message(STATUS "Found fmt: ${fmt_INCLUDE_DIR}")
    add_definitions(-DFMT_HEADER_ONLY)
    set(fmt_INCLUDE_DIRS ${fmt_INCLUDE_DIR})
    include_directories(${fmt_INCLUDE_DIR})

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_FMT")
    set(ULTRALISER_USE_FMT TRUE)
else(NOT ${fmt_INCLUDE_DIR})
    message (STATUS "FMT NOT Found")
    set(ULTRALISER_USE_FMT FALSE)
endif(NOT ${fmt_INCLUDE_DIR})


