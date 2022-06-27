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

# ZLIB
find_package( BZip2 REQUIRED )

if(BZIP2_FOUND)
    message(STATUS "Found ZLIB: ${BZIP2_INCLUDE_DIR}, ${BZIP2_LIBRARIES}")
    include_directories(${BZIP2_INCLUDE_DIR})
    link_libraries(${BZIP2_LIBRARIES})
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_ZLIB")
    #set(ULTRALISER_USE_ZLIB TRUE)
else(BZIP2_FOUND)
    message(STATUS "BZip2 NOT Found")
endif(BZIP2_FOUND)