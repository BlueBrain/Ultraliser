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

# HDF5
find_package(HDF5 COMPONENTS CXX HL)

if(HDF5_FOUND)
    message(STATUS "Found HDF5: ${HDF5_INCLUDE_DIRS}, ${HDF5_LIBRARIES}")

    # Library paths
    set(LIBRARY_PATHS
            /usr/lib
            /usr/local/lib
            /sw/lib
            /opt/local/lib)

    # Include directories
    include_directories(${HDF5_INCLUDE_DIRS})
    include_directories(/usr/include)
    include_directories(/usr/local/include)
    include_directories(/opt/local/include)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_H5")
    set(ULTRALISER_USE_H5 TRUE)

    link_libraries(${HDF5_LIBRARIES}
        ${HDF5_HL_LIBRARIES} ${HDF5_C_LIBRARIES} ${HDF5_CXX_LIBRARIES}
        ${HDF5_C_HL_LIBRARIES} ${HDF5_CXX_HL_LIBRARIES})
else(HDF5_FOUND)
     message(STATUS "HDF5 NOT Found")
     set(ULTRALISER_USE_H5 FALSE)
endif(HDF5_FOUND)

