####################################################################################################
# Copyright (c) 2016 - 2022
# Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
#
# Author(s)
#       Juan Jose Garcia Cantero <juanjose.garcia@epfl.ch>
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

# Eigen3
find_package(Eigen3 REQUIRED)

if(EIGEN3_FOUND)
    message(STATUS "Found EIGEN3: ${EIGEN3_INCLUDE_DIRS}")

    # Include directories
    include_directories(${EIGEN3_INCLUDE_DIRS})
    include_directories(/usr/include)
    include_directories(/usr/local/include)
    include_directories(/opt/local/include)

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DULTRALISER_USE_EIGEN3")
    set(ULTRALISER_USE_EIGEN3 TRUE)
else(EIGEN3_FOUND)
     message(STATUS "EIGEN3 NOT Found")
endif(EIGEN3_FOUND)