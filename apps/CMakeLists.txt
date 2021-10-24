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

# ultraMask2Mesh
if(${ULTRALISER_USE_TIFF})
    add_subdirectory(ultraMask2Mesh)
else(${ULTRALISER_USE_TIFF})
    message(STATUS "ultraMask2Mesh will not be compiled. Missing TIFF!")
endif(${ULTRALISER_USE_TIFF})

# ultraMesh2Mesh
add_subdirectory(ultraMesh2Mesh)

# ultraMeshes2Mesh
add_subdirectory(ultraMeshes2Mesh)

# ultraVolume2Mesh
add_subdirectory(ultraVolume2Mesh)

# ultraMesh2Volume
add_subdirectory(ultraMesh2Volume)

# ultraCleanMesh
add_subdirectory(ultraCleanMesh)

# ultraSplitPartitions
add_subdirectory(ultraSplitPartitions)

# ultraMeshQualityChecker
add_subdirectory(ultraMeshQualityChecker)

# ultraVasculatureMorphology2Mesh
if(${ULTRALISER_USE_H5})
    add_subdirectory(ultraVasculatureMorphology2Mesh)
else(${ULTRALISER_USE_H5})
    message(STATUS "ultraVasculatureMorphology2Mesh will not be compiled. Missing H5!")
endif(${ULTRALISER_USE_H5})



