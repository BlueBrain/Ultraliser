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

# ultraMask2Mesh
if(${ULTRALISER_USE_TIFF})
    add_subdirectory(ultraMask2Mesh)
else(${ULTRALISER_USE_TIFF})
    message(STATUS "ultraMask2Mesh will not be compiled. Missing the TIFF library!")
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

# ultraScaleMesh
add_subdirectory(ultraScaleMesh)

# ultraSplitPartitions
add_subdirectory(ultraSplitPartitions)

# ultraMeshQualityChecker
add_subdirectory(ultraMeshQualityChecker)

# ultraH5MeshConverter
if(${ULTRALISER_USE_H5})
add_subdirectory(ultraH5MeshConverter)
else(${ULTRALISER_USE_H5})
    message(STATUS "ultraH5MeshConverter will not be compiled. Missing H5!")
endif(${ULTRALISER_USE_H5})

# ultraVessMorpho2Mesh
add_subdirectory(ultraVessMorpho2Mesh)

# ultraAstroMorpho2Mesh
add_subdirectory(ultraAstroMorpho2Mesh)

# ultraNeuroMorpho2Mesh
add_subdirectory(ultraNeuroMorpho2Mesh)

# ultraSampleVascularMorphology
add_subdirectory(ultraSampleVascularMorphology)

# ultraAnalysis
add_subdirectory(ultraAnalysis)

# ultraResampling
add_subdirectory(ultraResampling)

include_directories(${CMAKE_BINARY_DIR}/include)
file(GLOB testEdgeCollapse_HEADER_FILES "*.h" "*.hh")
add_executable(testEdgeCollapse testEdgeCollapse.cpp ${testEdgeCollapse_HEADER_FILES})
target_link_libraries(testEdgeCollapse ultraliser)
add_custom_command(TARGET testEdgeCollapse POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:testEdgeCollapse>
                   ${CMAKE_BINARY_DIR}/bin)
