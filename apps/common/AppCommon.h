/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
 * MA 02111-1307, USA.
 * You can also find it on the GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <Ultraliser.h>
#include <AppOptions.h>

namespace Ultraliser
{

/**
 * @brief createVolumeGrid
 * Creates a volume grid compatiable with the input mesh.
 *
 * @param mesh
 * Input mesh.
 * @param options
 * System options.
 * @return
 * Volume base resolution (along the largest dimension).
 */
Volume* createVolumeGrid(Mesh *mesh, const AppOptions* options);

/**
 * @brief computeBoundingBoxForMeshes
 * Computes the bounding box for a group of meshes.
 *
 * @param boundsFile
 * Input bounding file containg a region of interest.
 * @param inputMeshesDirectory
 * A directory that contains all list of meshes.
 * @param meshFiles
 * A list of all the meshes.
 * @param pMax
 * Return pMax.
 * @param pMin
 * Return pMin.
 */
void computeBoundingBoxForMeshes(const std::string& boundsFile,
                                 const std::string& inputMeshesDirectory,
                                 std::vector< std::string > meshFiles,
                                 Vector3f& pMax, Vector3f& pMin,
                                 const float& xScale,const float& yScale,const float& zScale);

/**
 * @brief removeUnwantedPartitions
 * @param mesh
 * @param options
 * @return
 */
Mesh* removeUnwantedPartitions(Mesh* mesh, const AppOptions* options);

/**
 * @brief loadInputMesh
 * Loads an input mesh.
 * @param options
 * System options.
 * @return
 * Pointer to the loaded input mesh.
 */
Mesh* loadInputMesh(const AppOptions* options);

/**
 * @brief optimizeMesh
 * @param mesh
 * @param options
 */
void optimizeMesh(Mesh *mesh, const AppOptions* options);

/**
 * @brief ensureWatertightness
 * @param mesh
 * @param options
 */
void ensureWatertightness(Mesh *mesh, const AppOptions* options);

/**
 * @brief optimizeMeshWithPartitions
 * @param mesh
 * @param options
 */
void optimizeMeshWithPartitions(AdvancedMesh* mesh, const AppOptions* options);

/**
 * @brief applyLaplacianOperator
 * Apply the laplacian operator on the mesh.
 *
 * @param mesh
 * An input mesh that is guaranteed to be two-manifold.
 * @param options
 * User-defined options given to the executable.
 */
void applySmoothingOperator(Mesh *mesh, const AppOptions* options);

/**
 * @brief createMeshWithNoSelfIntersections
 * Process the two manifold mesh and create a watertight mesh using the AdvancedMesh processing.
 *
 * @param manifoldMesh
 * An input mesh that is guaranteed to be two-manifold but have some self intersections.
 * @param options
 * User-defined options given to the executable.
 */
void createWatertightMesh(const Mesh* manifoldMesh, const AppOptions* options);

/**
 * @brief optimizeMesh
 * Optimize the resulting mesh from the DMC algorithm.
 *
 * @param dmcMesh
 * The resulting mesh from the DMC algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateOptimizedMesh(Mesh *dmcMesh, const AppOptions* options);

/**
 * @brief generateOptimizedMeshWithROI
 * @param mesh
 * @param options
 * @param regions
 */
void generateOptimizedMeshWithROI(Mesh *mesh, const AppOptions* options, const ROIs& regions);

/**
 * @brief generateMarchingCubesMeshArtifacts
 * Writes the resulting mesh from the marching cubes algorithm.
 *
 * @param mesh
 * THe generetd mesh from the marching cubes algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateMarchingCubesMeshArtifacts(const Mesh *mesh, const AppOptions* options);

/**
 * @brief generateMeshArtifacts
 * Writes the resulting mesh artifacts.
 *
 * @param mesh
 * Reconstructed mesh from the marching cubes alghorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateReconstructedMeshArtifacts(Mesh* mesh, const AppOptions* options);

/**
 * @brief generateVolumeArtifacts
 * @param volume
 * @param options
 */
void generateVolumeArtifacts(Volume* volume, const AppOptions* options);

/**
 * @brief reconstructVolumeFromMesh
 * Reconstruc a volume from a given mesh.
 *
 * @param inputMesh
 * Given mesh
 * @param options
 * System options
 * @param releaseInputMesh
 * If this flag is set to true, the input mesh will be release in between to save memeory.
 * @return
 * A pointer to the reconstructed volume.
 */
Volume* reconstructVolumeFromMesh(Mesh* inputMesh, const AppOptions* options,
                                  const bool& releaseInputMesh = true);

/**
 * @brief reconstructMeshFromVolume
 * Reconstruct a mesh from a given volume.
 *
 * @param volume
 * Given volume
 * @param options
 * System options
 * @return
 * A pointer to the reconstructed mesh.
 */
Mesh* reconstructMeshFromVolume(Volume* volume, const AppOptions* options);

/**
 * @brief reconstructAdvancedMeshFromVolume
 * Reconstructs an advanced mesh from a given volume.
 *
 * @param volume
 * Input volume.
 * @param options
 * System options.
 * @return
 * A pointer to the advanced mesh reconstructed from the volume.
 */
AdvancedMesh* reconstructAdvancedMeshFromVolume(Volume* volume, const AppOptions* options);

}
