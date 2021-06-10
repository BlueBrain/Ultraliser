/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 **************************************************************************************************/

#include <Ultraliser.h>

namespace Ultraliser
{


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
                                 Vector3f& pMax, Vector3f& pMin);
/**
 * @brief applyLaplacianOperator
 * Apply the laplacian operator on the mesh.
 *
 * @param mesh
 * An input mesh that is guaranteed to be two-manifold.
 * @param options
 * User-defined options given to the executable.
 */
void applyLaplacianOperator(Mesh *mesh, const Options* options);

/**
 * @brief createMeshWithNoSelfIntersections
 * Process the two manifold mesh and create a watertight mesh using the AdvancedMesh processing.
 *
 * @param manifoldMesh
 * An input mesh that is guaranteed to be two-manifold but have some self intersections.
 * @param options
 * User-defined options given to the executable.
 */
void createWatertightMesh(const Mesh* manifoldMesh, const Options* options);

/**
 * @brief optimizeMesh
 * Optimize the resulting mesh from the DMC algorithm.
 *
 * @param dmcMesh
 * The resulting mesh from the DMC algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateOptimizedMesh(Mesh *dmcMesh, const Options* options);

/**
 * @brief generateReconstructedMeshArtifacts
 * Writes the resulting mesh from the marching cubes algorithm.
 *
 * @param mesh
 * THe generetd mesh from the marching cubes algorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateReconstructedMeshArtifacts(const Mesh *mesh, const Options* options);

/**
 * @brief generateMeshArtifacts
 * Writes the resulting mesh artifacts.
 *
 * @param mesh
 * Reconstructed mesh from the marching cubes alghorithm.
 * @param options
 * User-defined options given to the executable.
 */
void generateMeshArtifacts(Mesh* mesh, const Options* options);

/**
 * @brief generateVolumeArtifacts
 * @param volume
 * @param options
 */
void generateVolumeArtifacts(const Volume* volume, const Options* options);
}
