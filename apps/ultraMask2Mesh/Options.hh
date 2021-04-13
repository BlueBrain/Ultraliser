/*******************************************************************************
 * Copyright (c) 2016 - 2019
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 ******************************************************************************/

/*******************************************************************************
 * Copyright (c) 2016 - 2019
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License version 3.0 as published
 * by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this library; if not, write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 ******************************************************************************/

#ifndef OPTIONS_HH
#define OPTIONS_HH

#include <common/Common.h>

/**
 * @brief The Options struct
 */
struct Options
{
    /**
     * @brief maskDirectory
     * The input directory that contains the mask.
     */
    std::string maskDirectory;

    /**
     * @brief outputDirectory
     * The directory where the meshes (or data) will be created.
     * Output directory
     */
    std::string outputDirectory;

    /**
     * @brief maskWidth
     * The width of the mask.
     */
    int64_t maskWidth;

    /**
     * @brief maskHeight
     * The height of the mask.
     */
    int64_t maskHeight;

    /**
     * @brief volumeType
     * Use a specific volume for the voxelization process.
     */
    std::string volumeType;

    /**
     * @brief solid
     * Fill the interior of the volume using solid voxelization.
     */
    bool solid;

    /**
     * @brief projectXY
     * If this flag is set, the XY projection of the volume will be saved to a
     * PNG image. This flag is set to validate the output volume.
     */
    bool projectXY;

    /**
     * @brief projectZY
     * If this flag is set, the ZY projection of the volume will be saved to a
     * PNG image. This flag is set to validate the output volume.
     */
    bool projectZY;

    /**
     * @brief stackXY
     * Create an image stack along the XY plane.
     */
    bool stackXY;

    /**
     * @brief stackZY
     * Create an image stack along the ZY plane.
     */
    bool stackZY;

    /**
     * @brief writeBitVolume
     * If this flag is set, a binary volume will be created. This volume has
     * 1 bit per voxel.
     */
    bool writeBitVolume;

    /**
     * @brief writeByteVolume
     * If this flag is set, a default raw volume will be created.This volume
     * has 1 byte per voxel.
     */
    bool writeByteVolume;

    /**
     * @brief optimizeMesh
     * If this flag is set, the generated mesh will get optimised.
     */
    bool optimizeMesh;

    /**
     * @brief exportOFF
     * Export the reconstructed mesh as .OFF file.
     */
    bool exportOFF;

    /**
     * @brief exportOBJ
     * Export the reconstructed mesh as .OBJ file.
     */
    bool exportOBJ;

    /**
     * @brief xScale
     * X scaling factor for the mesh, default 1.0.
     */
    float xScale;

    /**
     * @brief yScale
     * Y scaling factor for the mesh, default 1.0.
     */
    float yScale;

    /**
     * @brief zScale
     * Z scaling factor for the mesh, default 1.0.
     */
    float zScale;

    /**
     * @brief smoothingIterations
     * Number of iterations required to smooth the optimized mesh, by default 10.
     */
    int64_t smoothingIterations;

    /**
     * @brief smoothingFactor
     * The rate at which the unnecessary vertices will be removed from the
     * optimized mesh. This fator has impact on the mesh size. The higher this
     * factor is the lower the mesh size becomes. By default 10.
     */
    float smoothingFactor;

    /**
     * @brief prefix
     * Just a prefix that will be used to label the output files. If this
     * is not given by the user, the name of the mesh file will be used.
     */
    std::string prefix;

    /**
     * @brief outputPrefix
     * Simply, the [OUTPUT_DIRECTORY]/[PREFIX]. This variable is just added to
     * make the code simpler.
     */
    std::string outputPrefix;
};

#endif // OPTIONS_HH
