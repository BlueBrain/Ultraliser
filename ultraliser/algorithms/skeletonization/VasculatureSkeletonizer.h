/***************************************************************************************************
 * Copyright (c) 2016 - 2023
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include <algorithms/skeletonization/Skeletonizer.h>

namespace Ultraliser
{

/**
 * @brief The VasculatureSkeletonizer class
 */
class VasculatureSkeletonizer : public Skeletonizer
{
public:

    /**
     * @brief VasculatureSkeletonizer
     * Constructor
     * @param mesh
     * Input mesh
     * @param volume
     */
    VasculatureSkeletonizer(Volume *volume, const Mesh *mesh);
    ~VasculatureSkeletonizer();

    /**
     * @brief segmentComponents
     * Segment the different components of the graph after having it reconstructed from the
     * voxel grid.
     * NOTE: Currently, vasculature datasets have only branches with fluctuating diameters.
     * Therefore, branches are only the available components to be segmented from the graph.
     * This is indeed unlike neurons that can have somata, branches and spines.
     */
    void segmentComponents() override;

    /**
     * @brief exportSkeletonVMV
     * Exports the morphology skeleton into the VessMorphoVis (or VMV) file format.
     * For further details on this file format, please visit the following Wiki page:
     * https://github.com/BlueBrain/VessMorphoVis/wiki/File-Formats
     * @param prefix
     * The file prefix.
     * @param fileName
     * The file name.
     */
    void exportSkeletonVMV(const std::string& prefix, const std::string& fileName);
};

}
