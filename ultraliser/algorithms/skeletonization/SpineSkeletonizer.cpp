/***************************************************************************************************
 * Copyright (c) 2016 - 2024
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

#include "SpineSkeletonizer.h"

namespace Ultraliser
{

SpineSkeletonizer::SpineSkeletonizer(Volume* spineVolume,
                                     const bool useAcceleration,
                                     const bool debugSkeleton,
                                     const std::string debuggingPrefix)
    : Skeletonizer(spineVolume, useAcceleration, debugSkeleton, debuggingPrefix)
{
    /// EMPTY CONSTRUCTOR
}

void SpineSkeletonizer::run(const bool verbose)
{
    // Initialize
    initialize(verbose);

    skeletonizeVolumeToCenterLines();

    /// Extract the nodes of the skeleton from the center-line "thinned" voxels and return a
    /// mapper that maps the indices of the voxels in the volume and the nodes in the skeleton
    auto indicesMapper = _extractNodesFromVoxels();

    /// Connect the nodes of the skeleton to construct its edges. This operation will not connect
    /// any gaps, it will just connect the nodes extracted from the voxels.
    _connectNodesToBuildEdges(indicesMapper);

    /// Inflate the nodes, i.e. adjust their radii
    _inflateNodes();

    /// Reconstruct the sections "or branches" from the nodes using the edges data
    _buildBranchesFromNodes(_nodes);


    _exportBranches(_debuggingPrefix);


    // Identify the connections at the terminals of each branch
    // identifyTerminalConnections(_branches);

    // Roots, terminals and others
    // confirmTerminalsBranches(_branches);
}

void SpineSkeletonizer::segmentComponents()
{

}


void SpineSkeletonizer::_exportBranches(const std::string& prefix)
{

    // Construct the file path
    std::string filePath = prefix + BRANCHES_EXTENSION;

    std::fstream stream;
    stream.open(filePath, std::ios::out);

    size_t progress = 0;
    for (size_t i = 0; i < _branches.size(); ++i)
    {
        // The @start marks a new branch in the file
        stream << "start " << _branches[i]->index << "\n";

        for (auto& node: _branches[i]->nodes)
        {
            stream << node->point.x() << " "
                   << node->point.y() << " "
                   << node->point.z() << " "
                   << node->radius << "\n";
        }
        // The @end marks the terminal sample of a branch
        stream << "end\n";

        ++progress;
    }

    // Close the file
    stream.close();
}
}
