/***************************************************************************************************
 * Copyright (c) 2016 - 2021
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

#include "PBRT.h"
#include <utilities/File.h>

#define DOF_EXTENSION std::string(".dof")

namespace Ultraliser
{
namespace PBRT
{

PBRT::Resolution computeFilmResolution(PBRT::PROJECTION projection,
                                       uint64_t baseResolution,
                                       Vector3f pMin, Vector3f pMax)
{
    Resolution resolution;
    if (projection == PBRT::PROJECTION::XY)
    {
        float deltaX = pMax.x() - pMin.x();
        float deltaY = pMax.y() - pMin.y();

        if (deltaX > deltaY)
        {
            resolution.x = int(baseResolution);
            resolution.y = int(baseResolution * deltaY / deltaX);
        }
        else
        {
            resolution.x = int(baseResolution * deltaX / deltaY);
            resolution.y = int(baseResolution);
        }
    }
    else
    {
        float deltaZ = pMax.z() - pMin.z();
        float deltaY = pMax.y() - pMin.y();

        if (deltaZ > deltaY)
        {
            resolution.x = int(baseResolution);
            resolution.y = int(baseResolution * deltaY / deltaZ);
        }
        else
        {
            resolution.x = int(baseResolution * deltaZ / deltaY);
            resolution.y = int(baseResolution);
        }
    }

    return resolution;
}

void createDepthOfFieldFile(const std::string dofDirectory,
                            const Neuron& neuron,
                            Vector3f pMax,
                            Vector3f pMin)
{
    // Write DOF file
    std::stringstream dofFile;
    dofFile << dofDirectory << "/" << neuron.gid << DOF_EXTENSION;

    std::fstream dofStream;
    dofStream.open(dofFile.str().c_str(), std::ios::out);
    dofStream << "X+ve: " << pMax.x() - neuron.somaPosition.x() << std::endl;
    dofStream << "X-ve: " << neuron.somaPosition.x() - pMin.x() << std::endl;
    dofStream << "Z+ve: " << pMax.z() - neuron.somaPosition.z() << std::endl;
    dofStream << "Z-ve: " << neuron.somaPosition.z() - pMin.z() << std::endl;
    dofStream.close();
}

}
}
