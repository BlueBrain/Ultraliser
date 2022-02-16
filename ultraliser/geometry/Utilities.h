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

#ifndef ULTRALISER_GEOMRTY_UTILITIES_H
#define ULTRALISER_GEOMRTY_UTILITIES_H

#include <common/Common.h>
#include <math/Vector.h>
#include <data/NeuronData.h>


namespace Ultraliser
{
namespace Utilities
{

struct NeuronBoundingBox
{
    Vector3f pMin;
    Vector3f pMax;
};

/**
 * @brief ComputeMaximumBoundingBox
 * @param boundingBoxList
 * @param pMin
 * @param pMax
 */
void computeMaximumBoundingBox(std::vector< NeuronBoundingBox> boundingBoxList,
                               Vector3f& pMin, Vector3f& pMax);
/**
 * @brief inRange
 * @param value
 * @param lowerBound
 * @param upperBound
 * @return
 */
template < class T >
bool inRange(T value, T lowerBound, T upperBound);

/**
 * @brief BBox
 * @param mesh
 * @param pMin
 * @param pMax
 */
// void BBox(const OriginalMesh & mesh, Vector3f & pMin, Vector3f & pMax);

/**
 * @brief BBox
 * @param points
 * @param pMin
 * @param pMax
 * @param init
 */
void BBox(const std::vector< Vector3f > & points,
          Vector3f & pMin, Vector3f & pMax,
          const bool initialized = false);

/**
 * @brief pInside
 * @param pMin
 * @param pMax
 * @param p
 * @return
 */
bool pInside(const Vector3f& pMin, const Vector3f& pMax, const Vector3f& p);


/**
 * @brief isNbr
 * @param a
 * @param b
 * @param vert
 * @return
 */
bool isNbr(const Vec3i_64 &a, const Vec3i_64 &b, int64_t vert);

/**
 * @brief adjlist
 * @param mesh
 * @param adjMat
 */
//void adjustList(const OriginalMesh & mesh, std::vector<std::vector<int32_t> > &adjMat);

/**
 * @brief printBoundingBoxData
 * @param pMin
 * @param pMax
 * @return
 */
void printBoundingBoxData(const Vector3f &pMin, const Vector3f &pMax,
                          std::string prefix = "volume");

/**
 * @brief getBoundingBox
 * @param options
 * @param neurons
 * @param pMax
 * @param pMin
 */
//void getBoundingBox(const NeuroVoxyOptions* options,
//                    Neurons neurons, Vector3f& pMax, Vector3f& pMin);
/**
 * @brief printVoxelData
 * @param pMin
 * @param pMax
 * @param volumeSize
 * @param prefix
 * @return
 */
void printVoxelData(const Vector3f &pMin, const Vector3f &pMax,
                    const Vector3f volumeSize,
                    std::string prefix = "volume");

}
}


#endif // ULTRALISER_GEOMRTY_UTILITIES_H
