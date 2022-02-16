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

#ifndef ULTRALISER_DATA_MESH_SIMPLE_MESH_STATISTICS_H
#define ULTRALISER_DATA_MESH_SIMPLE_MESH_STATISTICS_H

#include <common/Common.h>
#include <math/Math.h>
#include <data/meshes/simple/primitives/Primitives.h>
#include <data/meshes/simple/NeighbourTriangle.h>

namespace Ultraliser
{

/**
 * @brief The MeshStatistics class
 */
class MeshStatistics
{
public:

    /**
     * @brief MeshStatistics
     * @param vertices
     * @param triangles
     * @param numberVertices
     * @param numberTriangles
     */
    MeshStatistics(Vertex* vertices,
                   Triangle* triangles,
                   uint64_t numberVertices,
                   uint64_t numberTriangles)
        : _vertices(vertices)
        , _triangles(triangles)
        , _numberVertices(numberVertices)
        , _numberTriangles(numberTriangles)
    { }

    /**
     * @brief computeSurfaceAreaDistribution
     */
    void computeSurfaceAreaDistribution();

    /**
     * @brief computeTriangleAspectRatioDistribution
     * @return
     */
    std::vector< float > computeTriangleAspectRatioDistribution() const;

    /**
     * @brief computeTriangleRadiusRatioDistribution
     * @return
     */
    std::vector< float > computeTriangleRadiusRatioDistribution() const;

    /**
     * @brief computeTriangleEdgeRatioDistribution
     * @return
     */
    std::vector< float > computeTriangleEdgeRatioDistribution() const;

    /**
     * @brief computeTriangleRadiusToEdgeRatioDistribution
     * @return
     */
    std::vector< float > computeTriangleRadiusToEdgeRatioDistribution() const;

    /**
     * @brief computeTriangleMinAngleDistribution
     * @return
     */
    std::vector< float > computeTriangleMinAngleDistribution() const;

    /**
     * @brief computeTriangleMaxAngleDistribution
     * @return
     */
    std::vector< float > computeTriangleMaxAngleDistribution() const;

    /**
     * @brief computeTriangleShapeDistribution
     * @return
     */
    std::vector< float > computeTriangleShapeDistribution() const;

    /**
     * @brief computeTriangleShapeAndSizeDistribution
     * @return
     */
    std::vector< float > computeTriangleShapeAndSizeDistribution() const;

    /**
     * @brief computeTriangleRelativeSizeSquaredDistribution
     * @return
     */
    std::vector< float > computeTriangleRelativeSizeSquaredDistribution() const;

    /**
     * @brief computeTriangleConditionNumberDistribution
     * @return
     */
    std::vector< float > computeTriangleConditionNumberDistribution() const;

    /**
     * @brief computeTriangleDistortionDistribution
     * @return
     */
    std::vector< float > computeTriangleDistortionDistribution() const;

    /**
     * @brief computeTriangleAspectFrobeniusDistribution
     * @return
     */
    std::vector< float > computeTriangleAspectFrobeniusDistribution() const;

    /**
     * @brief computeTriangleScaledJacobianDistribution
     * @return
     */
    std::vector< float > computeTriangleScaledJacobianDistribution() const;

    /**
     * @brief writeStatsDistributions
     * @param reference
     */
    void writeStatsDistributions(const std::string &prefix);

private:

    /**
     * @brief _vertices
     */
    const Vertex* _vertices;

    /**
     * @brief _triangles
     */
    const Triangle* _triangles;

    /**
     * @brief _numberVertices
     */
    const uint64_t _numberVertices;

    /**
     * @brief _numberTriangles
     */
    const uint64_t _numberTriangles;

    /**
     * @brief _surfaceAreaDistribution
     */
    std::vector< float > _surfaceAreaDistribution;

    /**
     * @brief _averageSurfaceArea
     */
    float _averageSurfaceArea;
};

}

#endif // ULTRALISER_DATA_MESH_SIMPLE_MESH_STATISTICS_H
