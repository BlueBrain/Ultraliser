/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#include "MeshStatistics.h"
#include "MeshOperations.h"
#include "TriangleOperations.h"
#include <utilities/Utilities.h>
#include <math/Math.h>
#include <geometry/Intersection.h>
#include <common/Common.h>

namespace Ultraliser
{

void MeshStatistics::computeSurfaceAreaDistribution()
{
    // Resize
    _surfaceAreaDistribution.resize(_numberTriangles);
    uint64_t validTriangleCount = 0;

    TIMER_SET;
    LOOP_STARTS("Computing Surface Area Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        const float area = computeTriangleSurfaceArea(_vertices[triangle[0]],
                                                      _vertices[triangle[1]],
                                                      _vertices[triangle[2]]);

        if (area > 0.f)
        {
#ifdef ULTRALISER_USE_OPENMP

            #pragma omp atomic
            validTriangleCount++;

            #pragma omp atomic
            _averageSurfaceArea += area;
#else
            validTriangleCount++;
            _averageSurfaceArea += area;
#endif

            _surfaceAreaDistribution[i] = area;
        }
        else
        {
            _surfaceAreaDistribution[i] = -1.0;
        }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Normalize the average surface area
    _averageSurfaceArea /= validTriangleCount;
}
std::vector< float > MeshStatistics::computeTriangleAspectRatioDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Aspect Ratio Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleAspectRatio(_vertices[triangle[0]],
                                                         _vertices[triangle[1]],
                                                         _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float >
MeshStatistics::computeTriangleRadiusRatioDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Radius Ratio Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif

        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleRadiusRatio(_vertices[triangle[0]],
                                                         _vertices[triangle[1]],
                                                         _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleEdgeRatioDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Edge Ratio Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif

        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleEdgeRatio(_vertices[triangle[0]],
                                                       _vertices[triangle[1]],
                                                       _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleRadiusToEdgeRatioDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Radius to Edge Ratio Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
        Triangle triangle = _triangles[i];

        const float radiusRatio = computeTriangleRadiusRatio(
                    _vertices[triangle[0]], _vertices[triangle[1]],
                    _vertices[triangle[2]]);

        const float edgeRatio = computeTriangleEdgeRatio(
                    _vertices[triangle[0]], _vertices[triangle[1]],
                    _vertices[triangle[2]]);

        distribution[i] = radiusRatio / edgeRatio;

    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleMinAngleDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Min Angle Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif

        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleMinAngle(_vertices[triangle[0]],
                                                      _vertices[triangle[1]],
                                                      _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleMaxAngleDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Maximum Angle Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif

        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleMaxAngle(_vertices[triangle[0]],
                                                      _vertices[triangle[1]],
                                                      _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleShapeDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Shape Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleShape(_vertices[triangle[0]],
                    _vertices[triangle[1]], _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleShapeAndSizeDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Shape & Size Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleShapeAndSize(_vertices[triangle[0]],
                                                          _vertices[triangle[1]],
                                                          _vertices[triangle[2]],
                                                          _averageSurfaceArea);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleScaledJacobianDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Scaled Jacobian Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleScaledJacobian(_vertices[triangle[0]],
                _vertices[triangle[1]], _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleConditionNumberDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Condition Number Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleConditionNumber(_vertices[triangle[0]],
                _vertices[triangle[1]], _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleDistortionDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Distortion Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleDistrotion(_vertices[triangle[0]],
                _vertices[triangle[1]], _vertices[triangle[2]]);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MeshStatistics::computeTriangleRelativeSizeSquaredDistribution() const
{
    std::vector< float > distribution;
    distribution.resize(_numberTriangles);

    TIMER_SET;
    LOOP_STARTS("Computing Triangle Relative Size Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < _numberTriangles; ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, _numberTriangles);

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, _numberTriangles);
        ++progress;
#endif
        Triangle triangle = _triangles[i];
        if (_surfaceAreaDistribution[i] > 0.f)
            distribution[i] = computeTriangleRelativeSizeSquared(_vertices[triangle[0]],
                    _vertices[triangle[1]], _vertices[triangle[2]],
                    _averageSurfaceArea);
        else
            distribution[i] = -1.f;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

void MeshStatistics::writeStatsDistributions(const std::string &prefix)
{
    // Compute the surface area distribution
    computeSurfaceAreaDistribution();

    File::writeFloatDistributionToFile(
                prefix + RADIUS_RATIO_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleRadiusRatioDistribution());

    File::writeFloatDistributionToFile(
                prefix + ASPECT_RATIO_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleAspectRatioDistribution());

    File::writeFloatDistributionToFile(
                prefix + EDGE_RATIO_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleEdgeRatioDistribution());

    File::writeFloatDistributionToFile(
                prefix + RADIUS_TO_EDGE_RATIO_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleRadiusToEdgeRatioDistribution());

    File::writeFloatDistributionToFile(
                prefix + MIN_ANGLE_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleMinAngleDistribution());

    File::writeFloatDistributionToFile(
                prefix + MAX_ANGLE_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleMaxAngleDistribution());

    File::writeFloatDistributionToFile(
                prefix + TRIANGLE_SHAPE_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleShapeDistribution());

    File::writeFloatDistributionToFile(
                prefix + TRIANGLE_SHAPE_SIZE_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleShapeAndSizeDistribution());

    File::writeFloatDistributionToFile(
                prefix + SCALED_JACOBIAN_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleScaledJacobianDistribution());

    File::writeFloatDistributionToFile(
                prefix + CONDITION_NUMBER_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleConditionNumberDistribution());

    File::writeFloatDistributionToFile(
                prefix + DISTORTION_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleDistortionDistribution());

    File::writeFloatDistributionToFile(
                prefix + RELATIVE_SIZE_SUFFIX + DISTRIBUTION_EXTENSION,
                computeTriangleRelativeSizeSquaredDistribution());
}

}
