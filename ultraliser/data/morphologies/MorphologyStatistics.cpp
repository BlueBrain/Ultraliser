/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#include "MorphologyStatistics.h"
#include "MorphologyOpertions.h"

namespace Ultraliser
{

MorphologyStatistics::MorphologyStatistics(const Morphology* morphology)
    : _morphology(morphology)
{
    /// EMPTY CONSTRUCTOR
}

std::vector< float > MorphologyStatistics::computeSamplesRadiiDistribution() const
{
    std::vector< float > distribution;
    const auto samples = _morphology->getSamples();
    distribution.resize(samples.size());

    TIMER_SET;
    LOOP_STARTS("Computing Samples Radii Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < samples.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, samples.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, samples.size());
        ++progress;
#endif

        distribution[i] = samples[i]->getRadius();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< PerAxisAnalysisData* >
MorphologyStatistics::computeSampleRadiusDistributionAlongAxis() const
{
    std::vector< PerAxisAnalysisData* > distribution;
    const auto samples = _morphology->getSamples();
    distribution.resize(samples.size());

    TIMER_SET;
    LOOP_STARTS("Computing Samples Radii Along Axis Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < samples.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, samples.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, samples.size());
        ++progress;
#endif

        distribution[i] = new PerAxisAnalysisData(samples[i]->getRadius(), samples[i]->getPosition());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeSectionAverageRadiiDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Average Radii Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->computeAverageRadius();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeNumberSamplesPerSectionDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Number of Samples per Sections Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->getSamples().size();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeSegmentsLengthDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segments Length Distribution");
    size_t progress = 0;
    for (const auto section: sections)
    {
        LOOP_PROGRESS(progress++, sections.size());

        auto result = section->computeSegmentsLengthDistribution();
        distribution.insert(distribution.end(), result.begin(), result.end());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return  distribution;
}


std::vector< PerAxisAnalysisData* >
MorphologyStatistics::computeSegmentLengthDistributionAlongAxis() const
{
    std::vector< PerAxisAnalysisData* > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segment Length Along Axis Distribution");
    size_t progress = 0;

    for (size_t i = 0; i < sections.size(); ++i)
    {
        const auto samples = sections[i]->getSamples();

        for (size_t j = 0; j < samples.size() - 1; ++j)
        {
            const auto length = computeSegmentLength(samples[j], samples[j + 1]);
            const Vector3f position = 0.5 * (samples[j]->getPosition() + samples[j + 1]->getPosition());
            distribution.push_back(new PerAxisAnalysisData(length, position));
        }

        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector<PerAxisAnalysisData* >
MorphologyStatistics::computeSegmentSurfaceAreaDistributionAlongAxis() const
{
    std::vector< PerAxisAnalysisData* > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segment Surface Area Along Axis Distribution");
    size_t progress = 0;

    for (size_t i = 0; i < sections.size(); ++i)
    {
        const auto samples = sections[i]->getSamples();

        for (size_t j = 0; j < samples.size() - 1; ++j)
        {
            const auto area = computeSegmentSurfaceArea(samples[j], samples[j + 1]);
            const Vector3f position = 0.5 * (samples[j]->getPosition() + samples[j+ 1]->getPosition());
            distribution.push_back(new PerAxisAnalysisData(area, position));
        }

        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector<PerAxisAnalysisData* >
MorphologyStatistics::computeSegmentVolumeDistributionAlongAxis() const
{
    std::vector< PerAxisAnalysisData* > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segment Volume Along Axis Distribution");
    size_t progress = 0;

    for (size_t i = 0; i < sections.size(); ++i)
    {
        const auto samples = sections[i]->getSamples();

        for (size_t j = 0; j < samples.size() - 1; ++j)
        {
            const auto volume = computeSegmentVolume(samples[j], samples[j + 1]);
            const Vector3f position = 0.5 * (samples[j]->getPosition() + samples[j + 1]->getPosition());
            distribution.push_back(new PerAxisAnalysisData(volume, position));
        }

        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeSectionsLengthDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Length Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->computeLength();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeSegmentsSurfaceAreaDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segments Surface Area Distribution");
    size_t progress = 0;
    for (const auto section: sections)
    {
        LOOP_PROGRESS(progress++, sections.size());

        auto result = section->computeSegmentsSurfaceAreaDistribution();
        distribution.insert(distribution.end(), result.begin(), result.end());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return  distribution;
}

std::vector< float > MorphologyStatistics::computeSectionsSurfaceAreaDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Surface Area Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->computeSurfaceArea();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

std::vector< float > MorphologyStatistics::computeSegmentsVolumeDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();

    TIMER_SET;
    LOOP_STARTS("Computing Segments Volume Distribution");
    size_t progress = 0;
    for (const auto section: sections)
    {
        LOOP_PROGRESS(progress++, sections.size());

        auto result = section->computeSegmentsVolumeDistribution();
        distribution.insert(distribution.end(), result.begin(), result.end());
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return  distribution;
}

std::vector< float > MorphologyStatistics::computeSectionsVolumeDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Volume Distribution");
    size_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LOOP_PROGRESS_FRACTION(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->computeVolume();
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    return distribution;
}

void MorphologyStatistics::writeStatsDistributions(const std::string &prefix)
{
    File::writeFloatDistributionToFile(
                prefix + SAMPLES_RADII + DISTRIBUTION_EXTENSION,
                computeSamplesRadiiDistribution());

    File::writePerAxisFloatDistributionToFile(
                prefix + PER_AXIS_SAMPLES_RADII + DISTRIBUTION_EXTENSION,
                computeSampleRadiusDistributionAlongAxis());

    File::writePerAxisFloatDistributionToFile(
                prefix + PER_AXIS_SEGMENTS_LENGTH + DISTRIBUTION_EXTENSION,
                computeSegmentLengthDistributionAlongAxis());

    File::writePerAxisFloatDistributionToFile(
                prefix + PER_AXIS_SEGMENTS_SURFACE_AREA + DISTRIBUTION_EXTENSION,
                computeSegmentSurfaceAreaDistributionAlongAxis());

    File::writePerAxisFloatDistributionToFile(
                prefix + PER_AXIS_SEGMENTS_VOLUME + DISTRIBUTION_EXTENSION,
                computeSegmentVolumeDistributionAlongAxis());

    File::writeFloatDistributionToFile(
                prefix + SECTION_AVERAGE_RADIUS + DISTRIBUTION_EXTENSION,
                computeSectionAverageRadiiDistribution());

    File::writeFloatDistributionToFile(
                prefix + NUMBER_SAMPLES_PER_SECTION + DISTRIBUTION_EXTENSION,
                computeNumberSamplesPerSectionDistribution());

    File::writeFloatDistributionToFile(
                prefix + SEGMENTS_LENGTH + DISTRIBUTION_EXTENSION,
                computeSegmentsLengthDistribution());

    File::writeFloatDistributionToFile(
                prefix + SECTIONS_LENGTH + DISTRIBUTION_EXTENSION,
                computeSectionsLengthDistribution());

    File::writeFloatDistributionToFile(
                prefix + SEGMENTS_SURFACE_AREA + DISTRIBUTION_EXTENSION,
                computeSegmentsSurfaceAreaDistribution());

    File::writeFloatDistributionToFile(
                prefix + SECTIONS_SURFACE_AREA + DISTRIBUTION_EXTENSION,
                computeSectionsSurfaceAreaDistribution());

    File::writeFloatDistributionToFile(
                prefix + SEGMENTS_VOLUME + DISTRIBUTION_EXTENSION,
                computeSegmentsVolumeDistribution());

    File::writeFloatDistributionToFile(
                prefix + SECTIONS_VOLUME + DISTRIBUTION_EXTENSION,
                computeSectionsVolumeDistribution());
}

}
