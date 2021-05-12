#include "MorphologyStatistics.h"

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
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < samples.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, samples.size());

        #pragma omp atomic
        ++progress;
#else
        LONG_LOOP_PROGRESS(progress, samples.size());
        ++progress;
#endif

        distribution[i] = samples[i]->getRadius();
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
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LONG_LOOP_PROGRESS(progress, sections.size());
        ++progress;
#endif

        distribution[i] = sections[i]->computeAverageRadius();
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
    uint64_t progress = 0;
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

std::vector< float > MorphologyStatistics::computeSectionsLengthDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Length Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LONG_LOOP_PROGRESS(progress, sections.size());
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

}

std::vector< float > MorphologyStatistics::computeSectionsSurfaceAreaDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Surface Area Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LONG_LOOP_PROGRESS(progress, sections.size());
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

}

std::vector< float > MorphologyStatistics::computeSectionsVolumeDistribution() const
{
    std::vector< float > distribution;
    const auto sections = _morphology->getSections();
    distribution.resize(sections.size());

    TIMER_SET;
    LOOP_STARTS("Computing Sections Volume Distribution");
    uint64_t progress = 0;
#ifdef ULTRALISER_USE_OPENMP
    #pragma omp parallel for
#endif
    for (uint64_t i = 0; i < sections.size(); ++i)
    {
#ifdef ULTRALISER_USE_OPENMP
        if (omp_get_thread_num() == 0)
            LOOP_PROGRESS(progress, sections.size());

        #pragma omp atomic
        ++progress;
#else
        LONG_LOOP_PROGRESS(progress, sections.size());
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
    File::writeFloatDistributionToFile(prefix + SAMPLES_RADII + DISTRIBUTION_EXTENSION,
                                       computeSamplesRadiiDistribution());

    File::writeFloatDistributionToFile(prefix + SECTION_AVERAGE_RADIUS + DISTRIBUTION_EXTENSION,
                                       computeSectionAverageRadiiDistribution());

}

}
