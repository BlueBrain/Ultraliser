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

#ifndef ULTRALISER_COMMON_PROGRESS_H
#define ULTRALISER_COMMON_PROGRESS_H

#include <common/Headers.hh>
#include <common/Logging.h>

/**
 * @brief printProgressBar
 * Print the progress bar in a loop.
 *
 * @param current
 * Current count
 * @param total
 * Total count
 * @param barLength
 * The total langth of the progress bar.
 */
static inline void printProgressBar(
        const uint64_t& current, const uint64_t& total, const uint64_t barLength = 50)
{
    float percentage = ((100.f * current) / (1.f * total));
    float starts = std::floor((percentage * barLength) / 100.f);
    float spaces = std::floor(barLength - starts);
    std::string bar = "* Progress │";
    for(int i = 0; i < int(starts); i++) bar += "▒";
    for(int i = 0; i < int(spaces); i++) bar += " "; bar+= "│";
    printf("\r\t%s (%2.2f %%)", bar.c_str(), percentage);
    fflush(stdout);
}

/**
 * @brief printFractionProgressBar
 * Prints a very simple progress bar in a loop that only increments evey 10% of the loop.
 *
 * @param current
 * Current count
 * @param total
 * Total count
 * @param barLength
 * The total langth of the progress bar.
 */
static inline void printFractionProgressBar(
        const uint64_t& current, const uint64_t& total, const uint64_t barLength = 50)
{
    float percentage = ((100.f * current) / (1.f * total));

    if (static_cast< uint64_t >(percentage) % 10 == 0)
    {
        float starts = std::floor((percentage * barLength) / 100.f);
        float spaces = std::floor(barLength - starts);
        std::string bar = "* Progress │";
        for(int i = 0; i < int(starts); i++) bar += "▒";
        for(int i = 0; i < int(spaces); i++) bar += " "; bar+= "│";
        printf("\r\t%s (%2.2f %%)", bar.c_str(), percentage);
        fflush(stdout);
    }
}

/**
 * @brief progressUpdate
 * Update the progress of the current thread in a loop.
 *
 * @param progressValue
 * The progress value to be updated.
 */
static void inline progressUpdate(uint64_t& progressValue)
{
#ifdef ULTRALISER_USE_OPENMP
#pragma omp atomic
#endif
    ++progressValue;
}

// Setting a counter
#define LOOP_COUNTER_SET uint64_t COUNTER = 0
#define LOOP_COUNTER_RESET COUNTER = 0

// Prints a simple message before starting the loop
#define LOOP_STARTS(MESSAGE) printf("\t%s \n", MESSAGE);

// Print the progress in a loop
#define LOOP_PROGRESS(PROGRESS, TOTAL) printProgressBar(PROGRESS, TOTAL)

// Print the progress in a loop only in fractions
#define LOOP_PROGRESS_FRACTION(PROGRESS, TOTAL) printFractionProgressBar(PROGRESS, TOTAL)

// Print the status after the loop is done
#define LOOP_DONE { LOOP_PROGRESS(100, 100); printf(" \n"); }

// The progress variable itself
#define PROGRESS ULTRALISER_PROGRESS

// Set the progress to zero
#define PROGRESS_SET uint64_t ULTRALISER_PROGRESS = 0

// Reset the progress
#define PROGRESS_RESET ULTRALISER_PROGRESS = 0

// Update the progress bar
#define PROGRESS_UPDATE progressUpdate(ULTRALISER_PROGRESS)

#endif // ULTRALISER_COMMON_PROGRESS_H
