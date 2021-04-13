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

#ifndef ULTRALISER_COMMON_MACROS_H
#define ULTRALISER_COMMON_MACROS_H

#include <common/Headers.hh>
#include <common/Logging.h>

// Setting a counter
#define LOOP_COUNTER_SET uint64_t COUNTER = 0
#define LOOP_COUNTER_RESET COUNTER = 0

#define LOOP_STARTS(MESSAGE)                                                                        \
{                                                                                                   \
    printf("\t%s \n", MESSAGE);                                                                     \
}

// Print the progress in a loop
#ifdef ENABLE_PROGRESS_BAR
#define LOOP_PROGRESS(PROGRESS, TOTAL)                                                              \
{                                                                                                   \
    {                                                                                               \
        double PERCENTAGE = static_cast< double >                                                   \
            (100.0 * PROGRESS / static_cast< double >(TOTAL));                                      \
        double STARS = std::floor((PERCENTAGE * PROGRESS_BAR_LENGTH) / 100.0);                      \
        double SPACES = std::floor(PROGRESS_BAR_LENGTH - STARS);                                    \
        std::string BAR = "* Progress │";                                                           \
        for(int i = 0; i < int(STARS); i++) BAR += "▒";                                             \
        for(int i = 0; i < int(SPACES); i++) BAR += " "; BAR+= "│";                                 \
        printf("\r\t%s (%2.2f %%)", BAR.c_str(), PERCENTAGE);                                       \
        fflush(stdout);                                                                             \
    }                                                                                               \
}
#else
#define LOOP_PROGRESS(PROGRESS, TOTAL) { }
#endif

// Print the progress in a very long loop
#define LONG_LOOP_PROGRESS(PROGRESS, TOTAL)                                                         \
{                                                                                                   \
    {                                                                                               \
        if (((PROGRESS + 1) % 1000) == 0) LOOP_PROGRESS(PROGRESS, TOTAL)                            \
    }                                                                                               \
}

// Print the status after the loop is done
#define LOOP_DONE                                                                                   \
{                                                                                                   \
    LOOP_PROGRESS(100, 100);                                                                        \
    printf(" \n");                                                                                  \
}

#endif // ULTRALISER_COMMON_MACROS_H
