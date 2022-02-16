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

#ifndef ULTRALISER_UTILITIES_TIMER_H
#define ULTRALISER_UTILITIES_TIMER_H

#include <common/Common.h>
#include <chrono>

namespace Ultraliser
{
namespace Utilities
{

/**
 * @brief The Timer class
 */
class Timer
{
public:
    Timer();

    /**
     * @brief start
     * Start the timer.
     */
    void start();

    /**
     * @brief stop
     * Stops the timer.
     */
    void stop();

    /**
     * @brief elapsedTimeInSeconds
     * Get the elapsed time in seconds.
     * @return
     * Returns the elapsed time in seconds.
     */
    double elapsedTimeInSeconds();

    /**
     * @brief elapsedTimeInMilliSeconds
     * Get the elapsed time in milli-seconds.
     * @return
     * Returns the elapsed time in milli-seconds.
     */
    double elapsedTimeInMilliSeconds();

    /**
     * @brief elapsedTimeInMicroSeconds
     * Get the elapsed time in micro-seconds.
     * @return
     * Returns the elapsed time in micro-seconds.
     */
    double elapsedTimeInMicroSeconds();

private:

    /**
     * @brief _startingTime
     */
    std::chrono::time_point<std::chrono::high_resolution_clock> _startingTime;

    /**
     * @brief _endTime
     */
    std::chrono::time_point<std::chrono::high_resolution_clock> _endTime;

    /**
     * @brief _duration
     */
    double _duration;
};

}
}

// Sets a timer and starts it
#define TIMER_SET Ultraliser::Utilities::Timer timer; timer.start()

// Resets a previously declared timer
#define TIMER_RESET timer.start()

// Gets the time in seconds, TIMER_SET must be called before calling this
#define GET_TIME_SECONDS timer.elapsedTimeInSeconds()

// Gets the time in milli-seconds, TIMER_SET must be called before calling this
#define GET_TIME_MILLI_SECONDS timer.elapsedTimeInMilliSeconds()

// Gets the time in micro-seconds, TIMER_SET must be called before calling this
#define GET_TIME_MICRO_SECONDS timer.elapsedTimeInMicroSeconds()

#endif // ULTRALISER_UTILITIES_TIMER_H
