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

#include "Timer.h"

namespace Ultraliser
{
namespace Utilities
{

Timer::Timer()
{
    _duration = 0.0;
}

void Timer::start()
{
    // Clock
    _startingTime = std::chrono::high_resolution_clock::now();
}

void Timer::stop()
{
    // Clock
    _endTime = std::chrono::high_resolution_clock::now();
}

double Timer::elapsedTimeInSeconds()
{
    // Clock
    _endTime = std::chrono::high_resolution_clock::now();

    // Get the duration
    _duration = std::chrono::duration< double, std::milli >
            (_endTime-_startingTime).count();

    return _duration * 1e-3;
}

double Timer::elapsedTimeInMilliSeconds()
{
    // Clock
    _endTime = std::chrono::high_resolution_clock::now();

    // Get the duration
    _duration = std::chrono::duration< double, std::milli >
            (_endTime-_startingTime).count();

    return _duration;
}

double Timer::elapsedTimeInMicroSeconds()
{
    // Clock
    _endTime = std::chrono::high_resolution_clock::now();

    // Get the duration
    _duration = std::chrono::duration< double, std::micro >
            (_endTime-_startingTime).count();

    return _duration;
}

}
}
