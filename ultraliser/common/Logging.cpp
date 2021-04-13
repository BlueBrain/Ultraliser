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

#include "Logging.h"
#include "Logging.hh"
#include "Common.h"
#include "Colors.hh"

void log(const LOG_LEVEL& logLevel,
         const std::string &filePath,
         const int& lineNumber,
         const std::string &functionName,
         const char* lstring, ...)
{
    // Variable arguments information
    va_list argumentList;

    // Lof message
    char logMessage[1024];

    // Get the arguments and add them to the buffer
    va_start(argumentList, lstring);
    vsnprintf(logMessage, sizeof(logMessage), lstring, argumentList);

    // Get the time now
    time_t timeNow = time(nullptr);

    // Time stamp string
    char timeStamp[TIME_STAMP_CHAR_LENGTH];

    // Format the time string and get the local time
    strftime(timeStamp, sizeof(timeStamp),
              (char*)"%H:%M:%S", localtime(&timeNow));

    switch (logLevel)
    {
        case INFO:
#ifdef ULTRALISER_RELEASE
        printf(STD_WHITE "%s \n" STD_RESET,
               logMessage);
#else
        printf(STD_RED "[%d]" STD_RESET
               STD_YELLOW "[ %s ]" STD_RESET
               STD_CYAN " %s :" STD_RESET
               STD_GREEN "[%d]\n" STD_RESET
               STD_MAGENTA "\t* %s" STD_RESET
               STD_WHITE " %s \n" STD_RESET,
               (int) getpid(),
               timeStamp,
               filePath.c_str(),
               lineNumber,
               functionName.c_str(),
               logMessage);
#endif
        break;

    case DETAIL:
#ifdef ULTRALISER_RELEASE
    printf(STD_BOLD_CYAN "        * %s \n\n" STD_RESET,
           logMessage);
#else
    printf(STD_RED "[%d]" STD_RESET
           STD_YELLOW "[ %s ]" STD_RESET
           STD_CYAN " %s :" STD_RESET
           STD_GREEN "[%d]\n" STD_RESET
           STD_MAGENTA "\t* %s" STD_RESET
           STD_WHITE " %s \n" STD_RESET,
           (int) getpid(),
           timeStamp,
           filePath.c_str(),
           lineNumber,
           functionName.c_str(),
           logMessage);
#endif
    break;

    case STATUS:
#ifdef ULTRALISER_RELEASE
    printf(STD_WHITE "    %s \n" STD_RESET,
           logMessage);
#else
    printf(STD_RED "[%d]" STD_RESET
           STD_YELLOW "[ %s ]" STD_RESET
           STD_CYAN " %s :" STD_RESET
           STD_GREEN "[%d]\n" STD_RESET
           STD_MAGENTA "\t* %s" STD_RESET
           STD_WHITE " %s \n" STD_RESET,
           (int) getpid(),
           timeStamp,
           filePath.c_str(),
           lineNumber,
           functionName.c_str(),
           logMessage);
#endif
    break;

    case GREEN_STATUS:
#ifdef ULTRALISER_RELEASE
    printf(STD_GREEN "    %s \n" STD_RESET,
           logMessage);
#else
    printf(STD_RED "[%d]" STD_RESET
           STD_YELLOW "[ %s ]" STD_RESET
           STD_CYAN " %s :" STD_RESET
           STD_GREEN "[%d]\n" STD_RESET
           STD_MAGENTA "\t* %s" STD_RESET
           STD_WHITE " %s \n" STD_RESET,
           (int) getpid(),
           timeStamp,
           filePath.c_str(),
           lineNumber,
           functionName.c_str(),
           logMessage);
#endif
    break;

    case HEADER:
#ifdef ULTRALISER_RELEASE
    printf(STD_WHITE "%s \n" STD_RESET,
           logMessage);
#else
    printf(STD_RED "[%d]" STD_RESET
           STD_YELLOW "[ %s ]" STD_RESET
           STD_CYAN " %s :" STD_RESET
           STD_GREEN "[%d]\n" STD_RESET
           STD_MAGENTA "\t* %s" STD_RESET
           STD_WHITE " %s \n" STD_RESET,
           (int) getpid(),
           timeStamp,
           filePath.c_str(),
           lineNumber,
           functionName.c_str(),
           logMessage);
#endif
    break;

        case DEBUG:
#ifndef ULTRALISER_RELEASE
        printf(STD_RED "[%d]" STD_RESET
                STD_YELLOW "[ %s ]" STD_RESET
                STD_CYAN " %s :" STD_RESET
                STD_GREEN "[%d] \n" STD_RESET
                STD_MAGENTA "\t* %s" STD_RESET
                STD_BOLD_WHITE " %s \n" STD_RESET,
                (int) getpid(),
                timeStamp,
                filePath.c_str(),
                lineNumber,
                functionName.c_str(),
                logMessage);
#endif
            break;

        case WARNING:
#ifdef ULTRALISER_RELEASE
        printf(STD_BOLD_YELLOW "    WARNING: %s \n\n" STD_RESET,
               logMessage);
#else
        printf(STD_RED "[%d]" STD_RESET
               STD_YELLOW "[ %s ]" STD_RESET
               STD_CYAN " %s :" STD_RESET
               STD_GREEN "[%d] " STD_RESET
               STD_BOLD_YELLOW "[WARNGIN] \n" STD_RESET   // [WARNING]
               STD_MAGENTA "\t* %s" STD_RESET
               STD_BOLD_YELLOW " %s \n" STD_RESET,
               (int) getpid(),
               timeStamp,
               filePath.c_str(),
               lineNumber,
               functionName.c_str(),
               logMessage);
#endif
            break;

    case SUCCESS:
#ifdef ULTRALISER_RELEASE
    printf(STD_BOLD_GREEN "    NOTE: %s \n\n" STD_RESET,
           logMessage);
#else
    printf(STD_RED "[%d]" STD_RESET
           STD_YELLOW "[ %s ]" STD_RESET
           STD_CYAN " %s :" STD_RESET
           STD_GREEN "[%d] " STD_RESET
           STD_BOLD_YELLOW "[WARNGIN] \n" STD_RESET   // [WARNING]
           STD_MAGENTA "\t* %s" STD_RESET
           STD_BOLD_YELLOW " %s \n" STD_RESET,
           (int) getpid(),
           timeStamp,
           filePath.c_str(),
           lineNumber,
           functionName.c_str(),
           logMessage);
#endif
        break;

        case ERROR:
#ifdef ULTRALISER_RELEASE
        printf(STD_BOLD_RED "%s \n" STD_RESET,
               logMessage);
        exit(EXIT_SUCCESS);
#else
        printf(STD_RED "[%d]" STD_RESET
               STD_YELLOW "[ %s ]" STD_RESET
               STD_CYAN " %s :" STD_RESET
               STD_GREEN "[%d] " STD_RESET
               STD_BOLD_RED "[ERROR] \n" STD_RESET        // [ERROR]
               STD_MAGENTA "\t* %s" STD_RESET
               STD_BOLD_RED " %s \n" STD_RESET,
               (int) getpid(),
               timeStamp,
               filePath.c_str(),
               lineNumber,
               functionName.c_str(),
               logMessage);

#endif
            break;

        case EXIT:
        printf(STD_BOLD_RED
               "\t Exitting due to an error! ...\n " STD_RESET);
        exit(EXIT_FAILURE);
    }

    // Done
    va_end(argumentList);
}
