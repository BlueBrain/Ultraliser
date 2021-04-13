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

#include "String.h"

namespace Ultraliser
{
namespace String
{

void removeSubstring(std::string& mainString, const std::string subString)
{
    std::string::size_type n = subString.length();
    for(std::string::size_type i = mainString.find(subString);
        i != std::string::npos;
        i = mainString.find(subString))
        mainString.erase(i, n);
}

bool replaceSubstring(std::string& mainString,
                      const std::string &from,
                      const std::string &to)
{
    size_t start_pos = mainString.find(from);
    if (start_pos == std::string::npos)
        return false;
    mainString.replace(start_pos, from.length(), to);
    return true;
}

bool subStringFound(std::string& mainString,
                    const std::string &subString)
{
    size_t start_pos = mainString.find(subString);
    if (start_pos == std::string::npos)
        return false;
    return true;
}

std::vector< std::string > split(std::string string, char delimiter)
{
    std::vector< std::string > internal;
    std::stringstream stream(string);
    std::string token;

    while(getline(stream, token, delimiter))
    {
        internal.push_back(token);
    }

    return internal;
}

std::string formatStringToMultiLine(const std::string &input,
                                    const uint64_t limit,
                                    const bool& aligned)
{
    // Get a list of all the words in the help
    std::vector< std::string > words = String::split(input, ' ');

    std::string formated = "\t";
    uint64_t currentPosition = 0;
    for (uint64_t i = 0; i < words.size(); ++i)
    {
        std::string word = words.at(i);
        if (currentPosition + word.size() < limit)
        {
            formated += " " + word;
            currentPosition += word.size();
        }
        else
        {
            if (aligned)
            {
                formated += "\n\t " + word;
            }
            else
            {
                formated += "\n\t\t " + word;
            }

            currentPosition = 0;
        }

    }
    return formated;
}

bool sameString(char *a, char *b)
{
    int i = 0;
    while (a[i] != '\0' && b[i] != '\0')
    {
        if (tolower(a[i]) != tolower(b[i])) return false;
        i++;
    }

    return (a[i] == '\0' && b[i] == '\0');
}

bool sameString(std::string s1, std::string s2)
{
    const char* a = s1.c_str();
    const char* b = s2.c_str();

    int i = 0;
    while (a[i] != '\0' && b[i] != '\0')
    {
        if (tolower(a[i]) != tolower(b[i])) return false;
        i++;
    }

    return (a[i] == '\0' && b[i] == '\0');
}

}
}

