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

#include "Directory.h"
#include "String.h"
#include <common/Common.h>

namespace Ultraliser
{
namespace Directory
{

std::string getName(const std::string &path)
{
    // Split the path '/'
    std::vector< std::string > tokens = String::split(path, '/');

    // Return the lase
    return tokens.at(tokens.size() - 1);
}

void list(const std::string &directory,
                    std::vector< std::string >& files,
                    std::string extension)
{
    DIR* descriptor;
    struct dirent* dirNet;

    if ((descriptor = opendir(directory.c_str())) == nullptr)
    {
        LOG_ERROR("Canno't open the directory");
    }

    while ((dirNet = readdir(descriptor)) != nullptr)
    {
        std::string file = std::string(dirNet->d_name);
        if (file.find(extension) != std::string::npos)
            files.push_back(file);
    }

    closedir(descriptor);
}

void locateMeshFiles(const std::string &directory,
                     std::vector< std::string >& meshFiles)
{
    // List all the .obj files
    list(directory, meshFiles, ".obj");

    // List all the .ply files
    list(directory, meshFiles, ".ply");
}

bool exists(std::string path)
{
    DIR* directory = opendir(path.c_str());
    if (directory)
    {
        closedir(directory);
        return true;
    }

    return false;
}

}
}
