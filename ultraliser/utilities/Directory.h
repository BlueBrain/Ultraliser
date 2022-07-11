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

#pragma once

#include <common/Common.h>
#include <arguments/Arguments.h>

namespace Ultraliser
{
namespace Directory
{

/**
 * @brief getDirectoryName
 * Gets the name of a given directory.
 *
 * @param path
 * The path of the directory.
 * @return
 * Returns the name of the directory.
 */
std::string getName(const std::string &path);

/**
 * @brief listDirectory
 * Lists all the contents of the directory.
 *
 * @param directory
 * The path of the directory.
 * @param files
 * A list of all the found files.
 * @param extension
 * The extension of the files.
 * @return
 */
void list(const std::string &directory, std::vector< std::string >& files, std::string extension);

/**
 * @brief locateMeshFiles
 * Locates all the meshes of supported extenstions and returns in _meshFiles_ all the found images.
 *
 * @param directory
 * Input directory where the meshes will be searched.
 * @param meshFiles
 * The output list where the meshes will be returned.
 */
void locateMeshFiles(const std::string &directory, std::vector< std::string >& meshFiles);

/**
 * @brief directoryExists
 * Checks if the directory exists or not.
 *
 * @param path
 * The path of the directory.
 * @return
 * True or false.
 */
bool exists(std::string path);

}
}
