/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s): Marwan Abdellah <marwan.abdellah@epfl.ch>
 *
 * This file is part of Ultraliser <https://github.com/BlueBrain/Ultraliser>
 *
 * This library is free software; you can redistribute it and/or modify it under the terms of the
 * GNU Lesser General Public License version 3.0 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along with this library;
 * if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 **************************************************************************************************/

#ifndef APPARGUMENTS_H
#define APPARGUMENTS_H

#include <Ultraliser.h>

namespace Ultraliser
{

class AppArguments
{
public:

    /**
     * @brief AppArguments
     * @param argc
     * @param argv
     * @param help
     */
    AppArguments(const int& argc , const char** argv, const std::string &help);

public:

    /**
     * @brief addInputMeshArguments
     * @param args
     */
    void addInputMeshArguments();

    /**
     * @brief addInputMeshesDirectoryArguments
     */
    void addInputMeshesDirectoryArguments();

    /**
     * @brief addInputMorphologyArguments
     */
    void addInputMorphologyArguments();

    /**
     * @brief addInputMaskDirectoryArguments
     */
    void addInputMaskDirectoryArguments();

    /**
     * @brief addMaskArguments
     */
    void addMaskArguments();

    /**
     * @brief addInputVolumeArguments
     */
    void addInputVolumeArguments();

    /**
     * @brief addInputVolumeParametersArguments
     */
    void addInputVolumeParametersArguments();

    /**
     * @brief addOutputArguments
     * Apps output arguments.
     * @param args
     * Arguments structure.
     */
    void addOutputArguments();

    /**
     * @brief addVoxelizationArguments
     * @param args
     */
    void addVoxelizationArguments();

    /**
     * @brief addVolumeProjectionArguments
     * @param args
     */
    void addVolumeProjectionArguments();

    /**
     * @brief addVolumeExportArguments
     * @param args
     */
    void addVolumeExportArguments();


    /**
     * @brief addStacksArguments
     * @param args
     */
    void addStacksArguments();

    /**
     * @brief addMeshVoxelizationArgument
     */
    void addMeshVoxelizationArgument();
    /**
     * @brief addVolumeArguments
     * @param args
     */
    void addVolumeArguments();

    /**
     * @brief addMeshOptimizationArguments
     * @param args
     */
    void addMeshOptimizationArguments();

    /**
     * @brief addMeshExportArguments
     * @param args
     */
    void addMeshExportArguments();

    /**
     * @brief addLaplacianOperatorArguments
     * @param args
     */
    void addLaplacianOperatorArguments();

    /**
     * @brief addMeshScaleArguments
     */
    void addMeshScaleArguments();

    /**
     * @brief addMeshArguments
     * @param args
     */
    void addMeshArguments();

    /**
     * @brief addSuppressionArguments
     * @param args
     */
    void addSuppressionArguments();

    /**
     * @brief addDataArguments
     * @param args
     */
    void addDataArguments();

    /**
     * @brief getOptions
     * @return
     */
    Options* getOptions();

private:

    /**
     * @brief _args
     * Command line arguments handler.
     */
    Args* _args;

    /**
     * @brief _options
     * Systems options.
     */
    Options* _options;
};

}

#endif // APPARGUMENTS_H