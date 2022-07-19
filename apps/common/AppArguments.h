/***************************************************************************************************
 * Copyright (c) 2016 - 2022
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

#include <Ultraliser.h>
#include <AppOptions.h>

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
     * @brief addMorphologyExtractionArguments
     */
    void addMorphologyExtractionArguments();

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
     */
    void addOutputArguments();

    /**
     * @brief addSolidVoxelizationArguments
     */
    void addSolidVoxelizationArguments();

    /**
     * @brief addVoxelizationArguments
     */
    void addVoxelizationArguments();

    /**
     * @brief addVolumeProjectionArguments
     */
    void addVolumeProjectionArguments();

    /**
     * @brief addVolumeExportArguments
     */
    void addVolumeExportArguments();

    /**
     * @brief addStacksArguments
     */
    void addStacksArguments();

    /**
     * @brief addMeshVoxelizationArgument
     */
    void addMeshVoxelizationArgument();
    /**
     * @brief addVolumeArguments
     */
    void addVolumeArguments();

    /**
     * @brief addMeshJoiningArguments
     */
    void addMeshJoiningArguments();

    /**
     * @brief addMeshExtractionArguments
     */
    void addMeshExtractionArguments();

    /**
     * @brief addMeshOptimizationArguments
     */
    void addMeshOptimizationArguments();

    /**
     * @brief addMeshExportArguments
     */
    void addMeshExportArguments();

    /**
     * @brief addLaplacianOperatorArguments
     */
    void addLaplacianOperatorArguments();

    /**
     * @brief addMeshScaleArguments
     */
    void addMeshScaleArguments();

    /**
     * @brief addMeshArguments
     */
    void addMeshArguments();

    /**
     * @brief addSuppressionArguments
     */
    void addSuppressionArguments();

    /**
     * @brief addDataArguments
     */
    void addDataArguments();

    /**
     * @brief addAstrocyteSpecificArguments
     */
    void addAstrocyteSpecificArguments();

    /**
     * @brief addPackingAlgorithmArguments
     */
    void addPackingAlgorithmArguments();

    /**
     * @brief addNeuronMorphologyBranchOrderArguments
     */
    void addNeuronMorphologyBranchOrderArguments();

    /**
     * @brief getOptions
     * @return System option
     */
    AppOptions* getOptions();

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
    AppOptions* _options;
};

}
