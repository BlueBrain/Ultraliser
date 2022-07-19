/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
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

#ifdef ULTRALISER_USE_H5

#include "VasculatureH5Reader.h"
#include <iostream>

namespace Ultraliser
{

VasculatureH5Reader::VasculatureH5Reader(const std::string &h5MorphologyFilePath)
    : _h5MorphologyFilePath(h5MorphologyFilePath)
{
    // Read the file
    _h5MorphologyFile = new H5::H5File(_h5MorphologyFilePath, H5F_ACC_RDONLY);

    // Read the samples
    _readSamples();

    // Read the structure
    _readStructure();

    // Read the connectivity
    _readConnectivity();
}

VasculatureMorphology* VasculatureH5Reader::getMorphology()
{
    VasculatureMorphology* vasculatureMorphology =
            new VasculatureMorphology(_samples, _structure, _connectivity);

    // Return the pointer to the vasculature morphology
    return vasculatureMorphology;
}

void VasculatureH5Reader::_readSamples()
{
    // Read the points data set
    H5::DataSet pointsDataSet = _h5MorphologyFile->openDataSet("points");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = pointsDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _samples.resize(dimenions[0]);

    // Read the data
    pointsDataSet.read(_samples.data(), H5::PredType::NATIVE_FLOAT);

    // Close the dataset
    pointsDataSet.close();
}

void VasculatureH5Reader::_readStructure()
{
    // Read the structure data set
    H5::DataSet structureDataSet = _h5MorphologyFile->openDataSet("structure");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = structureDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _structure.resize(dimenions[0]);

    // Read the data
    structureDataSet.read(_structure.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    structureDataSet.close();
}

void VasculatureH5Reader::_readConnectivity()
{
    // Read the connectivity data set
    H5::DataSet connectivityDataSet =
            _h5MorphologyFile->openDataSet("connectivity");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = connectivityDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _connectivity.resize(dimenions[0]);

    // Read the data
    connectivityDataSet.read(_connectivity.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    connectivityDataSet.close();
}
}

#endif

