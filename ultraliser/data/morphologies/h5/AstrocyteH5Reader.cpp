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

#include "AstrocyteH5Reader.h"
#include <data/morphologies/EndfootPatch.hh>

namespace Ultraliser
{

AstrocyteH5Reader::AstrocyteH5Reader(const std::string &h5MorphologyFilePath)
    : _h5MorphologyFilePath(h5MorphologyFilePath)
{
    // Read the file
    _h5MorphologyFile = new H5::H5File(_h5MorphologyFilePath, H5F_ACC_RDONLY);

    // Read the astrocyte coordinates, or center
    _readCoordinates();

    // Read the samples, and DO NOT FORGET TO TRANSLATE THEM
    _readSamples();

    // Read the structure
    _readStructure();

    // Read the indices of the endfeet points
    _readEndfeetPointsIndices();

    // Read the endfeet points
    _readEndfeetPoints();

    // Read the indices of the triangles of the endfeet
    _readEndfeetPatchesIndices();

    // Read the endfeet triangles
    _readEndfeetPatches();

    // Construct the endfeet data
    _constructEndfeetData();
}

AstrocyteMorphology* AstrocyteH5Reader::getMorphology()
{
    AstrocyteMorphology* astrocyteMorphology = new AstrocyteMorphology(_skeletonSamples,
                                                                       _structure,
                                                                       _endfeetPatches,
                                                                       _coordinates);

    // Return the pointer to the vasculature morphology
    return astrocyteMorphology;
}

void AstrocyteH5Reader::_readSamples()
{
    // Read the points data set
    H5::DataSet pointsDataSet = _h5MorphologyFile->openDataSet("points");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = pointsDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _skeletonSamples.resize(dimenions[0]);

    // Read the data
    pointsDataSet.read(_skeletonSamples.data(), H5::PredType::NATIVE_FLOAT);

    // Close the dataset
    pointsDataSet.close();

    // Translate the samples to the global coordinates
    for (uint64_t i = 0; i < _skeletonSamples.size(); ++i)
    {
        _skeletonSamples[i].x += _coordinates.x();
        _skeletonSamples[i].y += _coordinates.y();
        _skeletonSamples[i].z += _coordinates.z();
    }

}

void AstrocyteH5Reader::_readStructure()
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

void AstrocyteH5Reader::_readCoordinates()
{
    // Read the coordinates dataset
    H5::DataSet coordinatesDataset = _h5MorphologyFile->openDataSet("coordinates");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = coordinatesDataset.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    std::vector< float > coordinates;
    coordinates.resize(dimenions[0]);

    // Read the data
    coordinatesDataset.read(coordinates.data(), H5::PredType::NATIVE_FLOAT);

    _coordinates.x() = coordinates[0];
    _coordinates.y() = coordinates[1];
    _coordinates.z() = coordinates[2];

    // Close the dataset
    coordinatesDataset.close();
}

void AstrocyteH5Reader::_readEndfeetPointsIndices()
{
    // Read the structure data set
    H5::DataSet endfeetPointsIndicesDataSet =
            _h5MorphologyFile->openDataSet("endfeet_vertex_indices");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = endfeetPointsIndicesDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _endfeetPointsIndices.resize(dimenions[0]);

    // Read the data
    endfeetPointsIndicesDataSet.read(_endfeetPointsIndices.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    endfeetPointsIndicesDataSet.close();
}

void AstrocyteH5Reader::_readEndfeetPoints()
{
    // Read the points data set
    H5::DataSet endfeetPointsDataSet = _h5MorphologyFile->openDataSet("endfeet_vertex_data");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = endfeetPointsDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _endfeetSamples.resize(dimenions[0]);

    // Read the data
    endfeetPointsDataSet.read(_endfeetSamples.data(), H5::PredType::NATIVE_FLOAT);

    // Close the dataset
    endfeetPointsDataSet.close();
}

void AstrocyteH5Reader::_readEndfeetPatchesIndices()
{
    // Read the structure data set
    H5::DataSet endfeetTrianglesIndicesDataSet =
            _h5MorphologyFile->openDataSet("endfeet_triangle_indices");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = endfeetTrianglesIndicesDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _endfeetTrianglesIndices.resize(dimenions[0]);

    // Read the data
    endfeetTrianglesIndicesDataSet.read(_endfeetTrianglesIndices.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    endfeetTrianglesIndicesDataSet.close();
}

void AstrocyteH5Reader::_readEndfeetPatches()
{
    // Read the structure data set
    H5::DataSet endfeetTrianglesDataSet = _h5MorphologyFile->openDataSet("endfeet_triangle_data");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataspace = endfeetTrianglesDataSet.getSpace();

    // Dataset dimenions
    hsize_t dimenions[2];
    dataspace.getSimpleExtentDims(dimenions, nullptr);

    // Resize the vector to contain the data
    _endfeetTriangles.resize(dimenions[0]);

    // Read the data
    endfeetTrianglesDataSet.read(_endfeetTriangles.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    endfeetTrianglesDataSet.close();
}

void AstrocyteH5Reader::_constructEndfeetData()
{
    for (uint64_t i = 0; i < _endfeetTrianglesIndices.size(); ++i)
    {
        // Get the triangle indices
        const auto firstTriangleIndex = _endfeetTrianglesIndices[i].firstIndex;
        const auto lastTriangleIndex = _endfeetTrianglesIndices[i].lastIndex;

        EndfootPatches endfootPatches;

        // Create the triangles
        for (uint64_t j = 0; j < _endfeetTriangles.size() -1 ; ++j)
        {
            // Get the incides of the vertices of the triangle
            uint64_t v0 = _endfeetTriangles[j].x();
            uint64_t v1 = _endfeetTriangles[j].y();
            uint64_t v2 = _endfeetTriangles[j].z();

            // Get the corresponding H5 sample
            H5Sample &s0 = _endfeetSamples.at(v0);
            H5Sample &s1 = _endfeetSamples.at(v1);
            H5Sample &s2 = _endfeetSamples.at(v2);

            // Convert it into a default sample, indicies values are not important in this context
            Sample* sample0 = new Sample(Vector3f(s0.x, s0.y, s0.z), s0.r * 0.5, 0);
            Sample* sample1 = new Sample(Vector3f(s1.x, s1.y, s1.z), s1.r * 0.5, 1);
            Sample* sample2 = new Sample(Vector3f(s2.x, s2.y, s2.z), s2.r * 0.5, 2);

            // Construct the endfeet patches
            EndfootPatch* patch = new EndfootPatch(sample0, sample1, sample2);

            // Append it to the patches list
            endfootPatches.push_back(patch);
        }

        _endfeetPatches.push_back(endfootPatches);
    }
}

}
