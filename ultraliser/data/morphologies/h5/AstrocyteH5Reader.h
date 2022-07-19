/***************************************************************************************************
 * Copyright (c) 2016 - 2022
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

#pragma once

#include <data/morphologies/h5/H5Sample.hh>
#include <data/morphologies/h5/H5Section.hh>
#include <data/morphologies/h5/EndfeetIndices.h>
#include <data/morphologies/AstrocyteMorphology.h>

namespace Ultraliser
{

/**
 * @brief The AstrocyteH5Reader class
 */
class AstrocyteH5Reader
{
public:

    /**
     * @brief AstrocyteH5Reader
     * Constructor.
     *
     * @param h5MorphologyFilePath
     * The path to the H5 morphology file.
     */
    AstrocyteH5Reader(const std::string &h5MorphologyFilePath);

public:

    /**
     * @brief getMorphology
     * Return a pointer to the vasculature morphology.
     * @return Return a pointer to the morphology.
     */
    AstrocyteMorphology* getMorphology();

private:

    /**
     * @brief _readSamples
     * Reads the samples from the morphology file.
     */
    void _readSamples();

    /**
     * @brief _readStructure
     * Reads the structure, or the sections, from the morphology file.
     */
    void _readStructure();

    /**
     * @brief _readCoordinates
     * Read the center of the morphology.
     */
    void _readCoordinates();

    /**
     * @brief _readEndfeetPointsIndices
     * Read the indices of the points of the endfeet.
     */
    void _readEndfeetPointsIndices();

    /**
     * @brief _readEndfeetPoints
     * Read the samples of the endfeet.
     */
    void _readEndfeetPoints();

    /**
     * @brief _readEndfeetPatchesIndices
     * Read the indices of the patches or triangles of the endfeet.
     */
    void _readEndfeetPatchesIndices();

    /**
     * @brief _readEndfeetPatches
     * Read the patches or the triangles of the endfeet.
     */
    void _readEndfeetPatches();

    /**
     * @brief _constructEndfeet
     * Construct the geometry of the endfeet from the patch data.
     */
    void _constructEndfeet();

private:

    /**
     * @brief _h5MorphologyFile
     * The path to the H5 morphology file
     */
    std::string _h5MorphologyFilePath;

    /**
     * @brief _h5MorphologyFile
     * A pointer to the .h5 file where we can access all of its data.
     */
    H5::H5File* _h5MorphologyFile;

    /**
     * @brief _samples
     * Vasculature samples.
     */
    H5Samples _skeletonSamples;

    /**
     * @brief _endfeetSamples
     * A list of the endfeet H5 samples.
     */
    H5Samples _endfeetSamples;

    /**
     * @brief _endfeetTriangles
     * A list of the endfeet triangles.
     */
    Triangles _endfeetTriangles;

    /**
     * @brief _endfeetPatches
     * A list of all the endfeet patches in the morphology.
     */
    EndfeetPatches _endfeetPatches;

    /**
     * @brief _structure
     * Endfeet structure, or sections.
     */
    H5Sections _structure;

    /**
     * @brief _coordinates
     * Astrocyte XYZ coordinates, or center.
     */
    Vector3f _coordinates;

    /**
     * @brief _endfeetPointsIndices
     * A list of the indices of the endfeet points.
     */
    EndfeetIndicesList _endfeetPointsIndices;

    /**
     * @brief _endfeetTrianglesIndices
     * A list of the indices of the endfeet triangles.
     */
    EndfeetIndicesList _endfeetTrianglesIndices;
};

}
#endif
