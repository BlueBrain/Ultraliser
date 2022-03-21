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

#ifndef VASCULATURE_VMV_READER_H
#define VASCULATURE_VMV_READER_H

#include <data/morphologies/h5/H5Sample.hh>
#include <data/morphologies/h5/VasculatureH5Section.hh>
#include <data/morphologies/h5/VasculatureH5Connectivity.hh>
#include <data/morphologies/VasculatureMorphology.h>

namespace Ultraliser
{

class VasculatureVMVReader
{
public:

    /**
     * @brief VasculatureVMVReader
     * Constructor
     * @param vmvMorphologyFilePath
     * The path to the VMV morphology file.
     */
    VasculatureVMVReader(const std::string &vmvMorphologyFilePath);

public:

    /**
     * @brief getMorphology
     * Return a pointer to the vasculature morphology.
     * @return Return a pointer to the morphology.
     */
    VasculatureMorphology* getMorphology();

private:

    /**
     * @brief _readAttributes
     * Reads the atttributes from the morphology file to be able to process the file easily.
     */
    void _readAttributes();

    /**
     * @brief _readSamples
     * Reads the samples from the morphology file.
     */
    void _readSamples();

    /**
     * @brief _readConnectivity
     * Reads the connectivity information form the morphology file.
     */
    void _readConnectivity();

private:

    /**
     * @brief _numberVerts
     * Number of samples in the morphology.
     */
    uint64_t _numberVerts;

    /**
     * @brief _numberStrands
     * Number of strands or edges in the morphology.
     */
    uint64_t _numberStrands;

    /**
     * @brief _numberAttributesPerVertex
     * Number of attributes per vertex.
     */
    uint64_t _numberAttributesPerVertex;

    /**
     * @brief _vmvMorphologyFile
     * The path to the VMV morphology file.
     */
    std::string _vmvMorphologyFile;

    /**
     * @brief _samples
     * Vasculature samples.
     */
    Samples _samples;

    /**
     * @brief _sections
     * Vasculature sections or strands.
     */
    Sections _sections;
};

}

#endif // VASCULATURE_VMV_READER_H
