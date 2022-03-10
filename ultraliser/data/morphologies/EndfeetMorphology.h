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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_ENDFEET_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_ENDFEET_MORPHOLOGY_H

#include <data/morphologies/Morphology.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

class EndfeetMorphology : public Morphology
{
public:

    /**
     * @brief EndfeetMorphology
     * Constructor
     *
     * @param samples
     * A list of samples.
     * @param sampleTriangles
     * A list of sample triangles.
     */
    EndfeetMorphology(Samples& samples, SampleTriangles& sampleTriangles);

public:


    /**
     * @brief getSampleTriangles
     * @return Returns a reference to the list of sample triangles.
     */
    SampleTriangles getSampleTriangles() const;

private:

    /**
     * @brief _samplesTriangles
     * A list of the actual sample triangles of the morphology.
     */
    SampleTriangles _sampleTriangles;
};

/**
 * @brief readEndfeetMorphology
 * Reads the endfeet morphology from a given file.
 * @param morphologyPath
 * The path to the morphology.
 * @return
 * EndfeetMorphology object.
 */
EndfeetMorphology* readEndfeetMorphology(std::string& morphologyPath);

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_NEURON_MORPHOLOGY_H
