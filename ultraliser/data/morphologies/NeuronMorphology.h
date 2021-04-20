/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_NEURON_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_NEURON_MORPHOLOGY_H

#include <data/morphologies/Morphology.h>
#include <data/morphologies/NeuronSWCSample.hh>

namespace Ultraliser
{

class NeuronMorphology: public Morphology
{
public:

    /**
     * @brief NeuronMorphology
     * Constructor
     * @param swcSamples
     * Samples list
     */
    NeuronMorphology(const NeuronSWCSamples& swcSamples);

public:

    /**
     * @brief getSomaSamples
     * @return Returns a reference to the list of soma samples.
     */
    Samples getSomaSamples() const;

    /**
     * @brief getFirstSections
     * @return Returns a reference to the list of first neurites sections in the morphology.
     */
    Sections getFirstSections() const;

private:

    /**
     * @brief _constructMorphology
     * Loads the morphology data from the param swcSamples.
     * @param swcSamples
     * Samples list
     */
    void _constructMorphology(const NeuronSWCSamples& swcSamples);

private:

    /**
     * @brief _somaSamples
     * A list of the actual samples of the morphology soma
     */
    Samples _somaSamples;

    /**
     * @brief _firstSecion
     * A list of the first neurites secitions of the morphology 
     */
    Sections _firstSections;
};

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_NEURON_MORPHOLOGY_H
