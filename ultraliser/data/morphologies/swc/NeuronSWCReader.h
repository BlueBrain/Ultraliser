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

#pragma once

#include <string>
#include <vector>
#include <data/morphologies/swc/NeuronSWCSample.hh>
#include <data/morphologies/NeuronMorphology.h>
#include <data/morphologies/AstrocyteMorphology.h>


namespace Ultraliser
{

class NeuronSWCReader
{
public:

    /**
     * @brief NeuronSWCReader
     * Constructor
     * @param swcMorphologyFilePath
     * The path to the morphology file.
     */
    NeuronSWCReader(const std::string &swcMorphologyFilePath);
    
    ~NeuronSWCReader();
    
public:

    /**
     * @brief getMorphology
     * Return a pointer to the neuron morphology.
     * @return Return a pointer to the morphology.
     */
    NeuronMorphology* getMorphology();

    /**
     * @brief getAstrocyteMorphology
     * @return
     */
    AstrocyteMorphology* getAstrocyteMorphology();

private:

    /**
     * @brief _readSamples
     * Reads the samples from the morphology file.
     */
    void _readSamples(const std::string &swcMorphologyFilePath);

private:

    /**
     * @brief _samples
     * Neuron samples.
     */
    NeuronSWCSamples _samples;
};

}
