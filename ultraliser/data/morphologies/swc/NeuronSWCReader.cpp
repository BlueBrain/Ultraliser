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

#include "NeuronSWCReader.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace Ultraliser
{

NeuronSWCReader::NeuronSWCReader(const std::string &swcMorphologyFilePath)
{
    // Read the samples
    _readSamples(swcMorphologyFilePath);
}

NeuronSWCReader::~NeuronSWCReader()
{
    for (auto sample: _samples)
        delete sample;
    _samples.clear();
}

NeuronMorphology* NeuronSWCReader::getMorphology()
{
    NeuronMorphology* neuronMorphology = new NeuronMorphology(_samples);

    // Return the pointer to the neuron morphology
    return neuronMorphology;
}

AstrocyteMorphology* NeuronSWCReader::getAstrocyteMorphology()
{
    AstrocyteMorphology* astrocyteMorphology = new AstrocyteMorphology(_samples);

    // Return the pointer to the astrocyte morphology
    return astrocyteMorphology;
}


void NeuronSWCReader::_readSamples(const std::string &swcMorphologyFilePath)
{
    _samples.clear();
    std::unordered_map< size_t, NeuronSWCSample* > samplesMap;
    try
    {
        std::ifstream swcFile(swcMorphologyFilePath.c_str());
        if (!swcFile.is_open())
        {
            LOG_WARNING("Unable to open an SWC morphology file [ %s ]",
                      swcMorphologyFilePath.c_str());
        }
        std::string line;
        while(!swcFile.eof())
        {
            std::getline(swcFile, line);

            // Remove whitespaces at the beginning of the line 
            line.erase(0, line.find_first_not_of(" "));
            if (line.empty() || line[0] == '#')
                continue;

            NeuronSWCSample* sample = new NeuronSWCSample();
            int64_t parentId, type;
            std::stringstream sstr(line);
            sstr >> sample->id >> type >>
                    sample->x >> sample->y >> sample->z >> sample->r >> parentId;

            switch(type)
            {
            case 1:
                sample->type = SWCSampleType::SOMA;
                break; 
            case 2:
                sample->type = SWCSampleType::AXON;
                break; 
            case 3:
                sample->type = SWCSampleType::BASAL;
                break; 
            case 4:
                sample->type = SWCSampleType::APICAL;
                break; 
            default:
                sample->type = UNKNOWN_SAMPLE;
                break;
            } 
            
            _samples.push_back(sample);

            samplesMap[sample->id] = sample;
            if (parentId >= 0)
            {
                auto search = samplesMap.find(parentId);
                if (search != samplesMap.end())
                {
                    auto parent = search->second;
                    parent->childrenSamples.push_back(sample);
                }
            }
            
        }
    }
    catch(const std::exception&)
    {
        LOG_WARNING("Unable to load an SWC morphology file [ %s ]", swcMorphologyFilePath.c_str());
    }
}

}
