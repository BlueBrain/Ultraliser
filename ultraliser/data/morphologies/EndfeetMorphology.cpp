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

#include "EndfeetMorphology.h"
#include <utilities/String.h>

namespace Ultraliser
{

EndfeetMorphology::EndfeetMorphology(Samples& samples, SampleTriangles& sampleTriangles)
{
    _samples = samples;
    _sampleTriangles = sampleTriangles; 

    _pMin = Vector3f(std::numeric_limits<float>::max());
    _pMax = Vector3f(std::numeric_limits<float>::lowest());
    for (auto sample: _samples)
    {
        Vector3f pMaxSample = sample->getPosition() + Vector3f(sample->getRadius());
        Vector3f pMinSample = sample->getPosition() - Vector3f(sample->getRadius());

        if (pMaxSample.x() > _pMax.x()) _pMax.x() = pMaxSample.x();
        if (pMaxSample.y() > _pMax.y()) _pMax.y() = pMaxSample.y();
        if (pMaxSample.z() > _pMax.z()) _pMax.z() = pMaxSample.z();

        if (pMinSample.x() < _pMin.x()) _pMin.x() = pMinSample.x();
        if (pMinSample.y() < _pMin.y()) _pMin.y() = pMinSample.y();
        if (pMinSample.z() < _pMin.z()) _pMin.z() = pMinSample.z();
    }
}


SampleTriangles EndfeetMorphology::getSampleTriangles() const
{
    return _sampleTriangles;
}

EndfeetMorphology* readEndfeetMorphology(std::string& morphologyPath)
{
    Samples samples;
    SampleTriangles sampleTriangles;

    samples.push_back(new Sample(Vector3f(0.0f, 1.0f, 0.5f), 0.01f, 0));

    samples.push_back(new Sample(Vector3f(-0.577f, 0.0f, 0.0f), 0.05f, 0));

    samples.push_back(new Sample(Vector3f(0.577f, 0.0f, 0.0f), 0.1f, 0));

    samples.push_back(new Sample(Vector3f(0.0f, -1.0f, 0.5f), 0.01f, 0));


    sampleTriangles.push_back(new SampleTriangle(samples[0], samples[1], samples[2]));
    sampleTriangles.push_back(new SampleTriangle(samples[1], samples[3], samples[2]));


    auto morphology = new EndfeetMorphology(samples, sampleTriangles);
    return morphology;



    // if (String::subStringFound(morphologyPath, std::string(".swc")))
    // {
    //     // Read the file
    //     // auto reader = std::make_unique< NeuronSWCReader >(morphologyPath);

    //     // Get a pointer to the morphology to start using it
    //     // return reader->getMorphology();
    // }
    // else
    // {
    //     LOG_ERROR("Unrecognized morphology file format [ %s ]", morphologyPath.c_str());
    // }

    // To avoid any warning issues.
    return nullptr;
}

}