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

#pragma once

#include <data/morphologies/Morphology.h>
#include <data/morphologies/swc/NeuronSWCSample.hh>
#include <data/morphologies/EndfootPatch.hh>
#include <data/morphologies/h5/H5Section.hh>
#include <data/morphologies/h5/H5Sample.hh>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/**
 * @brief The AstrocyteMorphology class
 */
class AstrocyteMorphology : public Morphology
{
public:

    /**
     * @brief AstrocyteMorphology
     * @param swcSamples
     */
    AstrocyteMorphology(const NeuronSWCSamples& swcSamples);

    /**
     * @brief EndfeetMorphology
     * Constructor
     *
     * @param samples
     * A list of samples.
     * @param EndfeetPatches
     * A list of sample triangles.
     */
    AstrocyteMorphology(Samples& samples, EndfeetPatches& endfeetPatches);

    /**
     * @brief AstrocyteMorphology
     * @param h5Samples
     * @param h5Sections
     * @param endfeetPatches
     * @param center
     */
    AstrocyteMorphology(const H5Samples& h5Samples,
                        const H5Sections& h5Sections,
                        const EndfeetPatches& endfeetPatches,
                        const Vector3f& center);

public:

    /**
     * @brief getEndfeetPatches
     * @return Returns a reference to the list of sample triangles.
     */
    EndfeetPatches getEndfeetPatches() const;

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

    /**
     * @brief getSomaCenter
     * Gets the center of the soma.
     *
     * @return
     * The Cartesian center of the soma.
     */
    Vector3f getSomaCenter() const;

    /**
     * @brief getSomaMeanRadius
     * Gets the average radius of the soma.
     * @return
     * The average radius of the soma.
     */
    float getSomaMeanRadius() const;

    /**
     * @brief getSomaMinRadius
     * Gets the minimum radius of the soma.
     * @return
     * The minimum radius of the soma.
     */
    float getSomaMinRadius() const;

    /**
     * @brief getSomaMaxRadius
     * Gets the maximum radius of the soma.
     * @return
     * The maximum radius of the soma.
     */
    float getSomaMaxRadius() const;

    /**
     * @brief getSmallestRadiusInMorphology
     * @return
     */
    float getSmallestRadiusInMorphology() const;

    /**
     * @brief getLargestRadiusInMorphology
     * @return
     */
    float getLargestRadiusInMorphology() const;

    void reIndexMorphology();

private:

    /**
     * @brief _constructMorphologyFromSWC
     * Loads the morphology data from the param swcSamples.
     * @param swcSamples
     * Samples list
     */
    void _constructMorphologyFromSWC(const NeuronSWCSamples& swcSamples);

private:

    /**
     * @brief _samplesTriangles
     * A list of the actual sample triangles of the morphology.
     */
    EndfeetPatches _endfeetPatches;

    /**
     * @brief _somaSamples
     * A list of the actual samples of the morphology soma.
     */
    Samples _somaSamples;

    /**
     * @brief _firstSecion
     * A list of the first neurites secitions of the morphology.
     */
    Sections _firstSections;

    /**
     * @brief _somaCenter
     * The center of the soma.
     */
    Vector3f _somaCenter;

    /**
     * @brief _somaMeanRadius
     * The mean radius of the soma.
     */
    float _somaMeanRadius;

    /**
     * @brief _somaMinRadius
     * The minimum radius of the soma.
     */
    float _somaMinRadius;

    /**
     * @brief _somaMaxRadius
     * The maximum radius of the soma.
     */
    float _somaMaxRadius;

    /**
     * @brief _smallestRadiusInMorphology
     */
    float _smallestRadiusInMorphology;

    /**
     * @brief _largestRadiusInMorphology
     */
    float _largestRadiusInMorphology;
};

/**
 * @brief readAstrocyteMorphology
 * Reads the endfeet morphology from a given file.
 * @param morphologyPath
 * The path to the morphology.
 * @return
 * EndfeetMorphology object.
 */
AstrocyteMorphology* readAstrocyteMorphology(std::string& morphologyPath);

}
