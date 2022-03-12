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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_ASTROCYTE_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_ASTROCYTE_MORPHOLOGY_H

#include <data/morphologies/Morphology.h>
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
     */
    AstrocyteMorphology(const H5Samples& h5Samples,
                        const H5Sections& h5Sections,
                        EndfeetPatches& endfeetPatches);

public:


    /**
     * @brief getEndfeetPatches
     * @return Returns a reference to the list of sample triangles.
     */
    EndfeetPatches getEndfeetPatches() const;

private:

    void _constructSkeleton();

private:

    /**
     * @brief _samplesTriangles
     * A list of the actual sample triangles of the morphology.
     */
    EndfeetPatches _endfeetPatches;
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

#endif // ULTRALISER_DATA_MORPHOLOGIES_ASTROCYTE_MORPHOLOGY_H
