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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H

#include <data/morphologies/Morphology.h>
#include <data/morphologies/h5/VasculatureH5Sample.hh>
#include <data/morphologies/h5/VasculatureH5Section.hh>
#include <data/morphologies/h5/VasculatureH5Connectivity.hh>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/**
 * @brief The VasculatureMorphology class
 * Vasculature morphology
 */
class VasculatureMorphology: public Morphology
{
public:

    /**
     * @brief VasculatureMorphology
     * Constructor
     * @param h5Samples
     * Samples list
     * @param h5Sections
     * Sections list
     * @param h5Connectivity
     * Connectivity list
     */
    VasculatureMorphology(const H5Samples& h5Samples,
                          const VasculatureH5Sections& h5Sections,
                          const VasculatureH5ConnectivityList& h5Connectivity);

    /**
     * @brief VasculatureMorphology
     * Constructor using a list of samples and sections only.
     * @param samples
     * A list of the samples in the morphology
     * @param sections
     * A list of the sections in the morphology.
     */
    VasculatureMorphology(Samples samples, Sections sections);

private:

    /**
     * @brief _constructSections
     * Loads the sections data from the _h5Samples and _h5Sections lists.
     */
    void _constructSections();

    /**
     * @brief _connectSections
     * Construct the connected graph of the morphology using the
     * _h5Connectivity list.
     */
    void _connectSections();

private:

    /**
     * @brief _h5Samples
     * Morphology samples as loaded from the .h5 morphology file.
     */
    const H5Samples _h5Samples;

    /**
     * @brief _h5Sections
     * Morphology sections as loaded from the .h5 morphology file.
     */
    const VasculatureH5Sections _h5Sections;

    /**
     * @brief _h5Connectivity
     * Morphology connectivity list as loaded from the .h5 morphology file.
     */
    const VasculatureH5ConnectivityList _h5Connectivity;
};

/**
 * @brief readVascularMorphology
 * Reads the vascular morphology from a given file.
 * @param morphologyPath
 * The path to the morphology.
 * @return
 * VasculatureMorphology object.
 */
VasculatureMorphology* readVascularMorphology(std::string& morphologyPath);

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H
