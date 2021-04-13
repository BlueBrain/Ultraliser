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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H

#include <data/morphologies/Sample.h>
#include <data/morphologies/Section.h>
#include <data/morphologies/VasculatureH5Sample.hh>
#include <data/morphologies/VasculatureH5Section.hh>
#include <data/morphologies/VasculatureH5Connectivity.hh>

namespace Ultraliser
{

class VasculatureMorphology
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
    VasculatureMorphology(const VasculatureH5Samples& h5Samples,
                          const VasculatureH5Sections& h5Sections,
                          const VasculatureH5ConnectivityList& h5Connectivity);

public:

    /**
     * @brief getSections
     * @return Returns a reference to the list of sections in the morphology.
     */
    Sections getSections() const;

    /**
     * @brief getConnectedPathsFromParentsToChildren
     * @param section
     * @return
     */
    Paths getConnectedPathsFromParentsToChildren(const Section* section) const;

    /**
     * @brief getBoundingBox
     * Returns the bounding box of the morphology.
     * @param pMin
     * The computed pMin of the bounding box.
     * @param pMax
     * The computed pMax of the bounding box.
     * @param bounds
     * The computed bounds of the bounding box.
     * @param center
     * The center of the bounding box.
     */
    void getBoundingBox(Vector3f& pMin, Vector3f& pMax,
                        Vector3f& bounds, Vector3f &center);

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

    /**
     * @brief getConnectedSamplesFromParentToChild
     * @return Returns a list of samples corresponding to each section in
     * the morphology starting from the parents to the children.
     */
    std::vector< Samples > _getConnectedSamplesFromParentToChild() const;

private:

    /**
     * @brief _h5Samples
     * Morphology samples as loaded from the .h5 morphology file.
     */
    const VasculatureH5Samples _h5Samples;

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

    /**
     * @brief _samples
     * A list of all the actual samples of the morphology.
     */
    Samples _samples;

    /**
     * @brief _sections
     * A list of all the actual sections of the morphology.
     */
    Sections _sections;

    /**
     * @brief _pMin
     */
    Vector3f _pMin;

    /**
     * @brief _pMax
     */
    Vector3f _pMax;
};

}

#endif // ULTRALISER_DATA_MORPHOLOGIES_VASCULATURE_MORPHOLOGY_H
