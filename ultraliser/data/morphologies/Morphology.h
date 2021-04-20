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

#ifndef ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H
#define ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H

#include <data/morphologies/Sample.h>
#include <data/morphologies/Section.h>

namespace Ultraliser
{

class Morphology
{
public:

    /**
     * @brief Morphology
     * Constructor
     */
    Morphology();

    /**
     * @brief Morphology
     * Destructor
     */
    virtual ~Morphology();

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

protected:

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

#endif // ULTRALISER_DATA_MORPHOLOGIES_MORPHOLOGY_H
