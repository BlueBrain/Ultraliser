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

#ifndef ULTRALISER_SECTION_GEOMETRY_H
#define ULTRALISER_SECTION_GEOMETRY_H

#include <data/meshes/simple/primitives/Primitives.h>
#include <data/morphologies/Section.h>

namespace Ultraliser
{

class CrossSectionGeometry;

/**
 * @brief The SectionGeometry class
 */
class SectionGeometry
{
public:

    /**
     * @brief CrossSectionGeometry
     * Constructor
     *
     * @param samples
     * Polyline samples to generate mesh from.
     * @param numberOfSides
     * Per sample geometry number of sides
     */
    SectionGeometry(const Samples& samples, uint64_t numberOfSides);

    /**
     * @brief CrossSectionGeometry
     * Constructor.
     * @param samples
     * polyline samples to generate mesh from.
     * @param numberOfSides
     * per sample geometry number of sides
     * @param samplesPerSegment
     * number of samples per polyline segment
     */
    SectionGeometry(const Samples& samples, uint64_t numberOfSides,
                    uint64_t samplesPerSegment);

    /**
     * @brief CrossSectionGeometry
     * Constructor.
     * @param samples
     * polyline samples to generate mesh from.
     * @param numberOfSides
     * per sample geometry number of sides
     * @param pSamplesPerSegment
     * numbers vector of samples per polyline segment
     */
    SectionGeometry(const Samples& samples, uint64_t numberOfSides,
                    uint64_t* pSamplesPerSegment);
    
    virtual ~SectionGeometry();

public: 

    /**
     * @brief vertices
     * A list of all the vertices in the section geometry.
     */
    Vertex* vertices;

    /**
     * @brief numVertices
     * Number of vertices of the section geometry.
     */
    uint64_t numVertices;

    /**
     * @brief triangles
     * A list of all the triangles in the section geometry.
     */
    Triangle* triangles;
    
    /**
     * @brief numTriangles
     * Number of triangles of the section geometry.
     */
    uint64_t numTriangles;
    
};

/**
 * @brief The CrossSectionGeometry class
 */
class CrossSectionGeometry
{
public:

    /**
     * @brief CrossSectionGeometry
     * Constructor.
     * @param numVertices
     * Number of the vertices of the cross section geometry.
     */
    CrossSectionGeometry(uint64_t numVertices);
    
    /**
     * @brief CrossSectionGeometry
     * Constructor.
     * @param numVertices
     * Number of the vertices of the cross section geometry.
     * @Param position
     * The position of the cross section geometry.
     * @param orientation
     * The orientation of the cross section geometry.
     * @param radius
     * The radius of the cross section geometry.
     */
    CrossSectionGeometry(uint64_t numVertices,
                         const Vector3f& position,
                         const Vector3f& orientation,
                         float radius);

    /**
     * @brief CrossSectionGeometry
     * Copy constructor.
     * @param input
     * Input cross section geometry to be copied.
     */
    CrossSectionGeometry(const CrossSectionGeometry& input);

    /**
     * @brief setPosition
     * Updates the position of the cross section geometry.
     * @param position
     * The new position.
     */
    void setPosition(const Vector3f& position);
    
    /**
     * @brief setOrientation
     * Updates the oriantation of the cross section geometry.
     * @param orientation
     * The new orientation.
     */
    void setOrientation(const Vector3f& orientation);

    /**
     * @brief setRadius
     * Updates the radius of the cross section geometry.
     * @param radius
     * The new radius.
     */
    void setRadius(const float& radius);

    /**
     * @brief setPosition
     * Updates the position, orientation and radius of the cross section
     * geometry.
     * @param position
     * The new position.
     * @param orientation
     * The new orientation.
     * @param radius
     * The new radius.
     */
    void setPositionOrientationRadius(const Vector3f& position,
                                      const Vector3f& orientation,
                                      float radius);
    
public:

    /**
     * @brief vertices
     * A list of vertices of the cross section geometry.
     */
    VertexPtr vertices;

    /**
     * @brief numVertices
     * Number of vertices of the cross section geometry.
     */
    uint64_t numVertices;

private:

    /**
     * @brief _position
     * The position of the cross section geometry.
     */
    Vector3f _position;
    
    /**
     * @brief _orientation
     * The orientation of the cross section geometry.
     */
    Vector3f _orientation;

    /**
     * @brief _radius
     * The radius of the cross section geometry.
     */
    float _radius;

public:
    ~CrossSectionGeometry(void);
};

}

#endif // ULTRALISER_SECTION_GEOMETRY_H
