/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Juan Jose Garcia Cantero < juanjose.garcia@epfl.ch>
 *      Marwan Abdellah < marwan.abdellah@epfl.ch >
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

#include "SectionGeometry.h"

namespace Ultraliser
{

SectionGeometry::SectionGeometry(const Samples& samples, uint64_t numberOfSides)
{
    uint64_t numSamples = samples.size();

    if (numSamples > 1)
    {

    // Tangent computation
    Vector3f tangents[numSamples];
    for (uint64_t i = 0; i < numSamples; ++i)
    {
        Vector3f tangent;
        if (i == 0)
            tangent = samples[1]->getPosition() - samples[0]->getPosition();
        else if (i == numSamples - 1)
            tangent = samples[i]->getPosition() - samples[i - 1]->getPosition();
        else 
            tangent = samples[i + 1]->getPosition() - samples[i - 1]->getPosition();
        tangents[i] = tangent.normalized();
    }

    // Per sample geometry generation
    numVertices = numberOfSides * numSamples + 2;
    vertices = new Vertex[numVertices];
    CrossSectionGeometry csg(numberOfSides);
    const float rollOffset = M_PI / numberOfSides;
    for (uint64_t i = 0; i < numSamples; ++i)
    {
        auto sample = samples[i];
        csg.setPositionOrientationRadius(sample->getPosition(), tangents[i], sample->getRadius());
        csg.setRoll(rollOffset*i);
        std::memcpy(vertices[i * numberOfSides], csg.vertices, numberOfSides*sizeof(Vertex));
    }
    vertices[numVertices-2] = samples[0]->getPosition();
    vertices[numVertices-1] = samples[numSamples-1]->getPosition();

    // Primitives assembly 
    numTriangles = numberOfSides * numSamples * 2;
    triangles = new Triangle[numTriangles];
    for (uint64_t i = 0; i < numSamples - 1; ++i)
    {
        for (uint64_t j = 0; j < numberOfSides; ++j)
        {
            uint64_t index0 = i * numberOfSides + j;
            uint64_t index1 = i * numberOfSides + ((j + 1) % numberOfSides);
            uint64_t index2 = (i + 1) * numberOfSides + j;
            uint64_t index3 = (i + 1) * numberOfSides + ((j + 1) % numberOfSides);
            triangles[(i * numberOfSides +j) * 2] =
                    Vec3i_64(index0, index1, index2);
            triangles[(i * numberOfSides + j)*2 + 1] =
                    Vec3i_64(index1, index3, index2);
        }
    }

    // Cap first and last samples
    uint64_t index0 = numVertices - 2;
    for (unsigned int i = 0; i < numberOfSides; ++i)
    {
        uint64_t index2 = i;
        uint64_t index1 = (i+1)%numberOfSides;
        triangles[numberOfSides*(numSamples-1)*2+i] =
                Vec3i_64(index0, index1, index2);
    }
    index0 = numVertices - 1;
    for (unsigned int i = 0; i < numberOfSides; ++i)
    {
        uint64_t index1 = ((numSamples-1)*numberOfSides) + i;
        uint64_t index2 = ((numSamples-1)*numberOfSides) + (i+1)%numberOfSides;
        triangles[numberOfSides*((numSamples-1)*2+1)+i] =
                Vec3i_64(index0, index1, index2);
    }}
}

CrossSectionGeometry::CrossSectionGeometry(uint64_t numVertices)
    : numVertices(numVertices)
    , _position(Vector3f(0.0f, 0.0f, 0.0f))
    , _orientation(Vector3f(0.0f, 0.0f, 1.0f))
    , _radius(1.0f)
    , _roll(0.0f)
{
    // Construct the vertices list
    vertices = new Vertex[numVertices];
    const float angleInc = 2 * ULTRALISER_PIF / numVertices;
    for (uint64_t i = 0; i < numVertices; ++i)
    {
        const float angle = angleInc * i;
        vertices[i] = Vertex(cos(angle), sin(angle), 0.0f);
    }
}

CrossSectionGeometry::CrossSectionGeometry(uint64_t numVertices,
                                           const Vector3f& position,
                                           const Vector3f& orientation,
                                           float radius)
    : CrossSectionGeometry(numVertices)
{
    setPositionOrientationRadius(position, orientation, radius);
}

CrossSectionGeometry::CrossSectionGeometry(const CrossSectionGeometry& input)
    : numVertices(input.numVertices)
    , _position(input._position)
    , _orientation(input._orientation)
    , _radius(input._radius)
{
    vertices = new Vertex[numVertices];
    std::memcpy(vertices, input.vertices, numVertices*sizeof(Vertex));
}

CrossSectionGeometry::~CrossSectionGeometry()
{
    delete vertices;
}

void CrossSectionGeometry::setPosition(const Vector3f& position)
{
    const Vector3f posInc = position - _position;
    for (uint64_t i = 0; i < numVertices; i++)
    {
        vertices[i] += posInc;
    }
    _position = position;
}

void CrossSectionGeometry::setOrientation(const Vector3f& orientation)
{
    const Quat4f q = Quat4f::fromTwoVectors(_orientation, orientation);
    for (uint64_t i = 0; i < numVertices; ++i)
    {
        vertices[i] = q.rotate(vertices[i] - _position) + _position;
    }
    _orientation = orientation;
}

void CrossSectionGeometry::setRadius(const float& radius)
{
    const float radiusRatio = radius / _radius;
    for (uint64_t i = 0; i < numVertices; ++i)
    {
        vertices[i] = (vertices[i] - _position) * radiusRatio + _position;
    }
    _radius = radius;
}

void CrossSectionGeometry::setRoll(const float& roll)
{
    const float rollDiff = roll - _roll;
    Quat4f q;
    q.setAxisAngle(rollDiff, _orientation);
    
    for (uint64_t i = 0; i < numVertices; ++i)
    {
        vertices[i] = q.rotate(vertices[i] - _position) + _position;
    }

    _roll = roll;
}

void CrossSectionGeometry::setPositionOrientationRadius(const Vector3f& position,
                                                        const Vector3f& orientation,
                                                        float radius)
{
    const Quat4f q = Quat4f::fromTwoVectors(_orientation, orientation);
    const float radiusRatio = radius / _radius;

    for (uint64_t i = 0; i < numVertices; ++i)
    {
        vertices[i] = q.rotate(vertices[i] - _position) * radiusRatio + position;
    }

    _position = position;
    _orientation = orientation;
    _radius = radius;
}

SectionGeometry::~SectionGeometry()
{
    /// EMPTY
}


}
