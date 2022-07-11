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
    std::vector<Vertex> newVertices;
    uint64_t numNewSamples = 0;
    CrossSectionGeometry csg(numberOfSides);
    const float rollOffset = M_PI / numberOfSides;
    for (uint64_t i = 0; i < numSamples - 1; ++i) 
    {
        // Linear interpolation per section
        Vector3f pos0 = samples[i]->getPosition();
        Vector3f pos1 = samples[i + 1]->getPosition();
        float radius0 = samples[i]->getRadius();
        float radius1 = samples[i + 1]->getRadius();
        Quat4f q0(tangents[i]);
        Quat4f q1(tangents[i + 1]);

        float meanRadius = (radius0 + radius1) * 0.5f;
        uint32_t numDivision = std::ceil((pos1 - pos0).abs() / (2.0f * M_PI * meanRadius /         
                                                                numberOfSides));
        float alphaIncr = 1.0f / numDivision;

        for (uint32_t j = 0; j < numDivision; ++j) 
        {
            float alpha = alphaIncr * j;
            Vector3f pos = Vector3f::lerp(pos0, pos1, alpha);
            float radius = radius0 * (1.0f - alpha) + radius1 * alpha;
            Vector3f tangent = (Quat4f::slerp(q0, q1, alpha)).xyz();
            csg.setPositionOrientationRadius(pos, tangent, radius);
            csg.setRoll(rollOffset * numNewSamples);
            ++numNewSamples;
            for (uint32_t vertexId = 0; vertexId < numberOfSides; ++vertexId)
            {
                newVertices.push_back(csg.vertices[vertexId]);
            }
        }
    }
    
    // Add the last section sample geometry
    auto sample = samples[numSamples - 1];
    csg.setPositionOrientationRadius(sample->getPosition(), tangents[numSamples - 1],           
                                     sample->getRadius());
    csg.setRoll(rollOffset * numNewSamples);
    ++numNewSamples;
    for (uint32_t vertexId = 0; vertexId < numberOfSides; ++vertexId)
      newVertices.push_back(csg.vertices[vertexId]);

    // Add last two vertices to cap begin and end of the section
    newVertices.push_back(samples[0]->getPosition());
    newVertices.push_back(samples[numSamples - 1]->getPosition());

    // Copy generated vertices to the final structure
    numVertices = newVertices.size();
    vertices = new Vertex[numVertices];
    std::copy(newVertices.begin(), newVertices.end(), vertices);
    newVertices.clear();

    // Primitives assembly 
    numTriangles = numberOfSides * numNewSamples * 2;
    triangles = new Triangle[numTriangles];
    for (uint64_t i = 0; i < numNewSamples - 1; ++i) 
    {
        for (uint64_t j = 0; j < numberOfSides; ++j) 
        {
            uint64_t index0 = i * numberOfSides + j;
            uint64_t index1 = i * numberOfSides + ((j + 1) % numberOfSides);
            uint64_t index2 = (i + 1) * numberOfSides + j;
            uint64_t index3 = (i + 1) * numberOfSides + ((j + 1) % numberOfSides);
            triangles[(i * numberOfSides + j) * 2] = Vec3i_64(index0, index1, index2);
            triangles[(i * numberOfSides + j) * 2 + 1] = Vec3i_64(index1, index3, index2);
        }
    }

    // Cap first and last samples
    uint64_t index0 = numVertices - 2;
    for (unsigned int i = 0; i < numberOfSides; ++i)
    {
        uint64_t index2 = i;
        uint64_t index1 = (i + 1) % numberOfSides;
        triangles[numberOfSides * (numNewSamples - 1) * 2 + i] = Vec3i_64(index0, index1, index2);
    }
    index0 = numVertices - 1;
    for (unsigned int i = 0; i < numberOfSides; ++i)
    {
        uint64_t index1 = ((numNewSamples - 1)* numberOfSides) + i;
        uint64_t index2 = ((numNewSamples - 1)* numberOfSides) + (i + 1) % numberOfSides;
        triangles[numberOfSides*((numNewSamples - 1) * 2 + 1) + i] = 
            Vec3i_64(index0, index1, index2);
    }
    }
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
