/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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
    uint64_t numSamplesPerSeg = samples.size() - 1;
    uint64_t pSamplesPerSeg[numSamplesPerSeg];
    for (uint64_t i = 0; i < numSamplesPerSeg; ++i)
    {
        float dist = (samples[i]->getRadius() + samples[i+1]->getRadius()) *
                ULTRALISER_PIF / numberOfSides;
        uint64_t samplesPerSeg = std::round((samples[i]->getPosition() -
                                    samples[i + 1]->getPosition()).abs() /
                                            dist) + 1;
        if (samplesPerSeg < 2 )
        {
            samplesPerSeg = 2;
        }
        pSamplesPerSeg[i] = samplesPerSeg;
    }
    *this = SectionGeometry(samples, numberOfSides, pSamplesPerSeg);
}

SectionGeometry::SectionGeometry(const Samples& samples, uint64_t numberOfSides,
                                 uint64_t samplesPerSegment)
{
    uint64_t numSamplesPerSeg = samples.size() - 1;
    uint64_t pSamplesPerSeg[numSamplesPerSeg];
    for (uint64_t i = 0; i < numSamplesPerSeg; ++i)
    {
        pSamplesPerSeg[i] = samplesPerSegment;
    }
    *this = SectionGeometry(samples, numberOfSides, pSamplesPerSeg);
}

SectionGeometry::SectionGeometry(const Samples& samples, uint64_t numberOfSides,
                                 uint64_t* pSamplesPerSegment)
{
    uint64_t numSamples = samples.size();
    uint64_t totalSamples = 1;

    // Tangent computation
    Vector3f tangents[numSamples];
    for (uint64_t i = 0; i < numSamples-1; ++i)
    {
        totalSamples += pSamplesPerSegment[i] - 1;
    }

    // Per segment geometry generation
    numVertices = numberOfSides * totalSamples + 2;
    vertices = new Vertex[numVertices];
    CrossSectionGeometry csg(numberOfSides);
    totalSamples = 0;
    for (uint64_t i = 0; i < numSamples - 1; ++i)
    {
        // Per segment tangents compute 
        Vector3f pos0;
        Vector3f pos1 = samples[i]->getPosition();
        Vector3f pos2 = samples[i+1]->getPosition();
        Vector3f pos3;
        if (i == 0)
        {
            pos0 = 2.0f * pos1 - pos2;
        }
        else
        {
            pos0 = samples[i - 1]->getPosition();
        }
        if (i == numSamples - 2)
        {
            pos3 = 2.0f * pos2 - pos1;
        }
        else
        {
            pos3 = samples[i + 2]->getPosition();
        }
        Vector3f p10 = pos1 - pos0;
        Vector3f p20 = pos2 - pos0;
        Vector3f p21 = pos2 - pos1;
        Vector3f p31 = pos3 - pos1;
        Vector3f p32 = pos3 - pos2;
        float t_1 = sqrt((p10).abs());
        float t_2 = t_1 + sqrt(p21.abs());
        float t_3 = t_2 + sqrt(p32.abs());
        float t_21 = t_2 - t_1;
        float t_31 = t_3 - t_1;
        float t_32 = t_3 - t_2;
        Vector3f m1 = (t_2 - t_1) * (p10 / t_1 - p20 / t_2 + p21 / t_21);
        Vector3f m2 = (t_2 - t_1) * (p21 / t_21 - p31 / t_31 + p32 / t_32);
        uint64_t samplesPerSeg = pSamplesPerSegment[i] - 1;
        float r1 = samples[i]->getRadius();
        float r2 = samples[i+1]->getRadius();
        float dt = 1.0f / samplesPerSeg;

        // Spline interpolation and vertices generation 
        if (i == numSamples - 2 )
        {
            ++samplesPerSeg;
        }
        for (uint64_t j = 0; j < samplesPerSeg; ++j)
        {
            float t = dt * j;
            float t2 = t * t;
            float t3 = t2 * t;
            Vector3f pos = (2 * t3 - 3 * t2 + 1) * pos1 +
                    (t3 - 2 * t2 + t) * m1 +
                    (-2 * t3 + 3 * t2) * pos2 +
                    (t3 - t2) * m2;
            Vector3f tangent =  (6 * t2 - 6 * t) * pos1 +
                    (3 * t2 - 4 * t + 1) * m1 +
                    (-6 * t2 + 6 * t) * pos2 +
                    (3 * t2 - 2* t) * m2;; 
            float radius = (2 * t3 - 3 * t2 + 1) * r1 +
                    (-2 * t3 + 3 * t2) * r2;

            csg.setPositionOrientationRadius(pos, tangent, radius);
            std::memcpy(vertices[(totalSamples + j) * numberOfSides],
                        csg.vertices, numberOfSides*sizeof(Vertex));
        }
        totalSamples += samplesPerSeg;
    }
    vertices[numVertices-2] = samples[0]->getPosition();
    vertices[numVertices-1] = samples[numSamples-1]->getPosition();

    // Primitives assembly 
    numTriangles = numberOfSides * (totalSamples) * 2;
    triangles = new Triangle[numTriangles];
    for (uint64_t i = 0; i < totalSamples - 1; ++i)
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
        triangles[numberOfSides*(totalSamples-1)*2+i] =
                Vec3i_64(index0, index1, index2);
    }
    index0 = numVertices - 1;
    for (unsigned int i = 0; i < numberOfSides; ++i)
    {
        uint64_t index1 = ((totalSamples-1)*numberOfSides) + i;
        uint64_t index2 = ((totalSamples-1)*numberOfSides) + (i+1)%numberOfSides;
        triangles[numberOfSides*((totalSamples-1)*2+1)+i] =
                Vec3i_64(index0, index1, index2);
    }
}

CrossSectionGeometry::CrossSectionGeometry(uint64_t numVertices)
    : numVertices(numVertices)
    , _position(Vector3f(0.0f, 0.0f, 0.0f))
    , _orientation(Vector3f(0.0f, 0.0f, 1.0f))
    , _radius(1.0)
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
