/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Marco Attene < IMATI-GE / CNR >
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
 *
 * The content of this file is based on MeshFix. The code has been modified under the terms of
 * the GNU General Public License as published by the Free Software Foundation either version 3 of
 * the License, or (at your option) any later version.
 * MeshFix has a dual license for free and commercial use. For further information, please refer
 * to the original repository at < https://github.com/MarcoAttene/MeshFix-V2.1>.
 **************************************************************************************************/

#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/meshes/simple/MeshOperations.h>
#include <data/meshes/advanced/Defines.hh>
#include <utilities/Timer.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

void AdvancedMesh::importOBJ(const std::string &fileName)
{
    Vertices vertices;
    Triangles triangles;

    // Import the file
    Ultraliser::importOBJ(fileName, vertices, triangles);

    // Generic
    AdvancedVertex* vertex;
    Node* node;

    // Fill the _vertices list
    for(uint64_t i = 0; i < vertices.size(); ++i)
    {
        _vertices.appendTail(newVertex(F2D(vertices[i].x()),
                                       F2D(vertices[i].y()),
                                       F2D(vertices[i].z())));
    }

    ExtendedVertex** vertexList = nullptr;
    vertexList =  static_cast< ExtendedVertex** >(malloc(sizeof(ExtendedVertex*) * vertices.size()));

    TIMER_SET;
    LOOP_COUNTER_SET;
    LOOP_STARTS("Creating Vertex List");
    FOR_EACH_VERTEX(vertex, node)
    {
        LOOP_PROGRESS_FRACTION(COUNTER, vertices.size());

        // Add the vertex to the vertex list
        vertexList[COUNTER++] = new ExtendedVertex(vertex);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Set the faces
    TIMER_RESET;
    LOOP_COUNTER_RESET;
    LOOP_STARTS("Creating Face List");
    for(uint64_t i = 0; i < triangles.size(); ++i)
    {
        LOOP_PROGRESS_FRACTION(i, triangles.size());

        // Get a reference to the triangle
        Triangle triangle = triangles[i];

        // Create an indexed triangle with edges stored
        createIndexedTriangle(vertexList, triangle[0], triangle[1], triangle[2]);
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Remove the extended vertices
    if (vertexList != nullptr)
    {
        for (uint64_t i = 0; i < vertices.size(); ++i)
        {
            delete(vertexList[i]);
        }
        free(vertexList);

        // Fix the connectivity of the mesh
        fixConnectivity();

        // Update the data
        eulerUpdate();
        _dBoundaries = _dHandles = _dShells = 1;
    }
}

void AdvancedMesh::exportOBJ(const std::string &filePath)
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Exporting OBJ Mesh : [ %s ]", filePath.c_str());

    // Open file for writing
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Cannot write [ %s ]! ", filePath.c_str());
        return;
    }

    // Generic
    Node *node;
    AdvancedVertex *vertex;

    // Write vertices
    uint64_t progress = 0;
    uint64_t numberVertices = _vertices.numberElements();
    LOOP_STARTS("Writing Vertices")
    FOR_EACH_VERTEX(vertex, node)
    {
        LOOP_PROGRESS(progress, numberVertices);
        ++progress;
        fprintf(filePointer, "v %f %f %f",
                NODE_TO_DOUBLE(vertex->x),
                NODE_TO_DOUBLE(vertex->y),
                NODE_TO_DOUBLE(vertex->z));
        fprintf(filePointer, "\n");
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Copy the vertices
    float *auxVertices = new float[I2UI64(_vertices.numberElements())];

    // Construct the faces from the vertices
    uint64_t i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        auxVertices[i++] = D2F(vertex->x);
    }

    // Reset the counter
    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = i++;
    }

    // Triangle data
    uint64_t numberTriangles = _triangles.numberElements();
    progress = 0;
    TIMER_RESET;
    LOOP_STARTS("Writing Triangles")
    FOR_EACH_NODE(_triangles, node)
    {
        LOOP_PROGRESS(progress, numberTriangles);
        ++progress;

        fprintf(filePointer, "f %d %d %d",
                VERTEX_1(node) + 1, VERTEX_2(node) + 1, VERTEX_3(node) + 1);
        fprintf(filePointer, "\n");
    }

    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    fclose(filePointer);

    // Reset the counter
    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = F2D(auxVertices[i++]);
    }

    // Free the auxiliary array
    delete[] auxVertices;
}

}
