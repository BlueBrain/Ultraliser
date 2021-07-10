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

#include <common/Headers.hh>
#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/meshes/simple/MeshOperations.h>
#include <data/meshes/advanced/Defines.hh>
#include <utilities/Timer.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

void AdvancedMesh::importPLY(const std::string &fileName)
{
    Vertices vertices;
    Triangles triangles;
    std::string filePath = std::string(fileName);

    // Import the file
    Ultraliser::importPLY(filePath, vertices, triangles);

    Utilities::Timer statsTimer;
    statsTimer.start();

    // Generic
    AdvancedVertex* vertex;
    Node* node;

    // Fill the _vertices list
    for(uint64_t i = 0; i < vertices.size(); ++i)
        _vertices.appendTail(newVertex(F2D(vertices[i].x()),
                                       F2D(vertices[i].y()),
                                       F2D(vertices[i].z())));

    ExtendedVertex** vertexList = nullptr;
    vertexList =  static_cast<ExtendedVertex**>
            (malloc(sizeof(ExtendedVertex*) * vertices.size()));

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

        Triangle triangle = triangles[i];
        if (createIndexedTriangle(vertexList, triangle[0], triangle[1], triangle[2])) { }
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Remove the extended vertices
    if (vertexList != nullptr)
    {
        TIMER_RESET;
        LOOP_STARTS("Cleaning Extended Data");

        for (uint64_t i = 0; i < vertices.size(); ++i)
            delete(vertexList[i]);
        free(vertexList);

        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);

        // Fix the connectivity of the mesh
        fixConnectivity();

        // Update the data
        eulerUpdate();
        _dBoundaries = _dHandles = _dShells = 1;
    }
}

void AdvancedMesh::exportPLY(const std::string &filePath, bool writeASCII)
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Exporting PLY Mesh : [ %s ]", filePath.c_str());

    // Open file for writing
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Cannot write [ %s ]! ", filePath.c_str());
        return;
    }

    // Header data
    fprintf(filePointer, "ply");
    fprintf(filePointer, "\n");

    if (writeASCII)
    {
        fprintf(filePointer, "format ascii 1.0");
        fprintf(filePointer, "\n");
    }
    else
    {
        fprintf(filePointer, "format binary_little_endian 1.0");
        fprintf(filePointer, "\n");
    }

    // Vertex count
    fprintf(filePointer, "element vertex %" PRIu64 "",_vertices.numberElements());
    fprintf(filePointer, "\n");

    fprintf(filePointer, "property float x");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property float y");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property float z");
    fprintf(filePointer, "\n");

    // Triangle count
    fprintf(filePointer, "element face %" PRIu64 "",_triangles.numberElements());
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property list uchar int vertex_indices");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "end_header");
    fprintf(filePointer, "\n");

    Node *node;
    AdvancedVertex *vertex;
    uint64_t numberVertices = _vertices.numberElements();
    uint64_t progress = 0;
    if (writeASCII)
    {
        LOOP_STARTS("Writing Vertices")
        FOR_EACH_VERTEX(vertex, node)
        {
            LOOP_PROGRESS(progress, numberVertices);
            ++progress;

            fprintf(filePointer, "%f %f %f",
                    NODE_TO_FLOAT(vertex->x), NODE_TO_FLOAT(vertex->y), NODE_TO_FLOAT(vertex->z));
             fprintf(filePointer, "\n");
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOOP_STARTS("Writing Vertices")
        FOR_EACH_VERTEX(vertex, node)
        {
            LOOP_PROGRESS(progress, numberVertices);
            ++progress;

            float floatCoordinates[3];
            floatCoordinates[0] = NODE_TO_FLOAT(vertex->x);
            floatCoordinates[1] = NODE_TO_FLOAT(vertex->y);
            floatCoordinates[2] = NODE_TO_FLOAT(vertex->z);

            fwrite(floatCoordinates, sizeof(float), 3, filePointer);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Auxiliary array to construct the triangles data
    double *auxiliaryVertices = new double[_vertices.numberElements()];

    uint64_t i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        auxiliaryVertices[i++] = vertex->x;
    }

    i = 0;
    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = i++;
    }

    // Writing the triangles
    uint64_t numberTriangles = _triangles.numberElements();
    if (writeASCII)
    {
        progress = 0;
        TIMER_RESET;
        LOOP_STARTS("Writing Triangles")
        FOR_EACH_NODE(_triangles, node)
        {
            LOOP_PROGRESS(progress, numberTriangles);
            ++progress;

            fprintf(filePointer, "3 %d %d %d", VERTEX_1(node), VERTEX_2(node), VERTEX_3(node));
            fprintf(filePointer, "\n");
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        progress = 0;
        TIMER_RESET;
        LOOP_STARTS("Writing Triangles")
        FOR_EACH_NODE(_triangles, node)
        {
            LOOP_PROGRESS(progress, numberTriangles);
            ++progress;

            int triangleVertexIndex[3];
            triangleVertexIndex[0] = VERTEX_1(node);
            triangleVertexIndex[1] = VERTEX_2(node);
            triangleVertexIndex[2] = VERTEX_3(node);

            const int elementCount = 3;
            fwrite(&elementCount, sizeof(unsigned char), 1, filePointer);
            fwrite(triangleVertexIndex, sizeof(int), 3, filePointer);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Close the file
    fclose(filePointer);

    i = 0;

    FOR_EACH_VERTEX(vertex, node)
    {
        vertex->x = auxiliaryVertices[i++];
    }

    // Clear the auxiliary array
    delete[] auxiliaryVertices;
}

}
