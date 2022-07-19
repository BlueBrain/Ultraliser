/***************************************************************************************************
 * Copyright (c) 2013
 * Consiglio Nazionale delle Ricerche, Istituto di Matematica Applicata e Tecnologie Informatiche
 * Sezione di Genova, IMATI-GE / CNR
 *
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
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

void AdvancedMesh::importOFF(const std::string &filePath)
{
    // Open the file
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "r")) == nullptr)
    {
        // File not there
        LOG_ERROR("Cannot load mesh file [ %s ]", filePath.c_str());
    }

    // The first line must have OFF keyword
    char streamData[256];
    if (fscanf(filePointer, "%255s", streamData)) { /* NOTHING TO DO */ }
    if (strcmp(streamData, "OFF") || feof(filePointer))
    {
        // File not there
        LOG_ERROR("The mesh file [ %s ] is corrupted!", filePath.c_str());
    }

    // Proceed to the first line
    char *line;
    do
    {
        line = File::readLineFromFile(filePointer);
    } while (line[0] == '#' || line[0] == '\0' || !sscanf(line, "%256s", streamData));

    // Got the first line, that has the data of the mesh
    size_t numberVertices, numberTriangles, numberEdges;
    if (sscanf(line, "%zu %zu %zu", &numberVertices, &numberTriangles, &numberEdges) < 3)
    {
        LOG_ERROR("The mesh file [ %s ] is corrupted!", filePath.c_str());
    }

    // The mesh must have at least three vertices
    if (numberVertices < 3)
    {
        LOG_ERROR("Cannot load meshes with less than 3 vertices!");
    }

    // The mesh must have at least one triangle
    if (numberTriangles < 1)
    {
        LOG_ERROR("Cannot load meshes with no triangles!");
    }

    // Skip and proceed to the actual meat
    File::skipCommentAndBlankLines(filePointer);

    // For all the vertices
    for (size_t i = 0; i < numberVertices; ++i)
    {
        // Vertex
        float x, y, z;
        if (fscanf(filePointer, "%f %f %f", &x, &y, &z) == 3)
        {
            _vertices.appendTail(newVertex(x, y, z));
        }
        else
        {
            LOG_ERROR("Could not read coordinates for vertex [ %" PRIu64 "]", i);
        }
    }

    // Create the extended vertices for the connectivity information
    ExtendedVertex **extendedVertices =
            (ExtendedVertex**) malloc(sizeof(ExtendedVertex*) * numberVertices);
    size_t vertexIndex = 0;
    Node *node;
    AdvancedVertex *vertex;
    FOR_EACH_VERTEX(vertex, node) extendedVertices[vertexIndex++] = new ExtendedVertex(vertex);

    File::skipCommentAndBlankLines(filePointer);

    // Triangles data
    bool triangulate = false;
    for (size_t i = 0; i < numberTriangles; ++i)
    {
        // Scan the triangle data
        size_t i1, i2, i3, i4;
        if (fscanf(filePointer,"%" PRIu64 " %" PRIu64 " %" PRIu64 " %" PRIu64 "",
                   &i4, &i1, &i2, &i3) == 4)
        {
            if (i1 < 0 || i2 < 0 || i3 < 0 || i4 < 3 ||
                i1 > (numberVertices - 1) || i2 > (numberVertices - 1) || i3 > (numberVertices - 1))
            {
                LOG_ERROR("Invalid index at face [ %" PRIu64 "]!", i);
            }

            for (size_t j = 3; j <= i4; j++)
            {
                if (i1 == i2 || i2 == i3 || i3 == i1)
                {
                    // Coincident indexes at triangle! Skipping.
                }
                else
                {
                    if (!createIndexedTriangle(extendedVertices, i1, i2, i3))
                    {
                        // Skipping this triangle.
                    }
                }

                i2 = i3;
                if (j < i4)
                {
                    if (fscanf(filePointer, "%" PRIu64 "", &i3) != 1)
                    {
                        LOG_ERROR("Could not read indexes for face [ %" PRIu64 "]", i);
                    }
                    else
                    {
                        triangulate = true;
                    }
                }
            }
        }
        else
        {
            LOG_ERROR("Could not read indexes for face #[ %" PRIu64 "]", i);
        }
    }

    // Reconstruct the connectivity data
    _closeLoadingSession(filePointer, numberTriangles -1, extendedVertices, (triangulate != 0));
}

void AdvancedMesh::exportOFF(const std::string &filePath)
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Exporting OFF Mesh : [ %s ]", filePath.c_str());

    // Open file for writing
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Cannot write [ %s ]! ", filePath.c_str());
        return;
    }

    // First line OFF
    fprintf(filePointer, "OFF");
    fprintf(filePointer, "\n");

    // Triangle and vertex count
    fprintf(filePointer, "%zu %zu 0", _vertices.numberElements(), _triangles.numberElements());
    fprintf(filePointer, "\n");

    // Vertex data
    Node *node;
    AdvancedVertex *vertex;
    size_t progress = 0;
    size_t numberVertices = _vertices.numberElements();
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

    double *ocds = new double[_vertices.numberElements()];

    int i = 0;
    FOR_EACH_VERTEX(vertex, node) ocds[i++] = vertex->x;

    i = 0;
    FOR_EACH_VERTEX(vertex, node) vertex->x = i++;

    // Triangle data
    size_t numberTriangles = _triangles.numberElements();
    progress = 0;
    TIMER_RESET;
    LOOP_STARTS("Writing Vertices")
    FOR_EACH_NODE(_triangles, node)
    {
        LOOP_PROGRESS(progress, numberTriangles);
        ++progress;

        fprintf(filePointer,"3 %d %d %d",VERTEX_1(node), VERTEX_2(node), VERTEX_3(node));
        fprintf(filePointer, "\n");
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the file
    fclose(filePointer);

    i=0;
    FOR_EACH_VERTEX(vertex, node) vertex->x = ocds[i++];
    delete[] ocds;
}

}
