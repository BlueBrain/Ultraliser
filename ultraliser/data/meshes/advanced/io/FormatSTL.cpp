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
#include <data/meshes/advanced/Defines.hh>
#include <data/meshes/advanced/AdvancedMesh.h>
#include <data/meshes/simple/MeshOperations.h>
#include <data/meshes/advanced/Defines.hh>
#include <utilities/Timer.h>
#include <utilities/Utilities.h>

namespace Ultraliser
{

/**
 * @brief AdvancedMesh::importSTL
 *
 * ASCII STL
 * An ASCII STL file begins with the line
 *      solid name
 * where name is optional string (though if name is omitted there must still be a space after solid).
 * The file continues with any number of triangles, each represented as follows:
 *
 * facet normal ni nj nk
 *      outer loop
 *          vertex v1x v1y v1z
 *          vertex v2x v2y v2z
 *          vertex v3x v3y v3z
 *      endloop
 *  endfacet
 * where each n or v is a floating-point number in sign-mantissa-"e"-sign-exponent format.
 * The file concludes with
 *  endsolid name
 *
 * Binary STL
 * UINT8[80]    – Header                 -     80 bytes
 * UINT32       – Number of triangles    -      4 bytes
 * foreach triangle                      -  50 bytes:
 *      REAL32[3] – Normal vector             - 12 bytes
 *      REAL32[3] – Vertex 1                  - 12 bytes
 *      REAL32[3] – Vertex 2                  - 12 bytes
 *      REAL32[3] – Vertex 3                  - 12 bytes
 *      UINT16    – Attribute byte count      -  2 bytes
 * end
 *
 * @param filePath
 * @return
 */
void AdvancedMesh::importSTL(const std::string &filePath)
{
    // Open the file
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "r")) == nullptr)
    {
        // File not there
        LOG_ERROR("Cannot load mesh file [ %s ]", filePath.c_str());
    }

    // Check if the STL file is BINARY OR ASCII encoded
    bool isBinary = false;
    char keyWord[64] = "";
    if (fscanf(filePointer, "%5s", keyWord)) { /* NOTHING TO DO */ }
    if (strcmp(keyWord, "solid"))
    {
        isBinary = true;
    }

    // In some cases, it might look ASCII just because it begins with 'solid', so double check
    if (!isBinary)
    {
        rewind(filePointer);

        char *line;
        if ((line = File::readLineFromFile(filePointer, 0)) == nullptr)
        {
            isBinary = true;
        }
        else if ((line = File::readLineFromFile(filePointer, 0)) == nullptr)
        {
            isBinary = true;
        }
        else
        {
            sscanf(line, "%64s", keyWord);
            if (strcmp(keyWord, "facet"))
            {
                isBinary = true;
            }
        }
    }

    // If the file is binary
    if (isBinary)
    {
        // Re-open it again with binary encoding
        filePointer = freopen(filePath.c_str(), "rb", filePointer);

        // Reposition the stream position indicator
        fseek(filePointer, 80, SEEK_SET);

        // Read the array (4 elements)
        int numberElements = 0;
        if (fread(&numberElements, 4, 1, filePointer)) { /*DUMMY TO AVOID THE WARNING */ }

        for (int i = 0; i < numberElements; ++i)
        {
            // Unexpected end of file
            char facet[50];
            if (!fread(facet, 50, 1, filePointer))
            {
                LOG_ERROR("Unexpected end of file [ %s ]!", filePath.c_str());
            }

            // This is the normal
            AdvancedPoint normal;
            normal.setValue((*((float*)(facet + 0))),
                            (*((float*)(facet + 4))),
                            (*((float*)(facet + 8))));

            // Vertices
            // AdvancedVertex *v1 = nullptr, *v2 = nullptr, *v3 = nullptr;
            AdvancedVertex *v1 = newVertex((*((float *)(facet + 12))),
                                           (*((float *)(facet + 16))),
                                           (*((float *)(facet + 20))));

            AdvancedVertex *v2 = newVertex((*((float *)(facet + 24))),
                                           (*((float *)(facet + 28))),
                                           (*((float *)(facet + 32))));

            AdvancedVertex *v3 = newVertex((*((float *)(facet + 36))),
                                           (*((float *)(facet + 40))),
                                           (*((float *)(facet + 44))));

            // Update the vertices list
            _vertices.appendHead(v1);
            _vertices.appendHead(v2);
            _vertices.appendHead(v3);

            AdvancedEdge *edge1 = createEdge(v1, v2);
            AdvancedEdge *edge2 = createEdge(v2, v3);
            AdvancedEdge *edge3 = createEdge(v3, v1);

            // Create the triangle
            AdvancedTriangle *triangle;
            if (AdvancedTriangle(edge1, edge2, edge3).getNormal() * normal < 0)
            {
                triangle = createTriangle(edge1, edge3, edge2);
            }
            else
            {
                triangle = createTriangle(edge1, edge2, edge3);
            }
        }
    }
    else
    {
        // Data
        AdvancedPoint normal;
        AdvancedVertex *currentVertex, *vertex1=nullptr, *vertex2=nullptr, *vertex3=nullptr;

        // Proceed with ASCII encoding
        char *line;
        while ((line = File::readLineFromFile(filePointer, 0)) != nullptr)
        {
            // Geth the coordinates
            float x, y, z;
            sscanf(line,"%64s %f %f %f", keyWord, &x, &y, &z);

            // Normal
            if (!strcmp(keyWord, "facet"))
            {
                char normalKeyWord[64] = "";
                sscanf(line,"%64s %64s %f %f %f", keyWord , normalKeyWord, &x, &y, &z);

                // Get the normal

                normal.setValue(x, y, z);
            }
            else
            {
                // Vertex
                if (!strcmp(keyWord, "vertex"))
                {
                    _vertices.appendHead((currentVertex = newVertex(x, y, z)));

                    if (vertex1 == nullptr)
                    {
                        vertex1 = currentVertex;
                    }
                    else if (vertex2 == nullptr)
                    {
                        vertex2 = currentVertex;
                    }
                    else if (vertex3 == nullptr)
                    {
                        vertex3 = currentVertex;

                        AdvancedEdge *edge1 = createEdge(vertex1, vertex2);
                        AdvancedEdge *edge2 = createEdge(vertex2, vertex3);
                        AdvancedEdge *edge3 = createEdge(vertex3, vertex1);

                        AdvancedTriangle *triangle;
                        if (AdvancedTriangle(edge1, edge2, edge3).getNormal() * normal < 0)
                        {
                            triangle = createTriangle(edge1, edge3, edge2);
                        }
                        else
                        {
                            triangle = createTriangle(edge1, edge2, edge3);
                        }

                        // Reset all the vertices
                        vertex1 = vertex2 = vertex3 = nullptr;
                    }
                }
            }
        }
    }

    // Close the file
    fclose(filePointer);

    LOG_INFO("Loaded [%d] vertices and [%d] faces.",
             _vertices.numberElements(), _triangles.numberElements());

    // Rebuild the connectivity
    if (_triangles.numberElements())
        rebuildConnectivity();

}

void AdvancedMesh::exportSTL(const std::string &filePath)
{
    // Start the timer
    TIMER_SET;

    LOG_STATUS("Exporting STL Mesh : [ %s ]", filePath.c_str());

    // If the path does not exist, just set a warning and return
    FILE *filePointer;
    if ((filePointer = fopen(filePath.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Can NOT write the mesh to the following file [ %s ]!", filePath.c_str());
        return;
    }

    // Header
    fprintf(filePointer, "solid Ultraliser");
    fprintf(filePointer, "\n");

    // Triangles and vertices
    Node *node;
    AdvancedTriangle *triangle;
    AdvancedPoint normal;
    uint64_t progress = 0;
    uint64_t numberTriangles = _triangles.numberElements();
    LOOP_STARTS("Writing Vertices, Normals and Triangles")
    FOR_EACH_TRIANGLE(triangle, node)
    {
        LOOP_PROGRESS(progress, numberTriangles);
        ++progress;

        // Get the normal
        normal = triangle->getNormal();

        // Face normal
        fprintf(filePointer, " facet normal %f %f %f", NODE_TO_FLOAT(normal.x),
                                                       NODE_TO_FLOAT(normal.y),
                                                       NODE_TO_FLOAT(normal.z));
        fprintf(filePointer, "\n");

        // Outer loop
        fprintf(filePointer, "  outer loop");
        fprintf(filePointer, "\n");

        // Vertex 1
        fprintf(filePointer, "   vertex %f %f %f", NODE_TO_FLOAT(triangle->v1()->x),
                                                   NODE_TO_FLOAT(triangle->v1()->y),
                                                   NODE_TO_FLOAT(triangle->v1()->z));
        fprintf(filePointer, "\n");

        // Vertex 2
        fprintf(filePointer, "   vertex %f %f %f", NODE_TO_FLOAT(triangle->v2()->x),
                                                   NODE_TO_FLOAT(triangle->v2()->y),
                                                   NODE_TO_FLOAT(triangle->v2()->z));
        fprintf(filePointer, "\n");

        // Vertex 3
        fprintf(filePointer, "   vertex %f %f %f", NODE_TO_FLOAT(triangle->v3()->x),
                                                   NODE_TO_FLOAT(triangle->v3()->y),
                                                   NODE_TO_FLOAT(triangle->v3()->z));
        fprintf(filePointer, "\n");

        // End loop
        fprintf(filePointer, "  endloop");
        fprintf(filePointer, "\n");

        // End face clause
        fprintf(filePointer, " endfacet");
        fprintf(filePointer, "\n");
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Just a terminator
    fprintf(filePointer, "endsolid Ultraliser");
    fprintf(filePointer, "\n");

    // Close the file
    fclose(filePointer);
}

}
