/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechnique Federale de Lausanne (EPFL)
 *
 * Author(s)
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

#include "MeshOperations.h"
#include "TriangleOperations.h"
#include <math/Math.h>
#include <common/Common.h>
#include <utilities/Utilities.h>
#include <utilities/TypeConversion.h>

namespace Ultraliser
{

void computeMeshBoundingBox(const Vertex* vertices,
                            const int64_t &numberVertices,
                            Vector3f& pMin, Vector3f& pMax,
                            const bool& verbose)
{
    // Starting the timer
    TIMER_SET;

    pMin.x() = std::numeric_limits<float>::max();
    pMin.y() = std::numeric_limits<float>::max();
    pMin.z() = std::numeric_limits<float>::max();

    pMax.x() = std::numeric_limits<float>::lowest();
    pMax.y() = std::numeric_limits<float>::lowest();
    pMax.z() = std::numeric_limits<float>::lowest();

    if (verbose) LOOP_STARTS("Computing Bounding Box");
    for (int64_t i = 0; i < numberVertices; ++i)
    {
        Vertex v = vertices[i];

        if (verbose) LOOP_PROGRESS_FRACTION(i, numberVertices);

        if (v.x() > pMax.x()) pMax.x() = v.x();
        if (v.y() > pMax.y()) pMax.y() = v.y();
        if (v.z() > pMax.z()) pMax.z() = v.z();

        if (v.x() < pMin.x()) pMin.x() = v.x();
        if (v.y() < pMin.y()) pMin.y() = v.y();
        if (v.z() < pMin.z()) pMin.z() = v.z();
    }
    if (verbose) LOOP_DONE;

    // Statistics
    if (verbose) LOG_STATS(GET_TIME_SECONDS);
}

float computeMeshVolume(const Vertex *vertices,
                        const Triangle* triangles,
                        const int64_t &numberTriangles,
                        const bool& verbose)
{
    // Starting the timer
    TIMER_SET;

    // Make an array to compute the area per triangle and then sum them up
    float* triangelsSignedVolumes = new float[I2UI64(numberTriangles)];

    if (verbose) LOOP_STARTS("Computing Mesh Volume");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0; i < numberTriangles; ++i)
    {
        // Update the progress bar
        if (verbose) LOOP_PROGRESS(PROGRESS, numberTriangles);
        PROGRESS_UPDATE;

        Vertex p1 = vertices[triangles[i][0]];
        Vertex p2 = vertices[triangles[i][1]];
        Vertex p3 = vertices[triangles[i][2]];

        triangelsSignedVolumes[i] = computeTriangleSignedVolume(p1, p2, p3);
    }
    if (verbose) LOOP_DONE;

    // Compute the total volume
    float meshVolume = 0.0f;
    for (int64_t i = 0; i < numberTriangles; ++i)
    {
        meshVolume += triangelsSignedVolumes[i];
    }

    // Statistics
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    return meshVolume;
}

float computeMeshSurfaceArea(const Vertex* vertices,
                             const Triangle* triangles,
                             const int64_t &numberTriangles,
                             const bool &verbose)
{
    // Starting the timer
    TIMER_SET;

    // Make an array to compute the area per triangle and then sum them up
    float* trianglesArea = new float[I2UI64(numberTriangles)];

    if (verbose) LOOP_STARTS("Computing Mesh Surface Area");
    PROGRESS_SET;
    OMP_PARALLEL_FOR
    for (int64_t i = 0; i < numberTriangles; ++i)
    {
        // Update the progress bar
        if (verbose) LOOP_PROGRESS(PROGRESS, numberTriangles);
        PROGRESS_UPDATE;

        // Get the triangle
        Vec3i_64 triangleIndices = triangles[i];

        // Get the points
        Vector3f p0 = vertices[triangleIndices[0]];
        Vector3f p1 = vertices[triangleIndices[1]];
        Vector3f p2 = vertices[triangleIndices[2]];

        trianglesArea[i] = computeTriangleSurfaceArea(p0, p1, p2);
    }
    if (verbose) LOOP_DONE;

    // Compute the total surface area
    float meshSurfaceArea = 0.0f;
    for (int64_t i = 0; i < numberTriangles; ++i)
    {
        meshSurfaceArea += trianglesArea[i];
    }

    // Statistics
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    return meshSurfaceArea;
}

void importOBJ(const std::string &filePath, Vertices& vertices, Triangles& triangles,
               const bool &verbose)
{
    // Start the timer
    TIMER_SET;

    // Get the number of files in the line
    if (verbose) LOOP_STARTS("Handling File");
    auto numberLines = File::getNumberLinesInFile(filePath);
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    std::ifstream file;
    file.open(filePath.c_str());

    // Parsing line by line
    std::string line;

    // Tokenize each line
    std::string token;

    // Parse the file
    TIMER_RESET;
    if (verbose) LOOP_STARTS("Loading Data")
    size_t progress = 0;
    for (size_t lineNumber = 0; lineNumber < numberLines; ++lineNumber)
    {
        if (verbose) LOOP_PROGRESS_FRACTION(progress, numberLines);

        std::getline(file, line);

        // End of file
        if (file.eof())
            break;

        // End of data
        if (line == OBJ_END_FLAG)
            break;

        // Wrong line
        if (line.size() < 3)
            continue;

        // Comment
        if (line.at(0) == OBJ_HASH)
            continue;

        // Get the line
        std::stringstream dataStream(line);
        dataStream >> token;

        // Vertex normal
        if (token == OBJ_VERTEX_NORMAL_FLAG)
            continue;

        // Vertex
        if (token == OBJ_VERTEX_FLAG)
        {
            Vertex vertex;
            dataStream >> vertex[0] >> vertex[1] >> vertex[2];
            vertices.push_back(vertex);
        }

        // Face
        else if (token == OBJ_FACE_FLAG)
        {
            // An item in the line
            std::string item;

            // Copy the line to a stream where we can process
            std::stringstream faceLineStream(line);

            // A list of face vertices (and only the vertices)
            std::vector< int64_t > faceVertices;

            // Consume the face 'f' token to proceed to the actual faces
            faceLineStream >> token;
            while (faceLineStream >> item)
            {
                // Replace all the '/' from the item with spaces
                std::replace(item.begin(), item.end(), C_BACK_SLASH, C_SPACE);
                std::stringstream itemStream(item);

                // The face vertex is simply the first element in the itemStream
                int64_t faceVertex;
                itemStream >> faceVertex;

                // Add the face vertex and proceed to the next item
                // NOTE: This step will ignore any normals and texture data
                faceVertices.push_back(faceVertex);
            }

            // Three vertices, and one triangle
            if (faceVertices.size() == 3)
            {
                Triangle triangle;
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);
            }

            // Four vertices, and two triangles
            else if (faceVertices.size() == 4)
            {
                Triangle triangle;

                // First triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);

                // Second triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[3] - 1;
                triangle[2] = faceVertices[0] - 1;
                triangles.push_back(triangle);
            }

            // Five vertices, and three triangles
            else if (faceVertices.size() == 5)
            {
                Triangle triangle;

                // First triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);

                // Second triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[3] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Third triangle
                triangle[0] = faceVertices[4] - 1;
                triangle[1] = faceVertices[0] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);
            }

            // Six vertices, and four triangles
            else if (faceVertices.size() == 6)
            {
                Triangle triangle;

                // First triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);

                // Second triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[3] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Third triangle
                triangle[0] = faceVertices[4] - 1;
                triangle[1] = faceVertices[5] - 1;
                triangle[2] = faceVertices[0] - 1;
                triangles.push_back(triangle);

                // Fourth triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[4] - 1;
                triangle[2] = faceVertices[0] - 1;
                triangles.push_back(triangle);
            }

            // Seven vertices, and five triangles
            else if (faceVertices.size() == 7)
            {
                Triangle triangle;

                // First triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);

                // Second triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[3] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Third triangle
                triangle[0] = faceVertices[4] - 1;
                triangle[1] = faceVertices[5] - 1;
                triangle[2] = faceVertices[6] - 1;
                triangles.push_back(triangle);

                // Fourth triangle
                triangle[0] = faceVertices[5] - 1;
                triangle[1] = faceVertices[0] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Fifth triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[2] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);
            }

            // Eight vertices, and six triangles
            else if (faceVertices.size() == 8)
            {
                Triangle triangle;

                // First triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[1] - 1;
                triangle[2] = faceVertices[2] - 1;
                triangles.push_back(triangle);

                // Second triangle
                triangle[0] = faceVertices[2] - 1;
                triangle[1] = faceVertices[3] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Third triangle
                triangle[0] = faceVertices[4] - 1;
                triangle[1] = faceVertices[5] - 1;
                triangle[2] = faceVertices[6] - 1;
                triangles.push_back(triangle);

                // Fourth triangle
                triangle[0] = faceVertices[6] - 1;
                triangle[1] = faceVertices[7] - 1;
                triangle[2] = faceVertices[0] - 1;
                triangles.push_back(triangle);

                // Fifth triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[2] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);

                // Sixth triangle
                triangle[0] = faceVertices[0] - 1;
                triangle[1] = faceVertices[6] - 1;
                triangle[2] = faceVertices[4] - 1;
                triangles.push_back(triangle);
            }

            // Otherwise, cannot handle
            else
            {
                LOG_ERROR("\nThe mesh [%s] has faces with N-gons (N > 8)!", filePath.c_str());
            }
        }

        // Texture, just ignore it
        else if (token == OBJ_TEXTURE_FLAG)
            continue;

        ++progress;
    }
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    file.close();
}

void importPLY(const std::string &filePath, Vertices& vertices, Triangles& triangles,
               const bool& verbose)
{
    // Start the timer
    TIMER_SET;

    std::ifstream file;
    file.open(filePath.c_str());

    // Line by line
    std::string line;
    while (true)
    {
        std::getline(file, line);
        if (std::string::npos != line.find(PLY_VERTEX_FLAG))
            break;
    }

    // Tokenize
    std::string token;

    // Start streaming the file
    std::stringstream stream(line);
    stream >> token >> token;
    size_t numberVertices;
    stream >> numberVertices;

    // Textures are not always there
    bool hasTexture = false;

    // Get all the elements
    while (true)
    {
        // Faces
        std::getline(file, line);
        if (std::string::npos != line.find(PLY_FACE_FLAG))
            break;

        // Textures
        if (std::string::npos != line.find(PLY_TEXTURE_FLAG))
            hasTexture = true;
    }

    std::stringstream lineStream(line);
    lineStream >> token >> token;
    size_t numberTriangles;
    lineStream >> numberTriangles;
    while (true)
    {
        // Done
        std::getline(file, line);
        if (std::string::npos != line.find(PLY_END_FLAG))
            break;
    }

    // Resize the vertex array
    vertices.resize(numberVertices);

    // Resize the texture array
    Textures textures;
    if (hasTexture)
        textures.resize(numberVertices);

    // Reading the vertices
    if (verbose) LOOP_STARTS("Loading Vertices");
    for (size_t i = 0; i < numberVertices; ++i)
    {
        if (verbose) LOOP_PROGRESS_FRACTION(i, numberVertices);

        std::getline(file, line);
        vertices[i] = Parsers::parseVector3f(line);

        // Texture included
        if (hasTexture)
        {
            std::getline(file, line);
            SimpleTexture texture = Parsers::parseVector2f(line);
            textures[i] = texture;
            textures[i][1] = 1 - textures[i][1];
        }
    }
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);


    // Read the triangles
    if (verbose) LOOP_STARTS("Loading Triangles");
    for (size_t i = 0; i < numberTriangles; ++i)
    {
        if (verbose) LOOP_PROGRESS_FRACTION(i, numberTriangles);

        std::getline(file, line);
        Parsers::parseFaces(line, triangles);
    }
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    // Statistics
    if (verbose) LOG_STATUS_IMPORTANT("Importing Stats.");
    if (verbose) LOG_STATS(GET_TIME_SECONDS);


    // Close the file
    file.close();
}

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
 * foreach triangle                      - 50 bytes:
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
void importSTL(const std::string &filePath, Vertices& vertices, Triangles& triangles,
               const bool &verbose)
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
    if (strcmp(keyWord, STL_SOLID_KEYWORD.c_str()))
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
            if (strcmp(keyWord, STL_FACET_KEYWORD.c_str()))
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
        if (fread(&numberElements, 4, 1, filePointer)) { /* NOTHING TO DO */ }

        size_t vertexIndex = 0;
        for (int i = 0; i < numberElements; ++i)
        {
            // Unexpected end of file
            char facet[50];
            if (!fread(facet, 50, 1, filePointer))
            {
                LOG_ERROR("Unexpected end of file [ %s ]!", filePath.c_str());
            }

            // Vertices
            Vertex v1((*((float *)(facet + 12))),
                      (*((float *)(facet + 16))),
                      (*((float *)(facet + 20))));

            Vertex v2((*((float *)(facet + 24))),
                      (*((float *)(facet + 28))),
                      (*((float *)(facet + 32))));

            Vertex v3((*((float *)(facet + 36))),
                      (*((float *)(facet + 40))),
                      (*((float *)(facet + 44))));

            // Construct the triangles with the vertices on the fly
            Triangle triangle;

            triangle[0] = vertexIndex;
            vertices.push_back(v1);
            vertexIndex++;

            triangle[1] = vertexIndex;
            vertices.push_back(v2);
            vertexIndex++;

            triangle[2] = vertexIndex;
            vertices.push_back(v3);
            vertexIndex++;

            // New triangle
            triangles.push_back(triangle);
        }
    }
    else
    {
        // Vertices flags
        bool v1Flag = false, v2Flag = false, v3Flag = false;

        // THe index that will be used to track the vertices.
        size_t vertexIndex = 0;

        // Proceed with ASCII encoding
        char *line;
        while ((line = File::readLineFromFile(filePointer, 0)) != nullptr)
        {
            // Geth the coordinates
            float x, y, z;
            sscanf(line,"%64s %f %f %f", keyWord, &x, &y, &z);

            // Normal
            if (strcmp(keyWord, STL_FACET_KEYWORD.c_str()) == 0)
            {
                char normalKeyWord[64] = "";
                sscanf(line,"%64s %64s %f %f %f", keyWord , normalKeyWord, &x, &y, &z);
            }
            else
            {
                // Vertex
                if (!strcmp(keyWord, STL_VERTEX_KEYWORD.c_str()))
                {
                    Vertex vertex(x, y, z);

                    if (!v1Flag)
                    {
                        v1Flag = true;
                        vertices.push_back(vertex);
                        vertexIndex++;
                    }
                    else if (!v2Flag)
                    {
                        v2Flag = true;
                        vertices.push_back(vertex);
                        vertexIndex++;
                    }
                    else if (!v3Flag)
                    {
                        v3Flag = true;
                        vertices.push_back(vertex);
                        vertexIndex++;

                        // Now all vertices of the triangle are reconstructed, building triangle
                        Triangle triangle;
                        triangle[0] = vertexIndex - 1;
                        triangle[1] = vertexIndex - 2;
                        triangle[2] = vertexIndex - 3;
                        triangles.push_back(triangle);

                        // Reset all the vertices
                        v1Flag = false;
                        v2Flag = false;
                        v3Flag = false;
                    }
                }
            }
        }
    }

    // Close the file
    fclose(filePointer);
}

void importOFF(const std::string &filePath, Vertices& vertices, Triangles& triangles,
               const bool &verbose)
{
    // Start the timer
    TIMER_SET;

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
    if (sscanf(line, "%" PRIu64 " %" PRIu64 " %" PRIu64 "",
               &numberVertices, &numberTriangles, &numberEdges) < 3)
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
    TIMER_RESET;
    if (verbose) LOOP_STARTS("Loading Vertices");
    for (size_t i = 0; i < numberVertices; ++i)
    {
        LOOP_PROGRESS(i, numberVertices);

        // Vertex
        float x, y, z;
        if (fscanf(filePointer, "%f %f %f", &x, &y, &z) == 3)
        {
            vertices.push_back(Vertex(x, y, z));
        }
        else
        {
            LOG_ERROR("Could not read coordinates for vertex [ %" PRIu64 "]", i);
        }
    }
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);

    File::skipCommentAndBlankLines(filePointer);

    // Triangles data
    bool triangulate = false;
    TIMER_RESET;
    if (verbose) LOOP_STARTS("Loading Triangles");
    for (size_t i = 0; i < numberTriangles; ++i)
    {
         LOOP_PROGRESS(i, numberTriangles);

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
                    Triangle triangle(i1, i2, i3);
                    triangles.push_back(triangle);
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
    if (verbose) LOOP_DONE;
    if (verbose) LOG_STATS(GET_TIME_SECONDS);
}

#ifdef ULTRALISER_USE_H5
Vertices extractVertexList(H5::H5File* h5File)
{
    // Read the vertices data set
    H5::DataSet dataset = h5File->openDataSet("vertices");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataSpace = dataset.getSpace();

    // Dataset dimenions
    hsize_t dimensions[2];
    dataSpace.getSimpleExtentDims(dimensions, nullptr);

    // Resize the vector to contain the data
    Vertices vertices;
    vertices.resize(dimensions[0]);

    // Read the data
    dataset.read(vertices.data(), H5::PredType::NATIVE_FLOAT);

    // Close the dataset
    dataset.close();

    // Return the result
    return vertices;
}

Triangles extractTriangleList(H5::H5File* h5File)
{
    // Read the faces data set
    H5::DataSet dataset = h5File->openDataSet("faces");

    // Get its data-space to be able to access its dimensionality
    H5::DataSpace dataSpace = dataset.getSpace();

    // Dataset dimenions
    hsize_t dimensions[2];
    dataSpace.getSimpleExtentDims(dimensions, nullptr);

    // Resize the vector to contain the data
    Triangles triangles;
    triangles.resize(dimensions[0]);

    // Read the data
    dataset.read(triangles.data(), H5::PredType::NATIVE_INT64);

    // Close the dataset
    dataset.close();

    // Return the result
    return triangles;
}
#endif


void importH5(const std::string &filePath, Vertices& vertices, Triangles& triangles,
              const bool &verbose)
{
#ifdef ULTRALISER_USE_H5

    H5::H5File* h5File = new H5::H5File(filePath, H5F_ACC_RDONLY);

    // Get the vertex list
    vertices = extractVertexList(h5File);

    // Get the triangle list
    triangles = extractTriangleList(h5File);

#else
    LOG_ERROR("Cannot load the file %s, HDF5 library is missing!", filePath.c_str());
#endif
}


void exportOBJ(const std::string &prefix,
               const Vertex *vertices,
               const size_t &numberVertices,
               const Triangle* triangles,
               const size_t &numberTriangles)
{
    // Open the file
    std::string fileName = prefix + OBJ_EXTENSION;
    std::ofstream stream(fileName.c_str());
    if (!stream.good())
    {
        LOG_ERROR("Cannot write mesh file [ %s ]", fileName.c_str());
    }

    LOG_STATUS("Exporting OBJ Mesh : [ %s ]", fileName.c_str());

    // Start the time
    TIMER_SET;

    // Write the vertices
    LOOP_STARTS("Writing Vertices");
    for (size_t i = 0; i < numberVertices; ++i)
    {
        LOOP_PROGRESS_FRACTION(i, numberVertices);

        stream << OBJ_VERTEX_FLAG << SPACE
               << vertices[i][0] << SPACE
               << vertices[i][1] << SPACE
               << vertices[i][2] << NEW_LINE;
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Reset timer
    TIMER_RESET;

    LOOP_STARTS("Writing Triangles");
    for (size_t i = 0; i < numberTriangles; ++i)
    {
        LOOP_PROGRESS_FRACTION(i, numberTriangles);
        stream << OBJ_FACE_FLAG << SPACE
               << triangles[i][0] + 1 << SPACE
               << triangles[i][1] + 1 << SPACE
               << triangles[i][2] + 1 << NEW_LINE;
    }
    LOOP_DONE;

    // Statistics
    LOG_STATS(GET_TIME_SECONDS);

    // Add the end flags
    stream<<OBJ_END_FLAG;

    // Close the file stream
    stream.close();
}

void exportOFF(const std::string &prefix,
               const Vertex *vertices,
               const size_t &numberVertices,
               const Triangle* triangles,
               const size_t &numberTriangles)
{
    // Open a stream to write the data
    std::string fileName = prefix + OFF_EXTENSION;

    std::ofstream outputStream(fileName.c_str());
    if (!outputStream.good())
    {
        LOG_ERROR("Cannot write mesh file [ %s ]", fileName.c_str());
    }

    LOG_STATUS("Exporting OFF Mesh : [ %s ]", fileName.c_str());

    // Header
    outputStream << "OFF" << std::endl;

    // Counts
    outputStream << numberVertices << SPACE
                 << numberTriangles << SPACE
                 << numberVertices  + numberTriangles - 2 << std::endl;

    // Write the vertices, with scientific precision
    LOOP_STARTS("Writing Vertices");
    outputStream << std::scientific;
    outputStream << std::setprecision(9);
    TIMER_SET;
    for (size_t i = 0; i < numberVertices; i++)
    {
        LOOP_PROGRESS_FRACTION(i, numberVertices);

        outputStream << vertices[i].x() << SPACE
                     << vertices[i].y() << SPACE
                     << vertices[i].z() << NEW_LINE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Write the indices
    LOOP_STARTS("Writing Triangles");
    TIMER_RESET;
    for (size_t i = 0; i < numberTriangles; i++)
    {
        LOOP_PROGRESS(i, numberTriangles);

        // Preprend with 3 for the edges (triangle and not quad mesh)
        outputStream << "3 "
                     << triangles[i][0] << SPACE
                     << triangles[i][1] << SPACE
                     << triangles[i][2] << NEW_LINE;
    }
    LOOP_DONE;
    LOG_STATS(GET_TIME_SECONDS);

    // Close the stream
    outputStream.close();
}

void exportSTL(const std::string &prefix,
               const Vertex *vertices,
               const size_t &,
               const Triangle* triangles,
               const size_t &numberTriangles)
{
    // Open a stream to write the data
    std::string fileName = prefix + STL_EXTENSION;

    LOG_STATUS("Exporting STL Mesh : [ %s ]", fileName.c_str());

    // If the path does not exist, just set a warning and return
    FILE *filePointer;
    if ((filePointer = fopen(fileName.c_str(), "w")) == nullptr)
    {
        LOG_WARNING("Can NOT write the mesh to the following file [ %s ]!", fileName.c_str());
        return;
    }

    // Header
    fprintf(filePointer, "solid Ultraliser");
    fprintf(filePointer, "\n");

    // Write the indices
    LOOP_STARTS("Writing Vertices and Triangles");
    TIMER_SET;
    for (size_t i = 0; i < numberTriangles; ++i)
    {
        LOOP_PROGRESS(i, numberTriangles);

        // Get the triangle
        Triangle triangle = triangles[i];

        // Get the vertices
        Vertex v0 = vertices[triangle[0]];
        Vertex v1 = vertices[triangle[1]];
        Vertex v2 = vertices[triangle[2]];

        // Compute the normal
        Vector3f normal = computeNormal(v0, v1, v2);

        // Face normal
        fprintf(filePointer, " facet normal %f %f %f", normal.x(), normal.y(), normal.z());
        fprintf(filePointer, "\n");

        // Outer loop
        fprintf(filePointer, "  outer loop");
        fprintf(filePointer, "\n");

        // V0
        fprintf(filePointer, "   vertex %f %f %f", v0.x(), v0.y(), v0.z());
        fprintf(filePointer, "\n");

        // V1
        fprintf(filePointer, "   vertex %f %f %f", v1.x(), v1.y(), v1.z());
        fprintf(filePointer, "\n");

        // V2
        fprintf(filePointer, "   vertex %f %f %f", v2.x(), v2.y(), v2.z());
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

void exportPLY(const std::string &prefix,
               const Vertex *vertices,
               const size_t &numberVertices,
               const Triangle* triangles,
               const size_t &numberTriangles,
               bool writeASCII)
{
    // Start the timer
    TIMER_SET;

    // Open a stream to write the data
    std::string filePath = prefix + PLY_EXTENSION;

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
    fprintf(filePointer, "element vertex %" PRIu64 "", numberVertices);
    fprintf(filePointer, "\n");

    fprintf(filePointer, "property float x");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property float y");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property float z");
    fprintf(filePointer, "\n");

    // Triangle count
    fprintf(filePointer, "element face %" PRIu64 "", numberTriangles);
    fprintf(filePointer, "\n");
    fprintf(filePointer, "property list uchar int vertex_indices");
    fprintf(filePointer, "\n");
    fprintf(filePointer, "end_header");
    fprintf(filePointer, "\n");

    size_t progress = 0;
    if (writeASCII)
    {
        LOOP_STARTS("Writing Vertices");
        for (size_t i = 0; i < numberVertices; i++)
        {
            LOOP_PROGRESS_FRACTION(i, numberVertices);
            fprintf(filePointer, "%f %f %f", vertices[i].x(), vertices[i].y(), vertices[i].z());
            fprintf(filePointer, "\n");
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        LOOP_STARTS("Writing Vertices");
        for (size_t i = 0; i < numberVertices; i++)
        {
            LOOP_PROGRESS_FRACTION(i, numberVertices);

            float floatCoordinates[3];
            floatCoordinates[0] = vertices[i].x();
            floatCoordinates[1] = vertices[i].y();
            floatCoordinates[2] = vertices[i].z();
            fwrite(floatCoordinates, sizeof(float), 3, filePointer);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Writing the triangles
    if (writeASCII)
    {
        progress = 0;
        TIMER_RESET;
        LOOP_STARTS("Writing Triangles");
        for (size_t i = 0; i < numberTriangles; i++)
        {
            LOOP_PROGRESS(i, numberTriangles);

            // Preprend with 3 for the edges (triangle and not quad mesh)
            fprintf(filePointer, "3 %ld %ld %ld", triangles[i][0], triangles[i][1], triangles[i][2]);
            fprintf(filePointer, "\n");
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }
    else
    {
        progress = 0;
        TIMER_RESET;
        LOOP_STARTS("Writing Triangles");
        for (size_t i = 0; i < numberTriangles; i++)
        {
            LOOP_PROGRESS(i, numberTriangles);

            int triangleVertexIndex[3];
            triangleVertexIndex[0] = triangles[i][0];
            triangleVertexIndex[1] = triangles[i][1];
            triangleVertexIndex[2] = triangles[i][2];

            const int elementCount = 3;
            fwrite(&elementCount, sizeof(unsigned char), 1, filePointer);
            fwrite(triangleVertexIndex, sizeof(int), 3, filePointer);
        }
        LOOP_DONE;
        LOG_STATS(GET_TIME_SECONDS);
    }

    // Close the file
    fclose(filePointer);
}

}
