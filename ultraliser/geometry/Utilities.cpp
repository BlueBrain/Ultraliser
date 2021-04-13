/***************************************************************************************************
 * Copyright (c) 2016 - 2021
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
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

#include <geometry/Utilities.h>
#include <utilities/File.h>
#include <common/Logging.h>

namespace Ultraliser
{
namespace Utilities
{

/**
 * @brief inRange
 * @param value
 * @param lowerBound
 * @param upperBound
 * @return
 */
template < class T >
bool inRange(T value, T lowerBound, T upperBound)
{
    return (value >= lowerBound && value <= upperBound);
}

template bool inRange(int8_t, int8_t, int8_t);
template bool inRange(int16_t, int16_t, int16_t);
template bool inRange(int32_t, int32_t, int32_t);
template bool inRange(int64_t, int64_t, int64_t);
template bool inRange(uint8_t, uint8_t, uint8_t);
template bool inRange(uint16_t, uint16_t, uint16_t);
template bool inRange(uint32_t, uint32_t, uint32_t);
template bool inRange(uint64_t, uint64_t, uint64_t);
template bool inRange(float, float, float);
template bool inRange(double, double, double);

/**
 * @brief computeMaximumBoundingBox
 * @param boundingBoxList
 * @param pMin
 * @param pMax
 */
void computeMaximumBoundingBox(std::vector< NeuronBoundingBox> boundingBoxList,
                               Vector3f& pMin, Vector3f& pMax)
{
    for (size_t i = 0; i < boundingBoxList.size(); ++i)
    {
        Vector3f pMinNeuron = boundingBoxList[i].pMin;
        Vector3f pMaxNeuron = boundingBoxList[i].pMax;

        for (int32_t iDim = 0; iDim < DIMENSIONS; iDim++)
        {
            if (pMinNeuron[iDim] < pMin[iDim])
                pMin[iDim] = pMinNeuron[iDim];

            if (pMaxNeuron[iDim] > pMax[iDim])
                pMax[iDim] = pMaxNeuron[iDim];
        }
    }
}

/**
 * @brief BBox
 * @param points
 * @param pMin
 * @param pMax
 * @param initialized
 */
void BBox(const std::vector< Vector3f >& points,
          Vector3f& pMin, Vector3f& pMax, const bool initialized)
{
    if (!initialized) { pMin = points[0]; pMax = points[0]; }

    for (size_t i = 1 ; i < points.size(); ++i)
    {
        for (int32_t iDim = 0; iDim < DIMENSIONS; iDim++)
        {
            if (points[i][iDim] < pMin[iDim])
                pMin[iDim] = points[i][iDim];

            if (points[i][iDim] > pMax[iDim])
                pMax[iDim] = points[i][iDim];
        }
    }
}

//void BBox(const OriginalMesh& mesh, Vector3f& pMin, Vector3f& pMax)
//{
//    BBox(mesh.vertices, pMin, pMax);
//}

bool pInside(const Vector3f& pMin, const Vector3f& pMax, const Vector3f& p)
{
    for (int32_t dim = 0; dim < DIMENSIONS; dim++)
        if (p[dim] <  pMin[dim] || p[dim] > pMax[dim])
            return false;
    return true;
}

//void makeCube(OriginalMesh& mesh, const Vector3f& pMin, const Vector3f& pMax)
//{
//    Vector3f diagonal = pMax - pMin;
//    Matrix3f transform = Matrix3f::identity();

//    for (int32_t i = 0; i < DIMENSIONS; ++i)
//        transform(i, i) = diagonal[i];

//    mesh = UNIT_CUBE;
//    for (size_t i = 0; i < mesh.vertices.size(); ++i)
//        mesh.vertices[i] = pMin + transform * mesh.vertices[i];
//}

bool isNbr(const Vec3i_64& a, const Vec3i_64& b, int64_t vert)
{
    for (int32_t i = 0; i < DIMENSIONS; ++i)
    {
        int64_t va = a[i];
        if (va <= vert)
            continue;

        for (int32_t j = 0; j < 3; j++)
        {
            int64_t vb = b[j];
            if (vb <= vert)
                continue;

            if (va == vb)
                return true;
        }
    }

    return false;
}

//void adjustList(const OriginalMesh& mesh, std::vector< std::vector< int32_t > >& adjMat)
//{
//    if (adjMat.size() == mesh.triangles.size())
//        return;

//    std::vector< std::vector< int64_t > > trigList;
//    trigList.resize(mesh.vertices.size());

//    for (uint64_t i = 0; i < mesh.triangles.size(); ++i)
//    {
//        for (uint64_t j = 0; j < DIMENSIONS; j++)
//        {
//            uint64_t vertexIdx = I2UI64(mesh.triangles[I2UI64(i)][I2I32(j)]);
//            trigList[I2UI64(vertexIdx)].push_back(I2I64(i));
//        }
//    }

//    adjMat.resize(mesh.triangles.size());
//    for (uint64_t i = 0; i < mesh.vertices.size(); ++i)
//    {
//        uint64_t nNumber = trigList[i].size();

//        for (uint64_t j = 0; j < nNumber; j++)
//        {
//            int64_t tj = trigList[I2UI64(i)][I2UI64(j)];
//            for (uint64_t k = (j + 1); k < nNumber; k++)
//            {
//                int64_t tk = trigList[i][k];
//                if (isNbr(mesh.triangles[I2UI64(tj)],
//                          mesh.triangles[I2UI64(tk)], I2I64(i)))
//                {
//                    adjMat[I2UI64(tj)].push_back(I2I32(tk));
//                    adjMat[I2UI64(tk)].push_back(I2I32(tj));
//                }
//            }
//        }
//    }
//}

void printBoundingBoxData(const Vector3f &pMin, const Vector3f &pMax,
                          std::string prefix)
{
    // Center the bounding box
    Vector3f bound = pMax - pMin;
    Vector3f pMinCenter = -bound / 2.0;
    Vector3f pMaxCenter = bound / 2.0;

    //    std::cout << "Center BB" << std::endl;
    //    std::cout << "\t"; pMinCenter.print();
    //    std::cout << "\t"; pMaxCenter.print();

    std::string fileName = prefix + std::string(BOUNDS_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    header << "Original BB:" << std::endl;
    header << "\t" << pMin.x() << " " << pMin.y() << " " << pMin.z()
           << std::endl;
    header << "\t" << pMax.x() << " " << pMax.y() << " " << pMax.z()
           << std::endl;

    header << "Center BB:" << std::endl;
    header << "\t"
           << pMinCenter.x() << " " << pMinCenter.y() << " " << pMinCenter.z()
           << std::endl;
    header << "\t"
           << pMaxCenter.x() << " " << pMaxCenter.y() << " " << pMaxCenter.z()
           << std::endl;
    header.close();
}

/**
 * @brief printVoxelData
 * @param pMin
 * @param pMax
 * @param volumeSize
 * @param prefix
 * @return
 */
void printVoxelData(const Vector3f &pMin, const Vector3f &pMax,
                    const Vector3f volumeSize,
                    std::string prefix)
{
    Vector3f voxelSize = (pMax - pMin) / volumeSize;

    std::string fileName = prefix + std::string(VOXEL_EXTENSION);
    std::fstream header;
    header.open(fileName.c_str(), std::ios::out);

    header << "Voxel Size:" << std::endl;
    header << voxelSize.x() << " "
           << voxelSize.y() << " "
           << voxelSize.z() << std::endl;
    header.close();
}

///**
// * @brief getBoundingBox
// * @param options
// * @param neurons
// * @param pMax
// * @param pMin
// */
//void getBoundingBox(const NeuroVoxyOptions* options, Neurons neurons,
//                     Vector3f& pMax, Vector3f& pMin)
//{
//    if (options->boundsFile == EMPTY)
//    {
//        printf("* The bounding box will be computed on-the-fly \n");

//        // Vectors containing all the pMin and pMax of all the neurons
//        std::vector< Vector3f > pMinVector, pMaxVector;

//        pMinVector.resize(neurons.size());
//        pMaxVector.resize(neurons.size());

//        /// Compute the bounding box of the entire volume that contains all the
//        /// loaded meshes
//        size_t loadedMeshCount = 0;
//        printf("* Loading Meshes ... \n");
//        #pragma omp parallel for schedule(dynamic, 1)
//        for (size_t iNeuron = 0; iNeuron < neurons.size(); iNeuron++)
//        {
//            // Create and load the mesh from the file
//            std::stringstream meshNameStream;
//            meshNameStream << "neuron_" <<  neurons[iNeuron].gid << ".ply";
//            std::string meshName = meshNameStream.str();
//            std::string meshFile = meshName;
//            if (options->inputDirectory != EMPTY)
//                meshFile = options->inputDirectory + "/" + meshName;

//            if (Utilities::fileExists(meshFile))
//            {
//                #pragma omp atomic
//                loadedMeshCount++;

//                // Load the mesh
//                Mesh* mesh = new Mesh(meshFile, meshName, false, true);

//                // Transform the neuron if required
//                if (options->transformNeuron)
//                {
//                    mesh->transform(neurons[iNeuron].localToGlobalTransform);
//                }

//                // Compute its bounding box
//                Vector3f pMinMesh, pMaxMesh;
//                mesh->computeBoundingBox(pMinMesh, pMaxMesh);
//                pMinVector[iNeuron] = pMinMesh;
//                pMaxVector[iNeuron] = pMaxMesh;

//                // Free the mesh
//                mesh->~Mesh();
//            }
//            else
//            {
//                std::cout << "Missing: " << meshFile << std::endl;
//            }
//        }
//        printf("** [%zu/%zu] meshes were loaded \n",
//               loadedMeshCount, neurons.size());

//        if (loadedMeshCount == 0)
//            LOG_ERROR("Zero loaded meshes");

//        // Compute the bounding box of the neuron group
//        printf("* Computing Bounding Box \n");
//        pMax[0] = std::numeric_limits< float >::min();
//        pMax[1] = std::numeric_limits< float >::min();
//        pMax[2] = std::numeric_limits< float >::min();

//        pMin[0] = std::numeric_limits< float >::max();
//        pMin[1] = std::numeric_limits< float >::max();
//        pMin[2] = std::numeric_limits< float >::max();

//        for (uint64_t iNeuron = 0; iNeuron < neurons.size(); iNeuron++)
//        {
//            Vector3f pMinNeuron = pMinVector[iNeuron];
//            Vector3f pMaxNeuron = pMaxVector[iNeuron];

//            if (pMinNeuron.x() < pMin.x()) pMin.x() = pMinNeuron.x();
//            if (pMinNeuron.y() < pMin.y()) pMin.y() = pMinNeuron.y();
//            if (pMinNeuron.z() < pMin.z()) pMin.z() = pMinNeuron.z();

//            if (pMaxNeuron.x() > pMax.x()) pMax.x() = pMaxNeuron.x();
//            if (pMaxNeuron.y() > pMax.y()) pMax.y() = pMaxNeuron.y();
//            if (pMaxNeuron.z() > pMax.z()) pMax.z() = pMaxNeuron.z();
//        }
//    }
//    else
//    {
//        printf("* Loading the bounding box from [ %s ] \n",
//                options->boundsFile.c_str());

//        // Verify the bounding box file
//        if (Utilities::fileExists(options->boundsFile))
//            Utilities::parseBoundsFile(options->boundsFile, pMin, pMax);
//        else
//            LOG_ERROR("No bounding box file is provided !");

//    }
//}

}
}

