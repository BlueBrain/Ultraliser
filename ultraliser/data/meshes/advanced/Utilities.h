#ifndef ULTRALISER_DATA_MESHES_ADVANCED_UTILITIES_H
#define ULTRALISER_DATA_MESHES_ADVANCED_UTILITIES_H

#include <data/meshes/advanced/primitives/AdvancedVertex.h>
#include <data/meshes/advanced/AdvancedMesh.h>

namespace Ultraliser
{

/**
 * @brief subSurfBetaLoop
 * @param k
 * @return
 */
double subSurfBetaLoop(int k);

/**
 * @brief loopRelaxOriginal
 * @param vertex
 */
void loopRelaxOriginal(AdvancedVertex *vertex);

/**
 * @brief remintsAppendCubeToList
 * @param inputTriangulation
 * @param inputList
 * @return
 */
bool remintsAppendCubeToList(AdvancedTriangle *inputTriangulation, List& inputList);
/**
 * @brief remintsIsVertexInCube
 * @param vertex
 * @param inputList
 * @return
 */
bool remintsIsVertexInCube(AdvancedVertex *vertex, List& inputList);

/**
 * @brief remintsSelectTrianglesInCubes
 * @param inputMesh
 */
void remintsSelectTrianglesInCubes(AdvancedMesh *inputMesh);

/**
 * @brief jitterIncrease
 * @param charArray
 */
void jitterIncrease(char *charArray);

/**
 * @brief jitterDecrease
 * @param charArray
 */
void jitterDecrease(char *charArray);

/**
 * @brief jitterCoordinate
 * @param coordinate
 * @param j
 */
void jitterCoordinate(double& coordinate, int j);

}

#endif // ULTRALISER_DATA_MESHES_ADVANCED_UTILITIES_H
