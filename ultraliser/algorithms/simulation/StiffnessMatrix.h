/***************************************************************************************************
 * Copyright (c) 2016 - 2022
 * Blue Brain Project (BBP) / Ecole Polytechniqe Federale de Lausanne (EPFL)
 *
 * Author(s)
 *      Juan Jose Garcia Cantero <juanjose.garcia@epfl.ch>
 *      Marwan Abdellah <marwan.abdellah@epfl.ch >
 *
 * This file is part of Ultraliser < https://github.com/BlueBrain/Ultraliser >
 *
 * This library is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License version 3.0 as published by the
 * Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the GNU General Public License along with
 * this library; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place - Suite 330, Boston, MA 02111-1307, USA. You can also find it on the
 * GNU web site < https://www.gnu.org/licenses/gpl-3.0.en.html >
 **************************************************************************************************/

#pragma once

#include "Spring.h"
#include "Tetrahedron.h"

#include <Eigen/Sparse>

namespace Ultraliser
{
namespace Simulation
{

/**
 * @brief The StiffnessMatrix class
 */
class StiffnessMatrix
{
public:

    /**
     * @brief StiffnessMatrix
     * Constructor
     *
     * @param nodes
     * Nodes conforming the stiffness matrix
     * @param tetrahedra
     * Tetrahedra conforming the stiffness matrix
     * @param stiffness
     * Young modulus
     * @param poissonRation
     * Poisson ratio
     * @param dt
     * Delta time
     */
    StiffnessMatrix(Nodes& nodes, 
                    Tetrahedra& tetrahedra, 
                    float stiffness = 1000.0f, 
                    float poissonRatio = 0.3f,
                    float dt = 0.01f);

    /**
     * @brief ~StiffnessMatrix
     * Destructor
     */
    ~StiffnessMatrix();
    
private:

    /**
     * @brief KMatrix
     */
    typedef struct _KMatrix
    {
        Eigen::Matrix3f k00;
        Eigen::Matrix3f k11;
        Eigen::Matrix3f k22;
        Eigen::Matrix3f k33;
        Eigen::Matrix3f k01;
        Eigen::Matrix3f k02;
        Eigen::Matrix3f k03;
        Eigen::Matrix3f k12;
        Eigen::Matrix3f k13;
        Eigen::Matrix3f k23;
    } _KMatrix;

    /**
     * @brief KMatrices
     */
    typedef std::vector<_KMatrix> _KMatrices;

    /**
     * @brief _compute
     * Compute the stifness matrices
     * 
     * @param nodes
     * Nodes conforming the stiffness matrix
     * @param tetrahedra
     * Tetrahedra conforming the stiffness matrix
     * @param stiffness
     * Young modulus
     * @param poissonRation
     * Poisson ratio
     * @param dt
     * Delta time
     */
    void _compute(Nodes& nodes, 
                  Tetrahedra& tetrahedra, 
                  float stiffness, 
                  float poissonRatio,
                  float dt);

    /**
     * @brief _computeTetrahedron
     * Compute tetrahedron stifness matrices
     * 
     * @param tetrahedron
     * Tetrahedron
     * @param kMatrix
     * Tetrahedron stiffness matrices
     * @param D
     * Strain matrix
     */
    void _computeTetrahedron(TetrahedronPtr tetrahedron, _KMatrix& kMatrix, Eigen::MatrixXf D);


    /**
     * @brief _computeTriplets
     * Compute tetrahedron stifness matrices
     * 
     * @param tetrahedra
     * Tetrahedra
     * @param kMatrices
     * Tetrahedra stiffness matrices
     * @param dt2
     * Squared delta time
     * @param kTriplets
     * stiffness matrix triplets 
     * @param aTriplets
     * stiffness and mass matrix triplets
     */
    void _computeTriplets(Tetrahedra& tetrahedra, 
                          _KMatrices& kMatrices, 
                          float dt2,
                          std::vector<Eigen::Triplet<float>>& kTriplets, 
                          std::vector<Eigen::Triplet<float>>& aTriplets);

    /**
     * @brief _addMatrixToTriplets
     * Add matrix to triplets
     * 
     * @param rowIndex
     * Row index
     * @param colIndex
     * Column index 
     * @param matrix
     * Matrix
     * @param triplets
     * Triplets 
     */
    static void _addMatrixToTriplets(size_t rowIndex,
                                     size_t colIndex,
                                     const Eigen::Matrix3f& matrix,
                                     std::vector<Eigen::Triplet<float>>& triplets);

    /**
     * @brief _addIdentityValueToTriplets
     * Add value to triplets formatted as identity matrix
     * 
     * @param index
     * Index
     * @param value
     * Value
     * @param triplets
     * Triplets 
     */
    static void _addIdentityValueToTriplets(size_t index,
                                            double value,
                                            std::vector<Eigen::Triplet<float>>& triplets);

public:

    /**
     * @brief _addVec3ToVec
     * Add Vector3f to Eigen::VectorXf 
     * 
     * @param index
     * Index
     * @param vec3
     * Vector3f
     * @param vec
     * Eigen::VectorXf 
     */
    static void addVec3ToVec(size_t index,
                             Vector3f& vec3,
                             Eigen::VectorXf& vec);

public:

    /**
     * @brief stiffnessMatrix
     */
    Eigen::SparseMatrix<float> stiffnessMatrix;
    
    /**
     * @brief aMatrixSolver
     */
    Eigen::ConjugateGradient<Eigen::SparseMatrix<float>> matrixSolver;

private:
    /**
     * @brief _aMatrix
     */
    Eigen::SparseMatrix<float> _aMatrix;

};

/**
 * @brief StiffnessMatrixPtr
 */
typedef StiffnessMatrix* StiffnessMatrixPtr;

}
}
