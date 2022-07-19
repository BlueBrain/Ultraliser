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

#include "StiffnessMatrix.h"

namespace Ultraliser
{
namespace Simulation
{

StiffnessMatrix::StiffnessMatrix(Nodes& nodes, 
                                 Tetrahedra& tetrahedra, 
                                 float stiffness, 
                                 float poissonRatio,
                                 float dt)
{
    _compute(nodes, tetrahedra, stiffness, poissonRatio, dt);
}

StiffnessMatrix::~StiffnessMatrix(){}

void StiffnessMatrix::_compute(Nodes& nodes, 
                          Tetrahedra& tetrahedra, 
                          float stiffness, 
                          float poissonRatio,
                          float dt)
{
    float d = stiffness / ((1.0f + poissonRatio) * (1.0f - 2.0f * poissonRatio));
    float d0 = d * (1.0f - poissonRatio);
    float d1 = d * poissonRatio;
    float d2 = d * (1.0f - 2.0f * poissonRatio) * 0.5f;
    float dt2 = dt * dt;

    // Compute strain matrix
    Eigen::MatrixXf D(6, 6);
    D << d0, d1, d1, .0f, .0f, .0f,
         d1, d0, d1, .0f, .0f, .0f,
         d1, d1, d0, .0f, .0f, .0f,
         .0f, .0f, .0f, d2, .0f, .0f,
         .0f, .0f, .0f, .0f, d2, .0f,
         .0f, .0f, .0f, .0f, .0f, d2;

    _KMatrices kMatrices(tetrahedra.size());
    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < nodes.size(); ++i) nodes[i]->index = i;

    // Compute tetrahedra stiffness matrices
    OMP_PARALLEL_FOR
    for (uint64_t i = 0; i < tetrahedra.size(); ++i)
        _computeTetrahedron(tetrahedra[i], kMatrices[i], D);

    // Compute stiffness matrices triplets
    std::vector<Eigen::Triplet<float>> kTriplets;
    std::vector<Eigen::Triplet<float>> aTriplets;
    _computeTriplets(tetrahedra, kMatrices, dt2, kTriplets, aTriplets);
    for (uint64_t i = 0; i < nodes.size(); ++i)
        _addIdentityValueToTriplets(i, nodes[i]->mass, aTriplets);

    // Compute stiffness matrices
    stiffnessMatrix.resize(nodes.size() * 3, nodes.size() * 3);
    _aMatrix.resize(nodes.size() * 3, nodes.size() * 3);
    stiffnessMatrix.resizeNonZeros(kTriplets.size());
    _aMatrix.resizeNonZeros(aTriplets.size());
    stiffnessMatrix.setFromTriplets(kTriplets.begin(), kTriplets.end());
    _aMatrix.setFromTriplets(aTriplets.begin(), aTriplets.end());

    kTriplets.clear();
    aTriplets.clear();
    matrixSolver.compute(_aMatrix);
}

void StiffnessMatrix::_computeTetrahedron(TetrahedronPtr tetrahedron, 
                                          _KMatrix& kMatrix, 
                                          Eigen::MatrixXf D)
{
    Vector3f x0 = tetrahedron->node0->initPosition();
    Vector3f x1 = tetrahedron->node1->initPosition();
    Vector3f x2 = tetrahedron->node2->initPosition();
    Vector3f x3 = tetrahedron->node3->initPosition();

    Matrix3f basis(x1 - x0, x2 - x0, x3 - x0, true);
    basis = basis.inverse();
    Vector3f b1 = basis.getRow(0);
    Vector3f b2 = basis.getRow(1);
    Vector3f b3 = basis.getRow(2);
    Vector3f b0 = -b1 - b2 - b3;

    Eigen::MatrixXf B0(6, 3);
    B0 << b0[0], 0.0, 0.0, 0.0, b0[1], 0.0, 0.0, 0.0, b0[2], b0[1], b0[0], 0.0,
        0.0, b0[2], b0[1], b0[2], 0.0, b0[0];
    Eigen::MatrixXf B0T = B0.transpose();

    Eigen::MatrixXf B1(6, 3);
    B1 << b1[0], 0.0, 0.0, 0.0, b1[1], 0.0, 0.0, 0.0, b1[2], b1[1], b1[0], 0.0,
        0.0, b1[2], b1[1], b1[2], 0.0, b1[0];
    Eigen::MatrixXf B1T = B1.transpose();

    Eigen::MatrixXf B2(6, 3);
    B2 << b2[0], 0.0, 0.0, 0.0, b2[1], 0.0, 0.0, 0.0, b2[2], b2[1], b2[0], 0.0,
        0.0, b2[2], b2[1], b2[2], 0.0, b2[0];
    Eigen::MatrixXf B2T = B2.transpose();

    Eigen::MatrixXf B3(6, 3);
    B3 << b3[0], 0.0, 0.0, 0.0, b3[1], 0.0, 0.0, 0.0, b3[2], b3[1], b3[0], 0.0,
        0.0, b3[2], b3[1], b3[2], 0.0, b3[0];
    Eigen::MatrixXf B3T = B3.transpose();

    float volume = tetrahedron->initVolume();

    kMatrix.k00 = B0T * D * B0 * volume;
    kMatrix.k11 = B1T * D * B1 * volume;
    kMatrix.k22 = B2T * D * B2 * volume;
    kMatrix.k33 = B3T * D * B3 * volume;
    kMatrix.k01 = B0T * D * B1 * volume;
    kMatrix.k02 = B0T * D * B2 * volume;
    kMatrix.k03 = B0T * D * B3 * volume;
    kMatrix.k12 = B1T * D * B2 * volume;
    kMatrix.k13 = B1T * D * B3 * volume;
    kMatrix.k23 = B2T * D * B3 * volume;
}


void StiffnessMatrix::_computeTriplets(Tetrahedra& tetrahedra, 
                                       _KMatrices& kMatrices, 
                                       float dt2,
                                       std::vector<Eigen::Triplet<float>>& kTriplets, 
                                       std::vector<Eigen::Triplet<float>>& aTriplets)
{
    for (uint64_t i = 0; i < tetrahedra.size(); ++i)
    {
        auto tet = tetrahedra[i];
        uint64_t id0 = tet->node0->index;
        uint64_t id1 = tet->node1->index;
        uint64_t id2 = tet->node2->index;
        uint64_t id3 = tet->node3->index;

        // row 0
        Eigen::Matrix3f k = kMatrices[i].k00;
        _addMatrixToTriplets(id0, id0, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id0, id0, k, aTriplets);
        k = kMatrices[i].k01;
        _addMatrixToTriplets(id0, id1, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id0, id1, k, aTriplets);
        k = kMatrices[i].k02;
        _addMatrixToTriplets(id0, id2, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id0, id2, k, aTriplets);
        k = kMatrices[i].k03;
        _addMatrixToTriplets(id0, id3, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id0, id3, k, aTriplets);
        // row 1
        k = kMatrices[i].k01.transpose();
        _addMatrixToTriplets(id1, id0, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id1, id0, k, aTriplets);
        k = kMatrices[i].k11;
        _addMatrixToTriplets(id1, id1, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id1, id1, k, aTriplets);
        k = kMatrices[i].k12;
        _addMatrixToTriplets(id1, id2, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id1, id2, k, aTriplets);
        k = kMatrices[i].k13;
        _addMatrixToTriplets(id1, id3, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id1, id3, k, aTriplets);
        // row 2
        k = kMatrices[i].k02.transpose();
        _addMatrixToTriplets(id2, id0, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id2, id0, k, aTriplets);
        k = kMatrices[i].k12.transpose();
        _addMatrixToTriplets(id2, id1, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id2, id1, k, aTriplets);
        k = kMatrices[i].k22;
        _addMatrixToTriplets(id2, id2, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id2, id2, k, aTriplets);
        k = kMatrices[i].k23;
        _addMatrixToTriplets(id2, id3, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id2, id3, k, aTriplets);
        // row 3
        k = kMatrices[i].k03.transpose();
        _addMatrixToTriplets(id3, id0, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id3, id0, k, aTriplets);
        k = kMatrices[i].k13.transpose();
        _addMatrixToTriplets(id3, id1, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id3, id1, k, aTriplets);
        k = kMatrices[i].k23.transpose();
        _addMatrixToTriplets(id3, id2, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id3, id2, k, aTriplets);
        k = kMatrices[i].k33;
        _addMatrixToTriplets(id3, id3, k, kTriplets);
        k *= dt2;
        _addMatrixToTriplets(id3, id3, k, aTriplets);
    }
}

void StiffnessMatrix::_addMatrixToTriplets(uint64_t rowIndex,
                                           uint64_t colIndex,
                                           const Eigen::Matrix3f& matrix,
                                           std::vector<Eigen::Triplet<float>>& triplets)
{
    rowIndex *= 3;
    colIndex *= 3;
    for (uint64_t i = 0; i < 3; ++i)
        for (uint64_t j = 0; j < 3; ++j)
            triplets.push_back(Eigen::Triplet<float>(rowIndex + i, colIndex + j, matrix(i, j)));
}

void StiffnessMatrix::_addIdentityValueToTriplets(uint64_t index,
                                                  double value,
                                                  std::vector<Eigen::Triplet<float>>& triplets)
{
    index *= 3;
    for (uint64_t i = 0; i < 3; ++i)
        triplets.push_back(Eigen::Triplet<float>(index + i, index + i, value));
}

void StiffnessMatrix::addVec3ToVec(uint64_t index,
                                   Vector3f& vec3,
                                   Eigen::VectorXf& vec)
{
    index *= 3;
    for (uint64_t i = 0; i < 3; ++i)
        vec[index + i] = vec3[i];
}

}
}
