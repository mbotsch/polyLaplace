//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#define PMP_SCALAR_TYPE_64

//=============================================================================

using namespace pmp;

//=============================================================================

void setup_diamond_mass_matrix(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &M);

void compute_primal_points(SurfaceMesh &mesh, int minpoint);

void setup_diamond_gradient_divergence_intrinsic(
    SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G,
    Eigen::SparseMatrix<double> &D);

void setup_diamond_gradient_divergence_intrinsic_mullification(
    SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G,
    Eigen::SparseMatrix<double> &D, double delta);

void setup_diamond_mass_matrix_intrinsic_mullification(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &M, double delta);
double mollification_epsilon( SurfaceMesh &mesh, double delta);


//----------------------------------------------------------------------------------
