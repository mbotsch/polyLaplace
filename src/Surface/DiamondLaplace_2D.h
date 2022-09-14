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

//=============================================================================

using namespace pmp;

//=============================================================================

void setup_diamond_mass_matrix(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &M);

void compute_primal_points(SurfaceMesh &mesh, int minpoint);

void setup_diamond_gradient_divergence_intrinsic(
    SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G,
    Eigen::SparseMatrix<double> &D);

//----------------------------------------------------------------------------------
