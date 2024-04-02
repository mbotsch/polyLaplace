//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/surface_mesh.h>
#include <pmp/algorithms/differential_geometry.h>
#include <unsupported/Eigen/NonLinearOptimization>
//=============================================================================

using namespace pmp;

//=============================================================================

void setup_stiffness_matrices(SurfaceMesh& mesh, Eigen::SparseMatrix<double>& S,
                              int Laplace, int minpoint = 0);

void setup_mass_matrices(SurfaceMesh& mesh, Eigen::SparseMatrix<double>& M,
                         int Laplace, int minpoint = 0, bool lumped = true);

void setup_gradient(SurfaceMesh& mesh, Eigen::SparseMatrix<double>& G,
                    int Laplace, int minpoint = 0);

void setup_divergence(SurfaceMesh& mesh, Eigen::SparseMatrix<double>& G,
                      int Laplace, int minpoint = 0);
