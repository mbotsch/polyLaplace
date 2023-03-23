//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/Triangulation.h>
#include <unsupported/Eigen/NonLinearOptimization>
//=============================================================================

using namespace pmp;

//=============================================================================


void setup_stiffness_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S, int Laplace, int minpoint = 0);

void setup_mass_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, int Laplace, int minpoint = 0, bool lumped = true);

void setup_gradient(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G, int Laplace, int minpoint = 0);

void setup_divergence(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G, int Laplace, int minpoint = 0);

