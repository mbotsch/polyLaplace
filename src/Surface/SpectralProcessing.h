//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "[AW11]Laplace.h"

//=============================================================================

using namespace pmp;

//=============================================================================

double solve_eigenvalue_problem(SurfaceMesh& mesh_, int laplace, int face_point,
                                const std::string& meshname = "default");

void analytic_eigenvalues_unitsphere(Eigen::VectorXd& eval, int n);

// spherical harmonic function at a point p for a band l and its range m
// see: http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
double sphericalHarmonic(pmp::Point p, int l, int m);

// renormalisation constant for poisson_SH function
double scale(int l, int m);

// evaluate an Associated Legendre Polynomial P(l,m,x) at x
double legendre_Polynomial(int l, int m, double x);

double rmse_sh(SurfaceMesh& mesh, int laplace, int min_point_,
               bool lumped = true);

double condition_number(pmp::SurfaceMesh& mesh, int laplace, int minpoint,
                        Eigen::Vector3d& values, bool generalized = false);
double get_condition_number(const Eigen::SparseMatrix<double>& M, bool firstEigZero = false);
//=============================================================================
