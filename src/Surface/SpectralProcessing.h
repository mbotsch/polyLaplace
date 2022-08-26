//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "AQAPoly_Laplacian.h"

//=============================================================================

using namespace pmp;

//=============================================================================



double solve_eigenvalue_problem(SurfaceMesh &mesh_, int laplace, int face_point, int degree = 2, CoarseDimension coarseningType = Edges,
                                std::string meshname = "default");

void analytic_eigenvalues_unitsphere(Eigen::VectorXd &eval, int n);


// spherical harmonic function at a point p for a band l and its range m
// see: http://silviojemma.com/public/papers/lighting/spherical-harmonic-lighting.pdf
double sphericalHarmonic(pmp::Point p, int l, int m);

// renormalisation constant for poisson_SH function
double scale(int l, int m);

// evaluate an Associated Legendre Polynomial P(l,m,x) at x
double legendre_Polynomial(int l, int m, double x);

double rmse_sh(SurfaceMesh &mesh, unsigned int laplace, unsigned int min_point_, bool lumped=true, int degree=1, CoarseDimension coarseningType =Edges);
//=============================================================================
