
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "LaplaceConstruction.h"
#include "DisneyLaplace.h"
//=============================================================================
using namespace pmp;


double solve_poisson_system(pmp::SurfaceMesh &mesh, int laplace, int minpoint,
                            int function, int degree = 1,CoarseDimension coarseningType=Edges, int l = 3,
                            int m = -1);

double poisson_function(pmp::Point &p, int function);

double laplace_of_poisson_function(pmp::Point &p, int function);

double franke_function(double x, double y);

double laplace_franke_function(double x, double y);

double spherical_harmonic_function(double x, double y, double z);

double spherical_harmonic_function_scaled(double x, double y, double z);

double inverse_mean_edgelenth(pmp::SurfaceMesh &mesh);

double write_poisson_convergence_data_csv(int function);

void solve_laplace_equation(pmp::SurfaceMesh &mesh, int laplace, int face_point,
                            int degree,CoarseDimension coarseningType=Edges);

void cut_to_franke_Halfsphere(pmp::SurfaceMesh &mesh);

double laplace_franke_halfsphere_function(double x, double y);

double franke_halfsphere_function(double x, double y);

double condition_number(pmp::SurfaceMesh &mesh, int laplace, int minpoint,
                           int degree = 1,CoarseDimension coarseningType=Edges);