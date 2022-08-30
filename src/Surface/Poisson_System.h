
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "[AW11]Laplace.h"
#include "LaplaceConstruction.h"
#include "[dGBD20]Laplace.h"
//=============================================================================
using namespace pmp;


double solve_poisson_system(pmp::SurfaceMesh &mesh, int laplace, int minpoint,
                            int function,int l = 3, int m = -1);

double poisson_function(pmp::Point &p, int function);

double laplace_of_poisson_function(pmp::Point &p, int function);

double franke_function(double x, double y);

double laplace_franke_function(double x, double y);

double spherical_harmonic_function(double x, double y, double z);

double spherical_harmonic_function_scaled(double x, double y, double z);

void solve_laplace_equation(pmp::SurfaceMesh &mesh, int laplace, int face_point);
