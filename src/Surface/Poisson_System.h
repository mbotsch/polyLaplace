
//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
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

double solve_poisson_system(pmp::SurfaceMesh& mesh, int laplace, int minpoint,
                            int function, int& iterations, double& condition_number, int l = 4, int m = 2);

double poisson_function(pmp::Point& p, int function);

double laplace_of_poisson_function(pmp::Point& p, int function);

double franke_function(double x, double y);

double laplace_franke_function(double x, double y);
