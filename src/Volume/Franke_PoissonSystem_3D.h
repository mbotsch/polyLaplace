
//=============================================================================

#pragma once

//=============================================================================

#include "../VolumeMesh/VolumeMesh.h"
#include "diffgeo_3D.h"

//=============================================================================

double franke3d(VolumeMesh::PointT vec);

double laplace_franke3d(VolumeMesh::PointT vec);

double solve_franke_poisson(VolumeMesh& mesh_,  int laplace, int face_point, int cell_point, int degree = 2, std::string meshname = "");

void solve_laplace_equation(VolumeMesh& mesh_, int laplace, int face_point, int cell_point, int degree=2);

