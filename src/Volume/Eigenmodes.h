
//=============================================================================

#pragma once

//=============================================================================

#include "../VolumeMesh/VolumeMesh.h"
#include "diffgeo_3D.h"

// https://www.dropbox.com/s/dzz7je2cbclq5gy/LB3D.pdf?dl=1 section 8.2
double solve_eigenvalue_problem(VolumeMesh &mesh_,Eigen::VectorXd &evalues,
                                int laplace, int face_point, int cell_point,
                                const std::string& meshname = "default");

void analytic_eigenvalues_unitBall(Eigen::VectorXd &eval, int n);
