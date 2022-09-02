#pragma once
//=============================================================================
#include "../VolumeMesh/VolumeMesh.h"

//=============================================================================

void setup_3D_polysimple_stiffness_matrix(VolumeMesh &mesh,
                                          Eigen::SparseMatrix<double> &S,
                                          Eigen::SparseMatrix<double> &S_tets,
                                          int face_point, int cell_point);

void setup_3D_polysimple_mass_matrix(VolumeMesh &mesh,
                                     Eigen::SparseMatrix<double> &M,
                                     int face_point, int cell_point);

