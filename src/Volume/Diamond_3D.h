//=============================================================================
#pragma once
//=============================================================================

#include "../VolumeMesh/VolumeMesh.h"


//=============================================================================

void setup_3D_diamond_mass_matrix(VolumeMesh &mesh,
                                  Eigen::SparseMatrix<double> &M);

void setup_3D_diamond_gradient(VolumeMesh &mesh, Eigen::SparseMatrix<double> &G,
                               Eigen::SparseMatrix<double> &V);

//=== Gradient ===================================================================

Eigen::MatrixXd cc_gradient_operator(const Eigen::MatrixXd& X);

Eigen::MatrixXd cc_bdy_gradient_operator(const Eigen::MatrixXd& X);

//-------------------------------------------------------------------------------
