//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include "PolyLaplace.h"

//=============================================================================

void setup_sec_Laplace_matrix(pmp::SurfaceMesh &mesh,
                              Eigen::SparseMatrix<double> &L);

void setup_sec_mass_matrix(pmp::SurfaceMesh &mesh,
                           Eigen::SparseMatrix<double> &M);

void setup_sec_gradient_operator(pmp::SurfaceMesh &mesh,
                                 Eigen::SparseMatrix<double> &G);

void setup_cc_P_matrix(pmp::SurfaceMesh &mesh_,
                         Eigen::SparseMatrix<double> &A0, int lvl,std::vector<Eigen::SparseMatrix<double>> &Al);

void setup_cc_retain_P_matrix(pmp::SurfaceMesh &mesh_,
                       Eigen::SparseMatrix<double> &A0, int lvl);

void setup_mod_butterfly_P_matrix(pmp::SurfaceMesh &mesh_,
                              Eigen::SparseMatrix<double> &A0, int lvl);

void setup_bdry_cc_A0_matrix(pmp::SurfaceMesh &mesh,
                         Eigen::SparseMatrix<double> &A0, int lvl, std::vector<Eigen::SparseMatrix<double>> &Al);

void setup_interpolating_sec_A0_matrix(pmp::SurfaceMesh &mesh,
                         Eigen::SparseMatrix<double> &A0, int lvl, std::vector<Eigen::SparseMatrix<double>> &Al);

void linear_interpolation_catmull_clark(pmp::SurfaceMesh &mesh_);

void setup_sec_divergence_operator(pmp::SurfaceMesh &mesh,
                                 Eigen::SparseMatrix<double> &D);

void normalize_sec_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
                             const Eigen::VectorXd &h);

void setup_smooth_laplace_matrices(pmp::SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M);


extern float sec_laplace_lambda_;

extern int sec_laplace_lvl;