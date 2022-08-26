//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include "PolyLaplace.h"

#define PMP_SCALAR_TYPE_64

//=============================================================================

void setup_subdiv_stiffness_matrix(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &S,
                                   int kind = 0);

void setup_subdiv_mass_matrix(pmp::SurfaceMesh &mesh,
                              Eigen::SparseMatrix<double> &M, int kind = 0);

void setup_subdiv_system_matrices(pmp::SurfaceMesh &mesh,
                                  Eigen::SparseMatrix<double> &S,
                                  Eigen::SparseMatrix<double> &M,
                                  Eigen::SparseMatrix<double> &P,
                                  int kind = 0);


double solve_poisson_non_dirichlet(pmp::SurfaceMesh &mesh, int function = 2,
                                   int kind = 0);

double curvature_non_dirichlet(pmp::SurfaceMesh &mesh, int kind = 0);

void setup_cc_prolongation(pmp::SurfaceMesh &mesh,
                           Eigen::SparseMatrix<double> &P,
                           int lvl);

void linear_interpolation_catmull_clark(pmp::SurfaceMesh &mesh_);

void setup_interpolating_prolongation(pmp::SurfaceMesh &mesh_,
                                      Eigen::SparseMatrix<double> &P, int lvl);

void setup_smooth_basis_matrix(pmp::SurfaceMesh &mesh, Eigen::MatrixXd &B, int kind);

void compute_reference_solution(pmp::SurfaceMesh &mesh, int kind);

extern int subdiv_lvl;