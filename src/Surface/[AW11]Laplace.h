//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

//=============================================================================

//!compute the Laplacian matrix for a polygonal mesh as described in "Discrete Laplacians on general polygonal meshes" [Alexa & Wardetzky, 2011]

void setup_E_and_B_perFace(pmp::SurfaceMesh &mesh, pmp::Face f,
                           Eigen::MatrixXd &E, Eigen::MatrixXd &B);

void setup_poly_Laplace_matrix(pmp::SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &L);

void setup_poly_mass_matrix(pmp::SurfaceMesh &mesh,
                            Eigen::SparseMatrix<double> &M);

void setup_poly_gradient_operator(pmp::SurfaceMesh &mesh,
                                  Eigen::SparseMatrix<double> &G);

void normalize_poly_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
                              const Eigen::VectorXd &h);

void setup_poly_divergence_operator(pmp::SurfaceMesh &mesh,
                                    Eigen::SparseMatrix<double> &D);

extern float poly_laplace_lambda_;
