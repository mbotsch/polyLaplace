//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/surface_mesh.h>

//=============================================================================

//!compute the Laplacian matrix for a polygonal mesh as described in "DDiscrete Differential Operators on Polygonal Meshes"

void setup_deGoes_gradient_operator(pmp::SurfaceMesh& mesh,
                                    Eigen::SparseMatrix<double>& G);

void setup_deGoes_divergence_operator(pmp::SurfaceMesh& mesh,
                                      Eigen::SparseMatrix<double>& D);

void setup_deGoes_sharp_operator_per_face(pmp::SurfaceMesh& mesh, pmp::Face f,
                                          Eigen::MatrixXd& U);

void setup_deGoes_discrete_flat_operator_per_face(pmp::SurfaceMesh& mesh,
                                                  pmp::Face f,
                                                  Eigen::MatrixXd& V);

void setup_deGoes_projection_operator_per_face(pmp::SurfaceMesh& mesh,
                                               pmp::Face f, Eigen::MatrixXd& P);

void setup_deGoes_inner_product_per_face(pmp::SurfaceMesh& mesh, pmp::Face f,
                                         Eigen::MatrixXd& M);

void setup_deGoes_laplace_operator(pmp::SurfaceMesh& mesh,
                                   Eigen::SparseMatrix<double>& L);

void setup_deGoes_mass_matrix(pmp::SurfaceMesh& mesh,
                              Eigen::SparseMatrix<double>& M);

void setup_deGoes_face_area_matrix(pmp::SurfaceMesh& mesh,
                                   Eigen::SparseMatrix<double>& A);

extern float deGoes_laplace_lambda_;
