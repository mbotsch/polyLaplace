//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

#define PMP_SCALAR_TYPE_64

//=============================================================================

//!compute the Laplacian matrix for a polygonal mesh as described in "DDiscrete Differential Operators on Polygonal Meshes"

void setup_disney_gradient_operator(pmp::SurfaceMesh &mesh,  Eigen::SparseMatrix<double>  &G);

void setup_disney_divergence_operator(pmp::SurfaceMesh &mesh,  Eigen::SparseMatrix<double>  &D);

void setup_disney_sharp_operator_per_face(pmp::SurfaceMesh &mesh, pmp::Face f,Eigen::MatrixXd &U);

void setup_disney_discrete_flat_operator_per_face(pmp::SurfaceMesh &mesh,pmp::Face f,Eigen::MatrixXd &V);

void setup_disney_projection_operator_per_face(pmp::SurfaceMesh &mesh,pmp::Face f, Eigen::MatrixXd &P);

void setup_disney_inner_product_per_face(pmp::SurfaceMesh &mesh,pmp::Face f, Eigen::MatrixXd &M);

void setup_disney_laplace_operator(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &L);

void setup_disney_mass_matrix(pmp::SurfaceMesh &mesh,
                              Eigen::SparseMatrix<double> &M);

void setup_disney_face_area_matrix(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &A);

void normalize_disney_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
const Eigen::VectorXd &h);

extern float disney_laplace_lambda_;
