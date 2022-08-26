//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>

//=============================================================================
using namespace pmp;

void setup_stiffness_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S, int minpoint = 0);

void setup_mass_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, int minpoint = 0);

void setup_sandwich_gradient_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G, int minpoint = 0);

void setup_sandwich_divergence_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &D, int minpoint = 0);

void localCotanMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                      Eigen::MatrixXd &L);

void localMassMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                     Eigen::MatrixXd &M);

void localGradientMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                         Eigen::MatrixXd &G);

//----------------------------------------------------------------------------------
