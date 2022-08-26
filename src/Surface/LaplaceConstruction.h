//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceTriangulation.h>
#include <unsupported/Eigen/NonLinearOptimization>
#include "AQAPoly_Laplacian.h"

//=============================================================================

using namespace pmp;
#define PMP_SCALAR_TYPE_64

//=============================================================================

void setup_sandwich_stiffness_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S, int minpoint = 0);

void setup_sandwich_mass_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, int minpoint = 0);

void setup_sandwich_gradient_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G, int minpoint = 0);

void setup_sandwich_divergence_matrix(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &D, int minpoint = 0);

void localCotanMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                      Eigen::MatrixXd &L);

void localMassMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                     Eigen::MatrixXd &M);

void testOperator(SurfaceMesh &mesh, unsigned int laplace, unsigned int min_point);

void localGradientMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min, Eigen::VectorXd &w,
                      Eigen::MatrixXd &G);

//----------------------------------------------------------------------------------

void setup_stiffness_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S, int Laplace, int minpoint = 0, int degree = 2, CoarseDimension coarseningType =Edges);

void setup_mass_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, int Laplace, int minpoint = 0,int degree = 2,CoarseDimension coarseningType=Edges,
                         bool lumped = true);
void setup_simple_2D_prolongation(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &P);