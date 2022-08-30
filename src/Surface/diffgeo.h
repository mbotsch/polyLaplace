//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceTriangulation.h>


//=============================================================================


//!compute the transformation matrix for an arbitrary mesh that inserts a chosen point per face .
void setup_prolongation_matrix(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &P);



//! Computes a diagonal matrix where the entries are the sums of the rows.
void lump_matrix(Eigen::SparseMatrix<double> &D);

//---------------------Gradient/ Divergence --------------------------------------------------------------

//! Computes the gradient on a triangle formed by the points i, j and k.
Eigen::Vector3d gradient_hat_function(pmp::Point i, pmp::Point j, pmp::Point k);

Eigen::Vector3d gradient_hat_function(Eigen::Vector3d i, Eigen::Vector3d j, Eigen::Vector3d k);

//! Computes the squared triangle area minimizing points and its convex combination weights
//! for each face and stores it in a prior defined property.
void setup_face_point_properties(pmp::SurfaceMesh &mesh, unsigned int min_point);


//------------------------Point and Weight minimizer -----------------------------------------------------------------

void find_area_minimizer_weights(const Eigen::MatrixXd &poly,
                       Eigen::VectorXd &weights);

void find_weights_for_point(const Eigen::MatrixXd &poly,
                            const Eigen::Vector3d &point,
                            Eigen::VectorXd &weights);

void absAreaJacobian(const Eigen::Matrix<double, -1, 3> &poly, const Eigen::Vector3d &a, Eigen::Vector3d &jac);

void absAreaHessian(const Eigen::Matrix<double, -1, 3> &poly, const Eigen::Vector3d &a, Eigen::Matrix3d &hess);

double absTriangleArea(const Eigen::MatrixXd &poly, const Eigen::VectorXd &a);

void optimizeAbsoluteTriangleArea(const Eigen::MatrixXd &poly, Eigen::Vector3d &x,
                                  const int maxIter = 20, const int maxLineSearchIter = 20);

//! Computes the circumcenter of a triangle formed by the points i, j and k.
Eigen::Vector3d triangle_circumcenter(const Eigen::MatrixXd &poly);

void barycentric_weights_triangle(const Eigen::MatrixXd &poly,
                                  const Eigen::Vector3d &point,
                                  Eigen::VectorXd &weights);
//=============================================================================