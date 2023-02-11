//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/Triangulation.h>


//=============================================================================


//!compute the transformation matrix for an arbitrary mesh that inserts a chosen point per face .
void setup_prolongation_matrix(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &P);

//! Computes a diagonal matrix where the entries are the sums of the rows.
void lump_matrix(Eigen::SparseMatrix<double> &D);

//! Computes the squared triangle area minimizing points and its convex combination weights
//! for each face and stores it in a prior defined property.
void setup_face_point_properties(pmp::SurfaceMesh &mesh, unsigned int min_point);

//------------------------Point and Weight minimizer -----------------------------------------------------------------

void find_area_minimizer_weights(const Eigen::MatrixXd &poly,
                       Eigen::VectorXd &weights);

//------------------------ Miscellaneous-----------------------------------------------------------------

void fit_plane_to_polygon(const Eigen::MatrixXd &poly, Eigen::Vector3d &normal, Eigen::Vector3d &origin);


//=============================================================================