//=============================================================================
#pragma once
//=============================================================================
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include <pmp/algorithms/DifferentialGeometry.h>
#include <pmp/algorithms/SurfaceTriangulation.h>

#define PMP_SCALAR_TYPE_64

enum CoarseDimension{
        Vertices,
        Edges,
        Virtual_Edges,
        Refined_mesh
};

enum MG_Solver
{
    RELAXER_JACOBI ,
    RELAXER_GAUSS_SEIDEL ,
    RELAXER_PARALLEL_GAUSS_SEIDEL ,
    RELAXER_COUNT
};//=============================================================================

//!compute the Laplacian matrix for a triangle mesh.
void setup_triangle_Laplace_matrix(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &L);

//!compute the transformation matrix for an arbitrary mesh that inserts a chosen point per face .
void setup_prolongation_matrix(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &A);

//! compute the mass matrix for a triangle mesh.
void setup_triangle_mass_matrix(pmp::SurfaceMesh &mesh,
                                Eigen::SparseMatrix<double> &M);

//! Computes a diagonal matrix where the entries are the sums of the rows.
void lump_matrix(Eigen::SparseMatrix<double> &D);


//! Computes a diagonal matrix where the entries are the sums of the absolute row entries.
void abs_lump_matrix(Eigen::SparseMatrix<double> &D);

//! Inverts a given diagonal matrix.
void invert_mass_matrix(Eigen::SparseMatrix<double> &M);

//! count number of negative entries in a sparse matrix
unsigned int num_negative(const Eigen::SparseMatrix<double> &_mat);

//! count number of nonzero entries in a sparse matrix without diagonal elements
unsigned int num_nonzero(const Eigen::SparseMatrix<double> &_mat);

//! count number of nonzero entries in a sparse matrix
unsigned int sparsity(const Eigen::SparseMatrix<double> &_mat);

//----------------------------------area computations------------------------------------------------------

//! barycenter/centroid of mesh, computed as area-weighted mean of vertices.
pmp::Point my_centroid(const pmp::SurfaceMesh &mesh);

//! computes the area of a face.
double face_area(const pmp::SurfaceMesh &mesh, pmp::Face f);

//! surface area of the mesh.
pmp::Scalar my_surface_area(const pmp::SurfaceMesh &mesh);

//---------------------Gradient/ Divergence --------------------------------------------------------------

//! Computes the gradient on a triangle formed by the points i, j and k.
Eigen::Vector3d gradient_hat_function(pmp::Point i, pmp::Point j, pmp::Point k);

Eigen::Vector3d gradient_hat_function(Eigen::Vector3d i, Eigen::Vector3d j, Eigen::Vector3d k);

//! Computes the Gradient matrix for any given mesh.
void setup_Gradient_Matrix(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G);

void setup_triangle_Gradient_Matrix(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G);

//! Computes the Divergence matrix for any given mesh.
void setup_Divergence_Matrix(pmp::SurfaceMesh &mesh,
                             Eigen::SparseMatrix<double> &Gt);

//! Computes the diagonal mass matrix W needed for L = DWG.
void setup_Gradient_Mass_Matrix(pmp::SurfaceMesh &mesh,
                                Eigen::SparseMatrix<double> &M);

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

void insert_points(pmp::SurfaceMesh& mesh_,unsigned int minpoint);

//=============================================================================