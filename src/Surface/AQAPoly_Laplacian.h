//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"

#define PMP_SCALAR_TYPE_64


enum MatrixType
{
    Mass,
    Stiffness
};

//=============================================================================

void setup_AQAPoly_matrices(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M, MatrixType matrix, CoarseDimension dof,
                               int minpoint, int degree=2);

template<unsigned int Degree>
void setup_local_AQAPoly_matrices(pmp::SurfaceMesh &mesh,std::vector<Eigen::Triplet<double>> &triplet, MatrixType matrix, CoarseDimension dof,
                                  int minpoint );

//---------------------------------------------------------------------------

void local_matrix_to_global(pmp::SurfaceMesh &mesh, pmp::Face &f, Eigen::MatrixXd &M,
                            std::vector<Eigen::Triplet<double>> &triplet, int degree);

void local_matrix_to_global(pmp::SurfaceMesh &mesh, pmp::Face &f, Eigen::MatrixXd &M,
                            std::vector<Eigen::Triplet<double>> &triplet,CoarseDimension dof);


int edge_node_to_idx(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he);

//int vertex_node_to_idx(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Vertex &v);
int vertex_node_to_idx(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he);

//========================================================================
void setup_triangle_FEM_stiffness_matrix(pmp::SurfaceMesh &mesh,
                                         Eigen::SparseMatrix<double> &S,
                                         bool linear);

void local_triangle_fem_matrix(pmp::SurfaceMesh &mesh, pmp::Face f,
                               std::vector<Eigen::Triplet<double>> &trip,
                               bool linear);

void setup_triangle_FEM_mass_matrix(pmp::SurfaceMesh &mesh,
                                    Eigen::SparseMatrix<double> &M);

void local_triangle_mass_fem_matrix(pmp::SurfaceMesh &mesh, pmp::Face f,
                                    std::vector<Eigen::Triplet<double>> &trip);

void project_triangle(const pmp::Point &p0, const pmp::Point &p1, const pmp::Point &p2,
                      pmp::vec2 &z0, pmp::vec2 &z1, pmp::vec2 &z2);

void problem_faces(pmp::SurfaceMesh &mesh, int minpoint);

template<unsigned int Degree>
double solve_2D_AQAPoly_Poisson(pmp::SurfaceMesh &mesh, CoarseDimension dof);

double solve_2D_AQAPoly_Poisson( pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree);

double solve_AQAPoly_Poisson_mg(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree, bool direct = true, MG_Solver solver = RELAXER_GAUSS_SEIDEL, int vcycles = 10, int iterations =10);

template<unsigned int Degree>
double solve_AQAPoly_Poisson_mg(pmp::SurfaceMesh &mesh, CoarseDimension dof, bool direct, MG_Solver solver, int vcycles, int iterations);

double AQAPoly_condition_nr(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree, std::vector<double> &condnr);

template<unsigned int Degree>
double AQAPoly_condition_nr(pmp::SurfaceMesh &mesh, CoarseDimension dof, std::vector<double> &condnr);

void project_triangle(Eigen::Matrix3d &Triangle3, Eigen::Matrix2d &Triangle2);
double TriangleJacobian(Eigen::Matrix3d Triangle, Eigen::Matrix2d &Jacobian);

// Hacky, dont use in general
double naive_Franke(pmp::SurfaceMesh &mesh,std::string &meshname);

template<unsigned int Degree>
double test_quadratic_reproduction(pmp::SurfaceMesh &mesh,CoarseDimension dof);

double test_quadratic_reproduction(pmp::SurfaceMesh &mesh, CoarseDimension dof, int degree);