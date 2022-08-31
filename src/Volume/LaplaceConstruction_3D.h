
//=============================================================================

#pragma once

//=============================================================================

#include <Eigen/SparseCore>
#include <Eigen/QR>
#include <Eigen/IterativeLinearSolvers>
#include "../VolumeMesh/VolumeMesh.h"

//=============================================================================

void setup_3D_cell_face_prolongation_matrix(VolumeMesh &mesh,
                                            Eigen::SparseMatrix<double> &P,
                                            Eigen::SparseMatrix<double> &Pc,
                                            Eigen::SparseMatrix<double> &Pf);

void setup_3D_stiffness_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S,
                               int Laplace, int face_point, int cell_point);

void setup_3D_mass_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M,
                          int Laplace, int face_point, int cell_point);

//--------------------------- Dual Laplacian----------------------------------------------

void circumcenter(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
                  const Eigen::Vector3d &c, Eigen::Vector3d &cc);

void circumcenter(const Eigen::Matrix<double, 4, 3> &t, Eigen::Vector3d &c);

double volume(const Eigen::Vector3d &a, const Eigen::Vector3d &b,
              const Eigen::Vector3d &c, const Eigen::Vector3d &d);

void dualLaplace(Eigen::MatrixXd &V, Eigen::MatrixXi &T,
                 Eigen::SparseMatrix<double> &L,
                 Eigen::SparseMatrix<double> &M);

void setup_3D_dual_laplace_libigl(VolumeMesh &mesh,
                                  Eigen::SparseMatrix<double> &S);

void setup_3D_dual_mass_matrix_libigl(VolumeMesh &mesh,
                                      Eigen::SparseMatrix<double> &M);

//--------------------------- Sandwich Laplacian -------------------------------------------

void setup_3D_sandwich_stiffness_matrix(VolumeMesh &mesh,
                                        Eigen::SparseMatrix<double> &S,
                                        Eigen::SparseMatrix<double> &S_tets,
                                        int face_point, int cell_point);

void setup_3D_sandwich_mass_matrix(VolumeMesh &mesh,
                                   Eigen::SparseMatrix<double> &M,
                                   int face_point, int cell_point);



