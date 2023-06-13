//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#pragma once

//=============================================================================

#include "../VolumeMesh/VolumeMesh.h"
#include "Diamond_3D.h"
#include "PolySimpleLaplace.h"

//=============================================================================

void setup_3D_cell_face_prolongation_matrix(VolumeMesh& mesh,
                                            Eigen::SparseMatrix<double>& P,
                                            Eigen::SparseMatrix<double>& Pc,
                                            Eigen::SparseMatrix<double>& Pf);

void setup_3D_stiffness_matrix(VolumeMesh& mesh, Eigen::SparseMatrix<double>& S,
                               int Laplace, int face_point, int cell_point);

void setup_3D_mass_matrix(VolumeMesh& mesh, Eigen::SparseMatrix<double>& M,
                          int Laplace, int face_point, int cell_point);

//-------------------------------------------------------------------------

double volume(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
              const Eigen::Vector3d& c, const Eigen::Vector3d& d);
