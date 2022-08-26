//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================
#include <Eigen/SparseCore>
#include <Eigen/QR>
#include <Eigen/IterativeLinearSolvers>
#include "../VolumeMesh/VolumeMesh.h"
#include "LaplaceConstruction_3D.h"
#include "../Surface/diffgeo.h"

//=============================================================================

void setup_3D_diamond_mass_matrix(VolumeMesh &mesh,
                                  Eigen::SparseMatrix<double> &M);

void setup_3D_diamond_gradient(VolumeMesh &mesh, Eigen::SparseMatrix<double> &G,
                               Eigen::SparseMatrix<double> &V);

//=== Gradient ===================================================================

Eigen::MatrixXd cc_gradient_operator(Eigen::MatrixXd X);

Eigen::MatrixXd cc_bdy_gradient_operator(Eigen::MatrixXd X);

//-------------------------------------------------------------------------------
