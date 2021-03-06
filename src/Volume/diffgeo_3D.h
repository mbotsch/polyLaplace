//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
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
#include <VolumeMesh.h>
#include "../Surface/diffgeo.h"
//=============================================================================

using namespace pmp;

//=============================================================================

void find_volume_minimizing_point(
    const std::vector<Eigen::Matrix3d,
                      Eigen::aligned_allocator<Eigen::Matrix3d>> &tetrahedrons,
    Eigen::Vector3d &point);

void find_weights_for_point_3d(
    const std::vector<Eigen::Vector3d,
                      Eigen::aligned_allocator<Eigen::Vector3d>> &poly,
    const Eigen::Vector3d &point, Eigen::VectorXd &weights);

// Computes virtual cell and face minimizer and their weights. Information is stored in properties
void compute_3D_virtual_points_and_weights(VolumeMesh &mesh, int face_point = 0,
                                           int cell_point = 0);

void compute_3D_virtual_points_and_weights_for_cells(VolumeMesh &mesh,
                                                     int cell_point = 0);

//Compute the  signed volume of the tetrahedron formed by the 4 given points
double volume_tetrahedron_signed(VolumeMesh::PointT i, VolumeMesh::PointT j,
                                 VolumeMesh::PointT k, VolumeMesh::PointT l);

//Divides Polyhedral Mesh into tetrahedron
VolumeMesh::PointT compute_triangle_normal(VolumeMesh::PointT a,
                                           VolumeMesh::PointT b,
                                           VolumeMesh::PointT c);

double compute_triangle_area(VolumeMesh::PointT p0, VolumeMesh::PointT p1,
                             VolumeMesh::PointT p2);