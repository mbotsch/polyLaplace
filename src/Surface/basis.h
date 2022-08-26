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

//=============================================================================

void linear_basis_gradient(Eigen::MatrixXd &basis);

void grad_phi_1(Eigen::Vector2d &val, double x, double y);
void grad_phi_2(Eigen::Vector2d &val, double x, double y);
void grad_phi_3(Eigen::Vector2d &val, double x, double y);
void grad_phi_4(Eigen::Vector2d &val, double x, double y);
void grad_phi_5(Eigen::Vector2d &val, double x, double y);
void grad_phi_6(Eigen::Vector2d &val, double x, double y);

double phi_1_lin(double x, double y);
double phi_2_lin(double x, double y);
double phi_3_lin(double x, double y);

double phi_1_quad(double x, double y);
double phi_2_quad(double x, double y);
double phi_3_quad(double x, double y);
double phi_4_quad(double x, double y);
double phi_5_quad(double x, double y);
double phi_6_quad(double x, double y);


int local_to_global_index(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he, int i);

double gauss_quadrature_gradient(Eigen::MatrixXd &jacobian, int bi, int bj);

double gauss_quadrature_mass(int bi, int bj);


void compute_basis_gradient(Eigen::Vector2d &grad,int base_idx, double x, double y);

double compute_triangle_basis(int basis_idx, double x, double y, bool lin);

double compute_quad_basis(int basis_idx, double x, double y, bool lin);


