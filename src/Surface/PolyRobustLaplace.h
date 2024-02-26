#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "diffgeo.h"
#include "LaplaceConstruction.h"

#ifndef GLOBALS_H
#define GLOBALS_H
const int penalty_ = 100000;
#endif

void find_trace_minimizer_weights(const Eigen::MatrixXd& polyVerts, Eigen::VectorXd &weights, bool use_fallback = true);

bool positive_diagonally_dominant(Eigen::Matrix2d& H, double _eps);

void project_positive_definite(Eigen::Matrix2d& H, double _eigenvalue_eps = 1e-9);

double getEnergy(const Eigen::Vector2d& vVert, const Eigen::MatrixXd& poly);
double cotan(Eigen::Vector2d& vL, Eigen::Vector2d& vC, Eigen::Vector2d& vR);
void computeHarmonicWeights(const Eigen::MatrixXd &poly, Eigen::VectorXd& weights, Eigen::Vector2d vV);
void getLeastSquaresPoint(const Eigen::MatrixXd& poly, Eigen::VectorXd& x);
void computeMaxAreaPolygonProjection(const Eigen::MatrixXd& poly, Eigen::Vector3d &n, Eigen::MatrixXd& pPoly, Eigen::MatrixXd& projMatrix);
double getStiffnessTrace(const Eigen::Vector3d& vVert,
                                          const Eigen::MatrixXd& poly,
                                          const Eigen::VectorXd& weights);
std::tuple<double, Eigen::Vector2d, Eigen::Matrix2d> getEnergyAndGradientAndHessian(const Eigen::Vector2d& vVert, const Eigen::MatrixXd& poly);

bool armijo_condition(
    double _f_curr,
    double _f_new,
    double _s,
    const Eigen::Vector2d& _d,
    const Eigen::Vector2d& _g,
    double _armijo_const = 1e-4);

Eigen::Vector2d line_search(
    const Eigen::Vector2d& _x0, const Eigen::Vector2d& _d, double _f,
    const Eigen::Vector2d& _g, const Eigen::MatrixXd& poly,
    double _s_max = 1.0, // Initial step size
    double _shrink = 0.8,
    int _max_iters = 64,
    double _armijo_const = 1e-4);