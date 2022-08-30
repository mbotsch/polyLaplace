#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <pmp/SurfaceMesh.h>

void buildStiffnessAndMass2d(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double>& S, Eigen::SparseMatrix<double>& M);

double solve_2D_Franke_harmonic(pmp::SurfaceMesh &mesh);

double solve_2D_harmonic_eigenvalue_problem(pmp::SurfaceMesh &mesh, Eigen::VectorXd &evalues, std::string& meshname_file);

