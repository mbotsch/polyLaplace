#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../HarmonicBase3D2D/PolyhedralMesh.hpp"
#include "../HarmonicBase3D2D/HarmonicPolyhedron.hpp"
#include "../HarmonicBase3D2D/ThreadTools.h"
#include "../HarmonicBase3D2D/WriteTetOVM.hpp"

void setup_3D_harmonic_matrices(PolyhedralMesh &mesh, Eigen::SparseMatrix<double> &S,Eigen::SparseMatrix<double> &M);

double solve_3D_Franke_harmonic(const std::string& meshname);

double solve_harmonic_eigenvalue_problem(std::string& meshname,Eigen::VectorXd &evalues, std::string& meshname_file);

void setup_3D_harmonic_stiffness_matrix(PolyhedralMesh &mesh, Eigen::SparseMatrix<double> &S);
