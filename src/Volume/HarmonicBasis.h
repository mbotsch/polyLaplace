#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include "../Harmonic/PolyhedralMesh.h"
#include "../Harmonic/HarmonicPolyhedron.h"
#include "../Harmonic/ThreadTools.h"
#include "../Harmonic/WriteTetOVM.h"

void setup_3D_harmonic_matrices(PolyhedralMesh& mesh,
                                Eigen::SparseMatrix<double>& S,
                                Eigen::SparseMatrix<double>& M);

double solve_3D_Franke_harmonic(const std::string& meshname);

double solve_harmonic_eigenvalue_problem(std::string& meshname,
                                         Eigen::VectorXd& evalues,
                                         std::string& meshname_file);

void setup_3D_harmonic_stiffness_matrix(PolyhedralMesh& mesh,
                                        Eigen::SparseMatrix<double>& S);
