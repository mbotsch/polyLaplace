#include <iostream>
#include <fstream>
#include <vector>
#include <functional>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/SparseCore>
#include <pmp/SurfaceMesh.h>

void buildStiffnessAndMass2d(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double>& K, Eigen::SparseMatrix<double>& M);

void buildStiffness2d(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double>& K);