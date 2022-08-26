//=============================================================================
#pragma once
//=============================================================================

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "../VolumeMesh/VolumeMesh.h"
#include "diffgeo_3D.h"
#include "../Surface/diffgeo.h"
//=============================================================================

template<unsigned int Degree>
double solve_3D_AQAPoly_Poisson(std::string &meshname, CoarseDimension dof, int function );

double solve_3D_AQAPoly_Poisson( std::string &meshname, CoarseDimension dof, int degree, int function=0);

template<unsigned int Degree>
double solve_3D_AQAPoly_EigenModes(std::string &meshname, CoarseDimension dof,std::string &tesselation);

double solve_3D_AQAPoly_EigenModes( std::string &meshname, CoarseDimension dof, int degree,std::string &tesselation);
//========================================================================

double solve_3D_AQAPoly_Poisson_mg(std::string &meshname, std::ofstream &timings_file,CoarseDimension dof, int degree, bool direct = true, MG_Solver solver = RELAXER_PARALLEL_GAUSS_SEIDEL, int vcycles = 10, int iterations =10);

template<unsigned int Degree>
double solve_3D_AQAPoly_Poisson_mg(std::string &meshname,std::ofstream &timings_file, CoarseDimension dof, bool direct, MG_Solver solver, int vcycles, int iterations);

template<unsigned int Degree>
double compare_3D_AQAPoly_Poisson_mg(std::string &meshname,std::ofstream &timings_file, CoarseDimension dof, bool direct, MG_Solver solver, int vcycles, int iterations,bool write_direct);

double solve_3D_AQAPoly_Poisson_mg_extended(std::string &meshname, std::ofstream &timings_file,CoarseDimension dof, int degree, bool direct, MG_Solver solver = RELAXER_PARALLEL_GAUSS_SEIDEL, int vcycles = 10, int iterations =10, bool write_direct = true);
