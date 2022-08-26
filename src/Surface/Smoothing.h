//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "LaplaceConstruction.h"

//=============================================================================

using namespace pmp;

//=============================================================================

class Smoothing
{
public:
    Smoothing(SurfaceMesh& mesh)
            : mesh_(mesh), vertices_(0), faces_(0), min_point_(0), laplace_(0)
    {
    }

    //! Perform implicit Laplacian smoothing with \c timestep.
//    void implicit_smoothing(Scalar timestep, unsigned int laplace,
//                            unsigned int min_point);

    void implicit_smoothing_misha(Scalar timestep, unsigned int laplace,
                                  unsigned int min_point);

private:
    void update_Laplaces(unsigned int laplace,unsigned int min_point)
    {
        if (mesh_.n_faces() != faces_ || mesh_.n_vertices() != vertices_
            || min_point_ != min_point || laplace_ != laplace)
        {
            vertices_ = mesh_.n_vertices();
            faces_ = mesh_.n_faces();
            min_point_ = min_point;
            laplace_ = laplace;

//            L_.resize(vertices_ + faces_, vertices_ + faces_);
//            L_triangle.resize(vertices_, vertices_);
            std::cout << "Stiffness matrix has been updated" << std::endl;
            setup_stiffness_matrices(mesh_,S_, laplace_, min_point_);

            // setup the new stiffness matrix
            Eigen::SparseMatrix<double> L_ , L;
//
//            setup_Laplace_matrix_local_Philipp(mesh_,L_,min_point);
//            setup_stiffness_matrices(mesh_,L_, 0, min_point_);
//            setup_stiffness_matrices(mesh_,L, 1, min_point_);
//            setup_poly_Laplace_matrix(mesh_, L);
//            std::cout << (L_-L).norm() << std::endl;
//            setup_triangle_Laplace_matrix(mesh_, L_triangle);
        }
    }

private:
    SurfaceMesh& mesh_;

//    Eigen::SparseMatrix<double> L_, L_poly, L_triangle;
    Eigen::SparseMatrix<double> S_;
    unsigned int vertices_, faces_, min_point_, laplace_;
};

//=============================================================================
