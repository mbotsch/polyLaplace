//=============================================================================
#pragma once
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "diffgeo.h"
#include "LaplaceConstruction.h"

//=============================================================================

using namespace pmp;

//=============================================================================

class Smoothing
{
public:
    Smoothing(SurfaceMesh& mesh)
        : mesh_(mesh), vertices_(0), faces_(0), min_point_(0), laplace_(-1)
    {
    }

    void implicit_smoothing(Scalar timestep, unsigned int laplace,
                            unsigned int min_point, bool rescale);

private:
    void update_Laplaces(unsigned int laplace, unsigned int min_point)
    {
        if (mesh_.n_faces() != faces_ || mesh_.n_vertices() != vertices_ ||
            min_point_ != min_point || laplace_ != laplace)
        {
            vertices_ = mesh_.n_vertices();
            faces_ = mesh_.n_faces();
            min_point_ = min_point;
            laplace_ = laplace;

            std::cout << "Stiffness matrix has been updated" << std::endl;
            setup_stiffness_matrices(mesh_, S_, laplace_, min_point_);
        }
    }

private:
    SurfaceMesh& mesh_;
    Eigen::SparseMatrix<double> S_;
    unsigned int vertices_, faces_, min_point_, laplace_;
};

//=============================================================================
