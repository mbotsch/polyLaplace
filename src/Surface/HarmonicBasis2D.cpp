//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch, Philipp Herholz.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#include "HarmonicBasis2D.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "../Harmonic/HarmonicPolygon.h"

void buildStiffnessAndMass2d(pmp::SurfaceMesh& mesh,
                             Eigen::SparseMatrix<double>& K,
                             Eigen::SparseMatrix<double>& M, int nKernel,
                             int nProbes)
{
    Eigen::MatrixXd V(mesh.n_vertices(), 3);
    std::vector<std::vector<int>> poly;

    for (auto v : mesh.vertices())
    {
        Eigen::Vector3d p(mesh.position(v)[0], mesh.position(v)[1],
                          mesh.position(v)[2]);
        V.row(v.idx()) = p;
    }
    for (auto f : mesh.faces())
    {
        std::vector<int> face;
        for (auto fv : mesh.vertices(f))
        {
            face.emplace_back(fv.idx());
        }
        poly.emplace_back(face);
    }

    const int nv = (int)V.rows();
    std::vector<Eigen::Triplet<double>> tripK, tripM;

    for (auto& p : poly)
    {
        Eigen::MatrixXd pts(p.size(), 3);
        for (int i = 0; i < (int)p.size(); ++i)
        {
            pts(i, 0) = V(p[i], 0);
            pts(i, 1) = V(p[i], 1);
            pts(i, 2) = V(p[i], 2);
        }

        HarmonicPolygon hp(pts, nKernel, nProbes);

        Eigen::MatrixXd Ki, Mi;

        hp.stiffnessMatrix(Ki);
        hp.massMatrix(Mi);

        for (int i = 0; i < (int)p.size(); ++i)
        {
            for (int j = 0; j < (int)p.size(); ++j)
            {
                tripK.emplace_back(p[i], p[j], Ki(i, j));
                tripM.emplace_back(p[i], p[j], Mi(i, j));
            }
        }
    }

    K.resize(nv, nv);
    K.setFromTriplets(tripK.begin(), tripK.end());

    M.resize(nv, nv);
    M.setFromTriplets(tripM.begin(), tripM.end());
}

void buildStiffness2d(pmp::SurfaceMesh& mesh, Eigen::SparseMatrix<double>& K,
                      int nKernel, int nProbes)
{
    Eigen::MatrixXd V(mesh.n_vertices(), 3);
    std::vector<std::vector<int>> poly;

    for (auto v : mesh.vertices())
    {
        Eigen::Vector3d p(mesh.position(v)[0], mesh.position(v)[1],
                          mesh.position(v)[2]);
        V.row(v.idx()) = p;
    }
    for (auto f : mesh.faces())
    {
        std::vector<int> face;
        for (auto fv : mesh.vertices(f))
        {
            face.emplace_back(fv.idx());
        }
        poly.emplace_back(face);
    }

    const int nv = (int)V.rows();
    std::vector<Eigen::Triplet<double>> tripK, tripM;

    for (auto& p : poly)
    {
        Eigen::MatrixXd pts(p.size(), 3);
        for (int i = 0; i < (int)p.size(); ++i)
        {
            pts(i, 0) = V(p[i], 0);
            pts(i, 1) = V(p[i], 1);
            pts(i, 2) = V(p[i], 2);
        }
        HarmonicPolygon hp(pts, nKernel, nProbes);

        Eigen::MatrixXd Ki, Mi;
        hp.stiffnessMatrix(Ki);

        for (int i = 0; i < (int)p.size(); ++i)
        {
            for (int j = 0; j < (int)p.size(); ++j)
            {
                tripK.emplace_back(p[i], p[j], Ki(i, j));
            }
        }
    }
    K.resize(nv, nv);
    K.setFromTriplets(tripK.begin(), tripK.end());
}
