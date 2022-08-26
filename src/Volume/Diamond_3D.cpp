//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "Diamond_3D.h"
#include <pmp/MatVec.h>
#include "diffgeo_3D.h"

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

void setup_3D_diamond_mass_matrix(VolumeMesh &mesh,
                                  Eigen::SparseMatrix<double> &M)
{

    std::vector<T> triplets;
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");

    triplets.reserve(mesh.n_cells() * 4 * 4);
    int n_v = mesh.n_vertices();
    int n_f = mesh.n_faces();
    int n_c = mesh.n_cells();

    VolumeMesh::PointT A, B, E_AB, E_AF, E_BF, F, K, L, E, D;
    Eigen::Vector3d N_kl, N_ab, N_ef;
    int l_idx, k_idx, a_idx, b_idx, f_idx;
    int f_count = 0;
    for (auto f : mesh.faces())
    {
        f_count++;
        // face midpoint
        F = f_prop[f];
        f_idx = n_v + f.idx();
        //cells adjacent to the face
        if (!mesh.is_boundary(f))
        {
            L = c_prop[mesh.face_cells(f)[1]];
            K = c_prop[mesh.face_cells(f)[0]];

            l_idx = n_v + n_f + mesh.face_cells(f)[1].idx();
            k_idx = n_v + n_f + mesh.face_cells(f)[0].idx();
        }
        else
        {
            for (auto hf : mesh.face_halffaces(f))
            {
                if (mesh.is_boundary(hf))
                {
                    L = F;
                    l_idx = f_idx;
                }
                else
                {
                    K = c_prop[mesh.incident_cell(hf)];
                    k_idx = n_v + n_f + mesh.incident_cell(hf).idx();
                }
            }
        }
        // iterate over all face edges to get diamonds
        auto hf = mesh.face_halffaces(f)[0];
        if (mesh.is_boundary(hf))
            hf = mesh.face_halffaces(f)[1];
        for (auto halfedge : mesh.halfface_halfedges(hf))
        {

            A = mesh.vertex(mesh.halfedge_vertices(halfedge)[1]);
            B = mesh.vertex(mesh.halfedge_vertices(halfedge)[0]);

            a_idx = mesh.halfedge_vertices(halfedge)[1].idx();
            b_idx = mesh.halfedge_vertices(halfedge)[0].idx();

            double volume_K = volume_tetrahedron_signed(B, A, F, K);
            double volume_L = volume_tetrahedron_signed(A, B, F, L);
            double vol = (volume_K + volume_L);

            triplets.emplace_back(l_idx, l_idx, vol / 2.0);
            triplets.emplace_back(k_idx, k_idx, vol / 2.0);

            triplets.emplace_back(a_idx, a_idx, vol / 3.0);
            triplets.emplace_back(b_idx, b_idx, vol / 3.0);
            triplets.emplace_back(f_idx, f_idx, vol / 3.0);
        }
    }
    M.resize(n_v + n_f + n_c, n_v + n_f + n_c);
    M.setFromTriplets(triplets.begin(), triplets.end());
}

//===========================Face and Volume Point Fitting=============================================================

void setup_3D_diamond_gradient(VolumeMesh &mesh, Eigen::SparseMatrix<double> &G,
                               Eigen::SparseMatrix<double> &V)
{

    std::vector<T> triplets, triplets_volume;

    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");

    triplets.reserve(mesh.n_cells() * 4 * 4);
    triplets_volume.reserve(mesh.n_cells() * 4 * 4);
    int n_v = mesh.n_vertices();
    int n_f = mesh.n_faces();
    int n_c = mesh.n_cells();
    int diamond_valence;
    VolumeMesh::PointT A, B, F, C0, C1;
    int c1_idx, c0_idx, a_idx, b_idx, f_idx;
    int row = 0;
    int nc;
    for (auto f : mesh.faces())
    {
        // face midpoint
        F = f_prop[f];
        f_idx = n_v + f.idx();
        //cells adjacent to the face
        if (!mesh.is_boundary(f))
        {
            C1 = c_prop[mesh.face_cells(f)[1]];
            C0 = c_prop[mesh.face_cells(f)[0]];

            c1_idx = n_v + n_f + mesh.face_cells(f)[1].idx();
            c0_idx = n_v + n_f + mesh.face_cells(f)[0].idx();

            //three face and two cell vertices
            diamond_valence = 5;
            //nr cell vertices
            nc = 2;
        }
        else
        {
            for (auto hf : mesh.face_halffaces(f))
            {
                if (mesh.is_boundary(hf))
                {
                    C1 = F;
                    c1_idx = f_idx;
                }
                else
                {
                    C0 = c_prop[mesh.incident_cell(hf)];
                    c0_idx = n_v + n_f + mesh.incident_cell(hf).idx();
                }
            }
            //three face and one cell vertex
            diamond_valence = 4;
            nc = 1;
        }
        // iterate over all face edges to get diamonds
        auto hf = mesh.face_halffaces(f)[0];
        if (mesh.is_boundary(hf))
            hf = mesh.face_halffaces(f)[1];
        for (auto halfedge : mesh.halfface_halfedges(hf))
        {

            A = mesh.vertex(mesh.halfedge_vertices(halfedge)[1]);
            B = mesh.vertex(mesh.halfedge_vertices(halfedge)[0]);

            a_idx = mesh.halfedge_vertices(halfedge)[1].idx();
            b_idx = mesh.halfedge_vertices(halfedge)[0].idx();

            Eigen::MatrixXd X, G_;
            X.resize(diamond_valence, 3);

            X.row(0) << C0[0], C0[1], C0[2];
            if (!mesh.is_boundary(f))
            {
                X.row(1) << C1[0], C1[1], C1[2];
            }
            X.row(nc) << A[0], A[1], A[2];
            X.row(nc + 1) << B[0], B[1], B[2];
            X.row(nc + 2) << F[0], F[1], F[2];

            if (mesh.is_boundary(f))
            {
                G_ = cc_bdy_gradient_operator(X);
            }
            else
            {
                G_ = cc_gradient_operator(X);
            }

            double volume_diamond = volume_tetrahedron_signed(B, A, F, C0) +
                                    volume_tetrahedron_signed(A, B, F, C1);

            // CO
            triplets.emplace_back(3 * row, c0_idx, G_(0, 0));
            triplets.emplace_back(3 * row + 1, c0_idx, G_(1, 0));
            triplets.emplace_back(3 * row + 2, c0_idx, G_(2, 0));

            //C1
            if (!mesh.is_boundary(f))
            {
                triplets.emplace_back(3 * row, c1_idx, G_(0, 1));
                triplets.emplace_back(3 * row + 1, c1_idx, G_(1, 1));
                triplets.emplace_back(3 * row + 2, c1_idx, G_(2, 1));
            }

            //Face Vertices

            int col = 0;
            for (int i = 0; i < 3; i++)
            {
                triplets.emplace_back(3 * row + i, a_idx, G_(i, nc + col));
            }
            col++;
            for (int i = 0; i < 3; i++)
            {
                triplets.emplace_back(3 * row + i, b_idx, G_(i, nc + col));
            }
            col++;
            for (int i = 0; i < 3; i++)
            {
                triplets.emplace_back(3 * row + i, f_idx, G_(i, nc + col));
            }

            for (int i = 0; i < 3; i++)
            {
                triplets_volume.emplace_back(3 * row + i, 3 * row + i,
                                             volume_diamond);
            }
            row++;
        }
    }
    V.resize(3 * row, 3 * row);
    V.setFromTriplets(triplets_volume.begin(), triplets_volume.end());
    G.resize(3 * row, n_v + n_f + n_c);
    G.setFromTriplets(triplets.begin(), triplets.end());
}

//===========================Gradient computations=============================================================

Eigen::MatrixXd cc_gradient_operator(Eigen::MatrixXd X)
{
    // We assume first two entries in X are the tips
    int n = X.rows() - 2;

    // local 3x3 matrices per diamond
    Eigen::Matrix3d A, scaledInverseA;

    // accumulate gradient and volumes
    double volumes = 0.0;
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3, n + 2);
    for (int i = 0; i < n; i++)
    {
        // indices of the two edge points
        int j1 = 2 + i;
        int j2 = 2 + (i + 1) % n;

        // build difference vectors from X
        Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, n + 2);
        D(0, 0) = -1.0;
        D(0, 1) = 1.0;
        D(1, 0) = -1.0;
        D(1, j1) = 1.0;
        D(2, 0) = -1.0;
        D(2, j2) = 1.0;

        // build 3x3 matrix
        A = D * X;

        // construct inverse(A)*determinant(A)
        scaledInverseA.col(0) = A.row(1).cross(A.row(2));
        scaledInverseA.col(1) = A.row(2).cross(A.row(0));
        scaledInverseA.col(2) = A.row(0).cross(A.row(1));

        // diamond volume == matrix determinant
        volumes += A.determinant();

        // accumulate gradient
        G += scaledInverseA * D;
    }

    G /= volumes;

    return G;
}

Eigen::MatrixXd cc_bdy_gradient_operator(Eigen::MatrixXd X)
{
    // We assuem first entry in in X is the tip
    int n = X.rows() - 1;

    // vector (0, 0, 1/n, ..., 1/n)
    Eigen::RowVectorXd one_over_n = Eigen::VectorXd::Constant(n + 1, 1.0 / n);
    one_over_n(0) = 0.0;

    // matrix for building difference vectors from X
    Eigen::MatrixXd D = Eigen::MatrixXd::Zero(3, n + 1);
    D.row(0) = one_over_n;
    D(0, 0) = -1.0;

    // local 3x3 matrices per diamond
    Eigen::Matrix3d A, scaledInverseA;

    // accumulate gradient and volumes
    double volumes = 0.0;
    Eigen::MatrixXd G = Eigen::MatrixXd::Zero(3, n + 1);
    for (int i = 0; i < n; i++)
    {
        // indices of the two edge points
        int j1 = 1 + i;
        int j2 = 1 + (i + 1) % n;

        // build difference vectors from X
        D.row(1) = -one_over_n;
        D(1, j1) += 1.0;
        D.row(2) = -one_over_n;
        D(2, j2) += 1.0;

        // build 3x3 matrix
        A = D * X;

        // construct inverse(A)*determinant(A)
        scaledInverseA.col(0) = A.row(1).cross(A.row(2));
        scaledInverseA.col(1) = A.row(2).cross(A.row(0));
        scaledInverseA.col(2) = A.row(0).cross(A.row(1));

        // diamond volume == matrix determinant
        volumes += A.determinant();

        // accumulate gradient
        G += scaledInverseA * D;
    }

    G /= volumes;

    return G;
}
