#include "LaplaceConstruction_3D.h"
#include <pmp/MatVec.h>
#include <cmath>
#include "diffgeo_3D.h"
#include "Diamond_3D.h"
 #include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include <random>
//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================
enum LaplaceMethods {
    Diamond = 0,
    Dual_Laplace = 1,
    PolySimpleLaplace = 2,
    Harmonic = 3
};

void setup_3D_stiffness_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S, int Laplace, int face_point,
                               int cell_point) {
    if (Laplace == Diamond) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        SparseMatrix G, V, Div, P, Pc, Pf;
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_diamond_gradient(mesh, G, V);
        Div = -G.transpose() * V;
        S = P.transpose() * Div * G * P;

    } else if (Laplace == PolySimpleLaplace) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        SparseMatrix Sf, G, Grad, V, Div, P, Pc, Pf;
        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_sandwich_stiffness_matrix(mesh, S, Sf, face_point, cell_point);
    } else if (Laplace == Dual_Laplace) {
        setup_3D_dual_laplace_libigl(mesh, S);
    }
}

//-----------------------------------------------------------------------------
void
setup_3D_mass_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M, int Laplace, int face_point,
                     int cell_point) {
    if (Laplace == Diamond) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        SparseMatrix P, Pc, Pf, M2;
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_diamond_mass_matrix(mesh, M2);
        M = P.transpose() * M2 * P;
        M /= 2.0;
    } else if (Laplace == PolySimpleLaplace) {
        setup_3D_sandwich_mass_matrix(mesh, M, face_point, cell_point);
        std::cout << "-----------------------" << std::endl;
        std::cout << "Sandwich mass: \n" << M << std::endl;
    } else if (Laplace == Dual_Laplace) {
        setup_3D_dual_mass_matrix_libigl(mesh, M);
    }
    double count = 0;
    for (int k = 0; k < M.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(M, k); it; ++it) {
            if (it.row() == it.col()) {
                double val = it.value();
            }
            count += it.value();
        }
    }
    std::cout << "Volume mass matrix: " << count << std::endl;
}

//-----------------------------------------------------------------------------

void setup_3D_dual_laplace_libigl(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S) {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> M;
    V.resize((int)mesh.n_vertices(), 3);
    F.resize((int)mesh.n_cells(), 4);
    for (auto v: mesh.vertices()) {
        VolumeMesh::PointT p = mesh.vertex(v);
        for (int i = 0; i < 3; i++) {
            V(v.idx(), i) = p[i];
        }
    }

    for (auto c: mesh.cells()) {
        int i = 0;
        for (auto v: mesh.cell_vertices(c)) {
            F(c.idx(), i) = v.idx();
            i++;
        }
    }

    dualLaplace(V, F, S, M);
}
//-----------------------------------------------------------------------------

void setup_3D_dual_mass_matrix_libigl(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M) {

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::SparseMatrix<double> S;
    V.resize((int)mesh.n_vertices(), 3);
    F.resize((int)mesh.n_cells(), 4);
    for (auto v: mesh.vertices()) {
        VolumeMesh::PointT p = mesh.vertex(v);
        for (int i = 0; i < 3; i++) {
            V(v.idx(), i) = p[i];
        }
    }

    for (auto c: mesh.cells()) {
        int i = 0;
        for (auto v: mesh.cell_vertices(c)) {
            F(c.idx(), i) = v.idx();
            i++;
        }
    }

    dualLaplace(V, F, S, M);
}

//-----------------------------------------------------------------------------

void setup_3D_cell_face_prolongation_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &P,
                                            Eigen::SparseMatrix<double> &P_cells,
                                            Eigen::SparseMatrix<double> &P_faces) {

    std::vector<T> triplet_face;
    std::vector<T> triplet_cells;

    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");


    int n_v = (int)mesh.n_vertices();
    int n_f = (int)mesh.n_faces();
    int n_c = (int)mesh.n_cells();


    for (auto v_it = mesh.v_iter(); v_it.valid(); ++v_it) {
        triplet_cells.emplace_back(T(v_it->idx(), v_it->idx(), 1.0));
        triplet_face.emplace_back(T(v_it->idx(), v_it->idx(), 1.0));
    }

    for (auto f_it = mesh.f_iter(); f_it.valid(); ++f_it) {
        Eigen::VectorXd face_w = f_w_prop[*f_it];
        auto f_v_it_pair = mesh.face_vertices(*f_it);
        int i = 0;
        for (auto f_v_it = f_v_it_pair.first;
             f_v_it != f_v_it_pair.second; ++f_v_it) {

            //add weight entries of the virtual face vertices in first prolongation matrix
            triplet_face.emplace_back(T(n_v + f_it->idx(), f_v_it->idx(), face_w(i)));
            i++;
        }
        //second prolongation matrix considers virtual face vertices as existent in the mesh
        triplet_cells.emplace_back(T(n_v + f_it->idx(), n_v + f_it->idx(), 1.0));
    }

    for (auto c_it = mesh.c_iter(); c_it.valid(); ++c_it) {

        Eigen::VectorXd cell_w = c_w_prop[*c_it];
        //face vertices are pushed back first during computation, so their weights make up the first part of the vector
        auto c_f_it_pair = mesh.cell_faces(*c_it);
        auto c_v_it_pair = mesh.cell_vertices(*c_it);

        int i = 0;
        for (auto c_f_it = c_f_it_pair.first;
             c_f_it != c_f_it_pair.second; ++c_f_it) {

            triplet_cells.emplace_back(T(n_v + n_f + c_it->idx(), n_v + c_f_it->idx(), cell_w(i)));
            i++;
        }
        for (auto c_v_it = c_v_it_pair.first;
             c_v_it != c_v_it_pair.second; ++c_v_it) {

            triplet_cells.emplace_back(T(n_v + n_f + c_it->idx(), c_v_it->idx(), cell_w(i)));
            i++;
        }

    }

    P_faces.resize(n_v + n_f, n_v);
    P_cells.resize(n_v + n_f + n_c, n_v + n_f);
    P_faces.setFromTriplets(triplet_face.begin(), triplet_face.end());
    P_cells.setFromTriplets(triplet_cells.begin(), triplet_cells.end());
    P = P_cells * P_faces;
}

//------------------------------------ Dual Laplacian--------------------------------------------------------

void
circumcenter(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, Eigen::Vector3d &cc) {
    const double l[3]{
            (b - c).squaredNorm(),
            (a - c).squaredNorm(),
            (a - b).squaredNorm()
    };

    const double ba[3]{l[0] * (l[1] + l[2] - l[0]), l[1] * (l[2] + l[0] - l[1]), l[2] * (l[0] + l[1] - l[2])};
    const double sum = ba[0] + ba[1] + ba[2];

    cc = (ba[0] / sum) * a + (ba[1] / sum) * b + (ba[2] / sum) * c;
}

void circumcenter(const Eigen::Matrix<double, 4, 3> &t, Eigen::Vector3d &c) {
    Eigen::Matrix3d A;
    Eigen::Vector3d b;

    const double n0 = t.row(0).squaredNorm();

    for (int k = 0; k < 3; ++k) {
        A.row(k) = t.row(k + 1) - t.row(0);
        b(k) = t.row(k + 1).squaredNorm() - n0;
    }

    c = 0.5 * A.fullPivHouseholderQr().solve(b);
}

double
volume(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Eigen::Vector3d &d) {
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;

    return A.determinant() / 6.;
}


void
dualLaplace(Eigen::MatrixXd &V, Eigen::MatrixXi &T, Eigen::SparseMatrix<double> &L,
            Eigen::SparseMatrix<double> &M) {

    const size_t nt = T.rows();
    const size_t nv = V.rows();

//     If volume negative switch orientation!
    for (int i = 0; i < (int)nt; i++) {
        Eigen::Matrix3d X;
        Eigen::Vector4i t = T.row(i);
        X.col(0) = V.row(t(1)) - V.row(t(0));
        X.col(1) = V.row(t(2)) - V.row(t(0));
        X.col(2) = V.row(t(3)) - V.row(t(0));

        if (X.determinant() < 0) {
            T(i, 0) = t(1);
            T(i, 1) = t(0);

            t = T.row(i);
            X.col(0) = V.row(t(1)) - V.row(t(0));
            X.col(1) = V.row(t(2)) - V.row(t(0));
            X.col(2) = V.row(t(3)) - V.row(t(0));
            if (X.determinant() < 0) {
                std::cout << " still zero!" << std::endl;
            }
        }
    }

    const int turn[4][4]
            {
                    {-1, 2,  3,  1},
                    {3,  -1, 0,  2},
                    {1,  3,  -1, 0},
                    {2,  0,  1,  -1}
            };

    auto getTet = [&](const int i, Eigen::Matrix<double, 4, 3> &t) {
        for (int k = 0; k < 4; ++k) {
            t.row(k) = V.row(T(i, k));
        }
    };

    std::vector<Eigen::Triplet<double>> tripL, tripM;

    Eigen::Vector3d cc;
    Eigen::Matrix<double, 4, 3> t;

    for (int k = 0; k < (int)nt; ++k) {
        getTet(k, t);
        circumcenter(t, cc);

        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (i != j) {
                    Eigen::Vector3d cf;
                    circumcenter(t.row(i), t.row(j), t.row(turn[i][j]), cf);

                    const Eigen::Vector3d ce = 0.5 * (t.row(i) + t.row(j));

                    const double vol = volume(t.row(i), ce, cf, cc);
                    const double wij = 6. * vol / (t.row(i) - t.row(j)).squaredNorm();

                    tripL.emplace_back(T(k, i), T(k, j), wij);
                    tripL.emplace_back(T(k, j), T(k, i), wij);

                    tripL.emplace_back(T(k, i), T(k, i), -wij);
                    tripL.emplace_back(T(k, j), T(k, j), -wij);

                    tripM.emplace_back(T(k, i), T(k, i), vol);
                    tripM.emplace_back(T(k, j), T(k, j), vol);
                }
            }
        }
    }

    L.resize((int)nv, (int)nv);
    M.resize((int)nv, (int)nv);

    L.setFromTriplets(tripL.begin(), tripL.end());
    M.setFromTriplets(tripM.begin(), tripM.end());

}


//--------------------------Sandwich Laplacian--------------------------------------------------------

void setup_3D_sandwich_stiffness_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S,
                                        Eigen::SparseMatrix<double> &S_tets, int face_point, int cell_point) {

    std::vector<T> triplet_list;
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");

    compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);

    triplet_list.reserve(mesh.n_cells() * 4 * 4);

    int n_v = (int)mesh.n_vertices();
    int n_f = (int)mesh.n_faces();
    int n_c = (int)mesh.n_cells();

    VolumeMesh::PointT i, j, k, l, n_ijk, n_ijl, n_jkl, n_ikl;

    for (auto c_it = mesh.c_iter(); c_it.valid(); ++c_it) {

        VolumeMesh::PointT c_bary = c_prop[*c_it];

        k = c_bary;
        auto c_f_it_pair = mesh.cell_faces(*c_it);
        for (auto c_f_it = c_f_it_pair.first; c_f_it != c_f_it_pair.second; ++c_f_it) {

            VolumeMesh::PointT f_bary = f_prop[*c_f_it];
            l = f_bary;
            auto f_he_it_pair = mesh.face_halfedges(*c_f_it);
            for (auto f_he_it = f_he_it_pair.first; f_he_it != f_he_it_pair.second; ++f_he_it) {

                auto v_from = mesh.from_vertex_handle(*f_he_it);
                auto v_to = mesh.to_vertex_handle(*f_he_it);

                i = mesh.vertex(v_from);
                j = mesh.vertex(v_to);

                n_ijk = compute_triangle_normal(i, j, k); // ijk inward
                n_ijl = compute_triangle_normal(i, l, j); // ilj inward
                n_jkl = compute_triangle_normal(j, l, k); // jlk inward
                n_ikl = compute_triangle_normal(i, k, l); // ikl inward

                double l_ij, l_ik, l_il, l_jk, l_kl, l_jl;
                double w_ij, w_ik, w_il, w_jk, w_kl, w_jl;

                int idx_i, idx_j, idx_k, idx_l;

                idx_i = v_from.idx();
                idx_j = v_to.idx();
                idx_l = c_f_it->idx() + n_v;
                idx_k = c_it->idx() + n_f + n_v;


                l_ij = (i - j).norm();
                l_ik = (i - k).norm();
                l_il = (i - l).norm();
                l_jk = (j - k).norm();
                l_jl = (j - l).norm();
                l_kl = (k - l).norm();


                w_ij = l_kl * (-n_ikl | n_jkl) / (n_ikl % n_jkl).norm();

                w_ik = l_jl * (-n_ijl | n_jkl) / (n_ijl % n_jkl).norm();
                w_il = l_jk * (-n_ijk | n_jkl) / (n_ijk % n_jkl).norm();

                w_jk = l_il * (-n_ijl | n_ikl) / (n_ijl % n_ikl).norm();
                w_jl = l_ik * (-n_ijk | n_ikl) / (n_ijk % n_ikl).norm();

                w_kl = l_ij * (-n_ijk | n_ijl) / (n_ijk % n_ijl).norm();


                triplet_list.emplace_back(T(idx_i, idx_j, w_ij));
                triplet_list.emplace_back(T(idx_j, idx_i, w_ij));

                triplet_list.emplace_back(T(idx_i, idx_k, w_ik));
                triplet_list.emplace_back(T(idx_k, idx_i, w_ik));

                triplet_list.emplace_back(T(idx_i, idx_l, w_il));
                triplet_list.emplace_back(T(idx_l, idx_i, w_il));

                triplet_list.emplace_back(T(idx_j, idx_k, w_jk));
                triplet_list.emplace_back(T(idx_k, idx_j, w_jk));

                triplet_list.emplace_back(T(idx_j, idx_l, w_jl));
                triplet_list.emplace_back(T(idx_l, idx_j, w_jl));

                triplet_list.emplace_back(T(idx_k, idx_l, w_kl));
                triplet_list.emplace_back(T(idx_l, idx_k, w_kl));


                triplet_list.emplace_back(T(idx_i, idx_i, -w_ij - w_ik - w_il));
                triplet_list.emplace_back(T(idx_j, idx_j, -w_ij - w_jk - w_jl));
                triplet_list.emplace_back(T(idx_k, idx_k, -w_ik - w_jk - w_kl));
                triplet_list.emplace_back(T(idx_l, idx_l, -w_il - w_jl - w_kl));
            }
        }
    }
    Eigen::SparseMatrix<double> P, Pt, P_faces, P_cells;
    setup_3D_cell_face_prolongation_matrix(mesh, P, P_cells, P_faces);

    S_tets.resize(n_v + n_f + n_c, n_v + n_f + n_c);
    S_tets.setFromTriplets(triplet_list.begin(), triplet_list.end());
    S_tets /= 6.0;

    Pt = P.transpose();
    S = Pt * S_tets * P;

}

//-----------------------------------------------------------------------------


void
setup_3D_sandwich_mass_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M, int face_point, int cell_point) {

    std::vector<T> triplet_list;
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");

    compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);

    triplet_list.reserve(mesh.n_cells() * 4 * 4);

    int n_v = (int)mesh.n_vertices();
    int n_f = (int)mesh.n_faces();
    int n_c = (int)mesh.n_cells();
    VolumeMesh::PointT i, j, k, l;

    for (auto c : mesh.cells()) {
        k = c_prop[c];
        for (auto f : mesh.cell_faces(c)) {
            l = f_prop[f];

            auto hf = mesh.face_halffaces(f)[0];
            if (mesh.is_boundary(hf)) {
                hf = mesh.face_halffaces(f)[1];
            }

            for (auto he : mesh.halfface_halfedges(hf)) {
                auto v_from = mesh.from_vertex_handle(he);
                auto v_to = mesh.to_vertex_handle(he);

                i = mesh.vertex(v_from);
                j = mesh.vertex(v_to);

                double volume;
                // is cell_point "left" or "right" from the face
                if (abs((k - c_prop[mesh.incident_cell(hf)]).norm()) < 0.00001) {
                    volume = volume_tetrahedron_signed(i, j, l, k);
                } else {
                    volume = volume_tetrahedron_signed(i, j, k, l);
                }
//                double volume = volume_tetrahedron_signed(i,j,l,k);
                if (volume < 0.0000001) std::cout << " negative Volume! \n";

                triplet_list.emplace_back(T(v_from.idx(), v_from.idx(), volume));
                triplet_list.emplace_back(T(v_to.idx(), v_to.idx(), volume));
                triplet_list.emplace_back(T(n_v + f.idx(), n_v + f.idx(), volume));
                triplet_list.emplace_back(T(n_v + n_f + c.idx(), n_v + n_f + c.idx(), volume));
            }
        }
    }

    Eigen::SparseMatrix<double> M_f, P, Pt, P_faces, P_cells;
    setup_3D_cell_face_prolongation_matrix(mesh, P, P_cells, P_faces);
    M_f.resize(n_v + n_f + n_c, n_v + n_f + n_c);
    M_f.setFromTriplets(triplet_list.begin(), triplet_list.end());
    M_f /= 4.0;
    Pt = P.transpose();
    M = Pt * M_f * P;

}
