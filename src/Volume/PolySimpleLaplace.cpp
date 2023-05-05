
#include "PolySimpleLaplace.h"
#include "diffgeo_3D.h"
#include "LaplaceConstruction_3D.h"

//--------------------------Polysimple Laplacian--------------------------------------------------------

void setup_3D_polysimple_stiffness_matrix(VolumeMesh& mesh,
                                          Eigen::SparseMatrix<double>& S,
                                          Eigen::SparseMatrix<double>& S_tets,
                                          int face_point, int cell_point)
{
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

    for (auto c_it = mesh.c_iter(); c_it.valid(); ++c_it)
    {
        VolumeMesh::PointT c_bary = c_prop[*c_it];

        k = c_bary;
        auto c_f_it_pair = mesh.cell_faces(*c_it);
        for (auto c_f_it = c_f_it_pair.first; c_f_it != c_f_it_pair.second;
             ++c_f_it)
        {
            VolumeMesh::PointT f_bary = f_prop[*c_f_it];
            l = f_bary;
            auto f_he_it_pair = mesh.face_halfedges(*c_f_it);
            for (auto f_he_it = f_he_it_pair.first;
                 f_he_it != f_he_it_pair.second; ++f_he_it)
            {
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

void setup_3D_polysimple_mass_matrix(VolumeMesh& mesh,
                                     Eigen::SparseMatrix<double>& M,
                                     int face_point, int cell_point)
{
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

    for (auto c : mesh.cells())
    {
        k = c_prop[c];
        for (auto f : mesh.cell_faces(c))
        {
            l = f_prop[f];

            auto hf = mesh.face_halffaces(f)[0];
            if (mesh.is_boundary(hf))
            {
                hf = mesh.face_halffaces(f)[1];
            }

            for (auto he : mesh.halfface_halfedges(hf))
            {
                auto v_from = mesh.from_vertex_handle(he);
                auto v_to = mesh.to_vertex_handle(he);

                i = mesh.vertex(v_from);
                j = mesh.vertex(v_to);

                double volume;
                // is cell_point "left" or "right" from the face
                if (abs((k - c_prop[mesh.incident_cell(hf)]).norm()) < 0.00001)
                {
                    volume = volume_tetrahedron_signed(i, j, l, k);
                }
                else
                {
                    volume = volume_tetrahedron_signed(i, j, k, l);
                }

                triplet_list.emplace_back(
                    T(v_from.idx(), v_from.idx(), volume));
                triplet_list.emplace_back(T(v_to.idx(), v_to.idx(), volume));
                triplet_list.emplace_back(
                    T(n_v + f.idx(), n_v + f.idx(), volume));
                triplet_list.emplace_back(
                    T(n_v + n_f + c.idx(), n_v + n_f + c.idx(), volume));
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
