#include "pmp/algorithms/SurfaceSubdivision.h"
#include "SECLaplace.h"
#include "DisneyLaplace.h"
#include "LaplaceConstruction.h"
#include "Diamond_2D.h"
#include "diffgeo.h"

//=============================================================================

using SurfaceMesh = pmp::SurfaceMesh;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

float sec_laplace_lambda_ = 1.0;

int sec_laplace_lvl = 2;
//=============================================================================
enum Basis
{

    CatmullClark = 0,
    AdaptedCC = 1,
    ReducingCC = 2,
    LinearInterpol = 3,
    Interpol_BF =4,
    Diamond_BF =5,

};

void setup_smooth_laplace_matrices(pmp::SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &S,
                               Eigen::SparseMatrix<double> &M)
{

    SurfaceMesh linear_interpolation_mesh = mesh;
    FaceProperty<pmp::Point> area_points;
    FaceProperty<Eigen::VectorXd> area_weights;
        if(!linear_interpolation_mesh.has_face_property("f:point")){
         area_points = linear_interpolation_mesh.add_face_property<pmp::Point>("f:point");
    }else{
        area_points = linear_interpolation_mesh.get_face_property<pmp::Point>("f:point");
    }

    if(!linear_interpolation_mesh.has_face_property("f:weights"))
    {
         area_weights = linear_interpolation_mesh.add_face_property<Eigen::VectorXd>("f:weights");
    }else{
         area_weights = linear_interpolation_mesh.get_face_property<Eigen::VectorXd>("f:weights");
    }


    setup_face_point_properties(linear_interpolation_mesh,2);

    //Final cc Prolongation matrix
    Eigen::SparseMatrix<double> P, P_cc, P_init, S_subdiv, M_subdiv;
    //Individual cc Prolongation matrices for each lvl
    // Ordering: 0-nv vertices, nv-nv+ne edge vertices, nv+ne+nf face vertices
     // lvl-1 since one subdivision step is always the interpolation
//
    setup_prolongation_matrix(linear_interpolation_mesh,P_init);
    insert_points(linear_interpolation_mesh, 2);
    setup_mod_butterfly_P_matrix(linear_interpolation_mesh, P_cc, sec_laplace_lvl-1 );
    P = P_cc*P_init;

//    setup_mod_butterfly_P_matrix(linear_interpolation_mesh, P_cc, sec_laplace_lvl );
//    P = P_cc;
//    setup_disney_laplace_operator(linear_interpolation_mesh,S_subdiv);
//    setup_disney_mass_matrix(linear_interpolation_mesh, M_subdiv);
//    setup_sandwich_stiffness_matrix(linear_interpolation_mesh, S_subdiv, 2);
//    setup_sandwich_mass_matrix(linear_interpolation_mesh, M_subdiv, 2);

    Eigen::SparseMatrix<double> G,D;

    //MULLIFICATION
    double delta = 0.0;
    for(auto e: linear_interpolation_mesh.edges()){
        delta += linear_interpolation_mesh.edge_length(e);
    }
    delta/= linear_interpolation_mesh.n_edges();
    delta*= 1e-5;

    setup_diamond_gradient_divergence_intrinsic_mullification(linear_interpolation_mesh, G, D, delta);
    setup_diamond_mass_matrix_intrinsic_mullification(linear_interpolation_mesh,M_subdiv, delta);

    S_subdiv = D*G;
    linear_interpolation_mesh.remove_face_property(area_points);
    linear_interpolation_mesh.remove_face_property(area_weights);


    S = P.transpose() * S_subdiv * P;
    M = P.transpose() * M_subdiv * P;
}
//=============================================================================

void setup_sec_Laplace_matrix(pmp::SurfaceMesh &mesh,
                              Eigen::SparseMatrix<double> &L)
{

    Eigen::MatrixXd Pos(mesh.n_vertices(), 3), Pos_subdiv;
    for (auto v : mesh.vertices())
    {
        pmp::Point p = mesh.position(v);
        Pos(v.idx(), 0) = p[0];
        Pos(v.idx(), 1) = p[1];
        Pos(v.idx(), 2) = p[2];
    }
    SurfaceMesh catmull_clark_mesh = mesh;
    SurfaceMesh linear_interpolation_mesh = mesh;
    SurfaceMesh cc_control_mesh = mesh;

    //        SurfaceSubdivision divider = SurfaceSubdivision(cc_control_mesh);
    for (int i = 0; i < sec_laplace_lvl; i++)
    {
        //                divider.catmull_clark();
        //        linear_interpolation_catmull_clark(cc_control_mesh);
        linear_interpolation_catmull_clark(linear_interpolation_mesh);
    }

    std::cout << "Subdiv mesh vertices: " << catmull_clark_mesh.n_vertices()
              << std::endl;
    std::cout << "Subdiv mesh control vertices: "
              << linear_interpolation_mesh.n_vertices() << std::endl;

    Eigen::SparseMatrix<double> L_subdiv, A_0;
    std::vector<Eigen::SparseMatrix<double>> Al;
    //
    setup_interpolating_sec_A0_matrix(catmull_clark_mesh, A_0, sec_laplace_lvl,
                                      Al);
    //            setup_bdry_cc_A0_matrix(catmull_clark_mesh,A_0,sec_laplace_lvl,Al);
    //        setup_sec_A0_matrix(catmull_clark_mesh, A_0, sec_laplace_lvl,Al);

    //    setup_sandwich_stiffness_matrix(catmull_clark_mesh, L_subdiv);

    //        setup_sandwich_stiffness_matrix(linear_interpolation_mesh, L_subdiv);

    //        setup_disney_laplace_operator(linear_interpolation_mesh, L_subdiv);
    //        setup_disney_laplace_operator(catmull_clark_mesh, L_subdiv);
    //    disney_laplace_lambda_ = sec_laplace_lambda_;
    L = A_0.transpose() * L_subdiv * A_0;

    //    Pos_subdiv = A_0 * Pos;
    //    for (auto v : cc_control_mesh.vertices())
    //    {
    //        Eigen::Vector3d sd, sd_test;
    //        sd_test << cc_control_mesh.position(v)[0],
    //            cc_control_mesh.position(v)[1], cc_control_mesh.position(v)[2];
    //        sd << Pos_subdiv(v.idx(), 0), Pos_subdiv(v.idx(), 1),
    //            Pos_subdiv(v.idx(), 2);
    //        if (abs((sd - sd_test).norm()) > 0.000001)
    //        {
    //            std::cout << "Subdiv vertex " << v.idx() << " : ("
    //                      << cc_control_mesh.position(v)[0] << ","
    //                      << cc_control_mesh.position(v)[1] << ","
    //                      << cc_control_mesh.position(v)[2]
    //                      << ")  prolongated Vertex : (" << Pos_subdiv(v.idx(), 0)
    //                      << "," << Pos_subdiv(v.idx(), 1) << ","
    //                      << Pos_subdiv(v.idx(), 2) << ")" << std::endl;
    //        }
    //    }
}

void setup_sec_mass_matrix(pmp::SurfaceMesh &mesh,
                           Eigen::SparseMatrix<double> &M)
{

    SurfaceMesh subdiv_mesh = mesh;
    SurfaceMesh lin_interpolation_mesh = mesh;
    SurfaceMesh cc_control_mesh = mesh;
    Eigen::MatrixXd Pos(mesh.n_vertices(), 3), Pos_subdiv;
    for (auto v : mesh.vertices())
    {
        pmp::Point p = mesh.position(v);
        Pos(v.idx(), 0) = p[0];
        Pos(v.idx(), 1) = p[1];
        Pos(v.idx(), 2) = p[2];
    }

    Eigen::SparseMatrix<double> M_subdiv, A_0, P, M_;
    //    for (int i = 0; i < sec_laplace_lvl; i++) {
    SurfaceSubdivision divider = SurfaceSubdivision(cc_control_mesh);
    for (int i = 0; i < sec_laplace_lvl; i++)
    {
        divider.catmull_clark();
        linear_interpolation_catmull_clark(lin_interpolation_mesh);
    }
    std::cout << sec_laplace_lvl << std::endl;
    std::vector<Eigen::SparseMatrix<double>> Al;
    //    setup_bdry_cc_A0_matrix(subdiv_mesh,A_0,sec_laplace_lvl,Al);
    setup_interpolating_sec_A0_matrix(subdiv_mesh, A_0, sec_laplace_lvl, Al);
    //            setup_sec_A0_matrix(subdiv_mesh, A_0, sec_laplace_lvl,Al);
    //    setup_sandwich_mass_matrix(subdiv_mesh, M_subdiv, 0);

    //     setup_poly_mass_matrix(lin_interpolation_mesh, M_subdiv);
    setup_sandwich_mass_matrix(lin_interpolation_mesh, M_subdiv, 0);
    //    setup_sandwich_mass_matrix(cc_control_mesh, M_subdiv, 0);

    //    std::cout << A_0 << std::endl;
    //    if (sec_laplace_lvl == 1)
    //    {
    //        lump_matrix(M_subdiv);
    //    }
    //    lump_matrix(M_subdiv);

    Eigen::SparseMatrix<double> M_l, M_l1, Slt, Sl;
    //    M_l = M_subdiv;
    //    for (int i = 0; i < sec_laplace_lvl; i++)
    //    {
    //        Slt = Al[sec_laplace_lvl - 1 - i].transpose();
    //        Sl = Al[sec_laplace_lvl - 1 - i];
    //        M_l1 = Slt * M_l * Sl;
    //        //        std::cout << Slt << std::endl;
    //        //        std::cout << M_l1 << std::endl;
    ////                lump_matrix(M_l1);
    //        M_l = M_l1;
    //        if (sec_laplace_lvl > 2 &&
    //            i % (sec_laplace_lvl - 1) == sec_laplace_lvl - 1)
    //        {
    //            lump_matrix(M_l);
    //        }
    //        else if (sec_laplace_lvl == 2 && i == 0)
    //        {
    //            lump_matrix(M_l);
    //        }
    //    }
    M = A_0.transpose() * M_subdiv * A_0;
    //    M = M_l;
    lump_matrix(M);
    //    std::cout << M << std::endl;
    //    Pos_subdiv = A_0 * Pos;
    //        for (auto v : cc_control_mesh.vertices())
    //        {
    //            Eigen::Vector3d sd, sd_test;
    //            sd_test << cc_control_mesh.position(v)[0],
    //                cc_control_mesh.position(v)[1], cc_control_mesh.position(v)[2];
    //            sd << Pos_subdiv(v.idx(), 0), Pos_subdiv(v.idx(), 1),
    //                Pos_subdiv(v.idx(), 2);
    //                    if (abs((sd - sd_test).norm()) > 0.000001) {
    //            std::cout << "Subdiv vertex " << v.idx() << " : ("
    //                      << cc_control_mesh.position(v)[0] << ","
    //                      << cc_control_mesh.position(v)[1] << ","
    //                      << cc_control_mesh.position(v)[2]
    //                      << ")  prolongated Vertex : (" << Pos_subdiv(v.idx(), 0)
    //                      << "," << Pos_subdiv(v.idx(), 1) << ","
    //                      << Pos_subdiv(v.idx(), 2) << ")" << std::endl;
    //                    }
    //        }
}

void setup_sec_gradient_operator(pmp::SurfaceMesh &mesh,
                                 Eigen::SparseMatrix<double> &G)
{
    SurfaceMesh subdiv_mesh = mesh;
    Eigen::SparseMatrix<double> Grad_subdiv, A_0;
    std::vector<Eigen::SparseMatrix<double>> Al;
    setup_cc_P_matrix(subdiv_mesh, A_0, sec_laplace_lvl, Al);
    poly_laplace_lambda_ = sec_laplace_lambda_;
    setup_poly_gradient_operator(subdiv_mesh, Grad_subdiv);
    G = 0.5 * Grad_subdiv * A_0;
}

void setup_cc_P_matrix(pmp::SurfaceMesh &mesh_, Eigen::SparseMatrix<double> &A0,
                       int lvl, std::vector<Eigen::SparseMatrix<double>> &Al)
{

    Eigen::SparseMatrix<double> S_l;
    Al.resize(lvl);
    for (int i = 0; i < lvl; i++)
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.n_vertices() * 4 * 4);

        pmp::VertexProperty<pmp::Point> points_ =
            mesh_.vertex_property<pmp::Point>("v:point");

        // reserve memory
        size_t nv = mesh_.n_vertices();
        size_t ne = mesh_.n_edges();
        size_t nf = mesh_.n_faces();
        mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);
        S_l.resize(nv + ne + nf, nv);

        // get properties
        auto vpoint = mesh_.add_vertex_property<pmp::Point>("catmull:vpoint");
        auto epoint = mesh_.add_edge_property<pmp::Point>("catmull:epoint");
        auto fpoint = mesh_.add_face_property<pmp::Point>("catmull:fpoint");

        // compute face vertices
        for (auto f : mesh_.faces())
        {
            pmp::Point p(0.0, 0.0, 0.0);
            pmp::Scalar c(0.0);
            double val_f = mesh_.valence(f);
            for (auto v : mesh_.vertices(f))
            {
                p += points_[v];
                ++c;
                triplets.emplace_back(nv + ne + f.idx(), v.idx(),
                                      (1.0 / val_f));
            }
            p /= c;
            fpoint[f] = p;
        }

        // compute edge vertices
        for (auto e : mesh_.edges())
        {
            // boundary or feature edge?
            if (mesh_.is_boundary(e))
            {
                epoint[e] = 0.5 * (points_[mesh_.vertex(e, 0)] +
                                   points_[mesh_.vertex(e, 1)]);

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.5);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.5);
            }

            // interior edge
            else
            {
                pmp::Point p(0.0, 0.0, 0.0);
                p += points_[mesh_.vertex(e, 0)];
                p += points_[mesh_.vertex(e, 1)];
                p += fpoint[mesh_.face(e, 0)];
                p += fpoint[mesh_.face(e, 1)];
                p *= 0.25;
                epoint[e] = p;

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                auto f0 = mesh_.face(e, 0);
                auto f1 = mesh_.face(e, 1);
                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.25);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.25);

                double val_f0 = mesh_.valence(f0);
                double val_f1 = mesh_.valence(f1);

                for (auto v : mesh_.vertices(f0))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f0) * 0.25);
                }
                for (auto v : mesh_.vertices(f1))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f1) * 0.25);
                }
            }
        }

        // compute new positions for old vertices
        for (auto v : mesh_.vertices())
        {
            // isolated vertex?
            if (mesh_.is_isolated(v))
            {
                vpoint[v] = points_[v];

                triplets.emplace_back(v.idx(), v.idx(), 1.0);
            }

            // boundary vertex?
            else if (mesh_.is_boundary(v))
            {
                auto h1 = mesh_.halfedge(v);
                auto h0 = mesh_.prev_halfedge(h1);

                pmp::Point p = points_[v];
                p *= 6.0;
                p += points_[mesh_.to_vertex(h1)];
                p += points_[mesh_.from_vertex(h0)];
                p *= 0.125;

                vpoint[v] = p;

                auto v_to = mesh_.to_vertex(h1);
                auto v_from = mesh_.from_vertex(h0);

                triplets.emplace_back(v.idx(), v.idx(), 6.0 * 0.125);
                triplets.emplace_back(v.idx(), v_to.idx(), 0.125);
                triplets.emplace_back(v.idx(), v_from.idx(), 0.125);
            }

            // interior vertex
            else
            {
                // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

                const double k = mesh_.valence(v);
                pmp::Point p(0, 0, 0);

                for (auto vv : mesh_.vertices(v))
                {
                    p += points_[vv];
                    triplets.emplace_back(v.idx(), vv.idx(), 1.0 / (k * k));
                }

                for (auto f : mesh_.faces(v))
                {
                    p += fpoint[f];
                    double f_val = mesh_.valence(f);
                    for (auto vv : mesh_.vertices(f))
                    {
                        triplets.emplace_back(v.idx(), vv.idx(),
                                              1.0 / (k * k * f_val));
                    }
                }

                p /= (k * k);

                p += ((k - 2.0) / k) * points_[v];

                triplets.emplace_back(v.idx(), v.idx(), (k - 2.0) / k);

                vpoint[v] = p;
            }
        }

        // assign new positions to old vertices
        for (auto v : mesh_.vertices())
        {
            points_[v] = vpoint[v];
        }

        // split edges
        for (auto e : mesh_.edges())
        {
            // normal edge
            mesh_.insert_vertex(e, epoint[e]);
        }

        // split faces
        for (auto f : mesh_.faces())
        {
            auto h0 = mesh_.halfedge(f);
            mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

            auto h1 = mesh_.next_halfedge(h0);
            mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }

        // clean-up properties
        mesh_.remove_vertex_property(vpoint);
        mesh_.remove_edge_property(epoint);
        mesh_.remove_face_property(fpoint);

        S_l.setFromTriplets(triplets.begin(), triplets.end());
        Al[i] = S_l;
        if (i == 0)
        {
            A0 = S_l;
        }
        else
        {
            std::cout << "A0: " << A0.rows() << "x" << A0.cols() << std::endl;
            std::cout << "Sl: " << S_l.rows() << "x" << S_l.cols() << std::endl;

            Eigen::SparseMatrix<double> S_l1 = S_l * A0;
            A0 = S_l1;
        }
    }
}

void setup_cc_retain_P_matrix(pmp::SurfaceMesh &mesh_,
                              Eigen::SparseMatrix<double> &A0, int lvl)
{
    Eigen::SparseMatrix<double> S_l;
    for (int i = 0; i < lvl; i++)
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.n_vertices() * 4 * 4);

        pmp::VertexProperty<pmp::Point> points_ =
            mesh_.vertex_property<pmp::Point>("v:point");

        // reserve memory
        size_t nv = mesh_.n_vertices();
        size_t ne = mesh_.n_edges();
        size_t nf = mesh_.n_faces();
        mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);
        S_l.resize(nv + ne + nf, nv);

        // get properties
        auto vpoint = mesh_.add_vertex_property<pmp::Point>("catmull:vpoint");
        auto epoint = mesh_.add_edge_property<pmp::Point>("catmull:epoint");
        auto fpoint = mesh_.add_face_property<pmp::Point>("catmull:fpoint");
        auto v_orig = mesh_.get_vertex_property<bool>("v:feature");

        // compute face vertices
        for (auto f : mesh_.faces())
        {
            pmp::Point p(0.0, 0.0, 0.0);
            pmp::Scalar c(0.0);
            double val_f = mesh_.valence(f);
            for (auto v : mesh_.vertices(f))
            {
                p += points_[v];
                ++c;
                triplets.emplace_back(nv + ne + f.idx(), v.idx(),
                                      (1.0 / val_f));
            }
            p /= c;
            fpoint[f] = p;
        }

        // compute edge vertices
        for (auto e : mesh_.edges())
        {
            // boundary or feature edge?
            if (mesh_.is_boundary(e))
            {
                epoint[e] = 0.5 * (points_[mesh_.vertex(e, 0)] +
                                   points_[mesh_.vertex(e, 1)]);

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.5);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.5);
            }

            // interior edge
            else
            {
                pmp::Point p(0.0, 0.0, 0.0);
                p += points_[mesh_.vertex(e, 0)];
                p += points_[mesh_.vertex(e, 1)];
                p += fpoint[mesh_.face(e, 0)];
                p += fpoint[mesh_.face(e, 1)];
                p *= 0.25;
                epoint[e] = p;

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                auto f0 = mesh_.face(e, 0);
                auto f1 = mesh_.face(e, 1);
                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.25);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.25);

                double val_f0 = mesh_.valence(f0);
                double val_f1 = mesh_.valence(f1);

                for (auto v : mesh_.vertices(f0))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f0) * 0.25);
                }
                for (auto v : mesh_.vertices(f1))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f1) * 0.25);
                }
            }
        }

        // compute new positions for old vertices
        for (auto v : mesh_.vertices())
        {
            // isolated vertex?
            if (mesh_.is_isolated(v) || v_orig[v])
            {
                vpoint[v] = points_[v];
                triplets.emplace_back(v.idx(), v.idx(), 1.0);
            }

            // boundary vertex?
            else if (mesh_.is_boundary(v))
            {
                auto h1 = mesh_.halfedge(v);
                auto h0 = mesh_.prev_halfedge(h1);

                pmp::Point p = points_[v];
                p *= 6.0;
                p += points_[mesh_.to_vertex(h1)];
                p += points_[mesh_.from_vertex(h0)];
                p *= 0.125;

                vpoint[v] = p;

                auto v_to = mesh_.to_vertex(h1);
                auto v_from = mesh_.from_vertex(h0);

                triplets.emplace_back(v.idx(), v.idx(), 6.0 * 0.125);
                triplets.emplace_back(v.idx(), v_to.idx(), 0.125);
                triplets.emplace_back(v.idx(), v_from.idx(), 0.125);
            }

            // interior vertex
            else
            {
                // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

                const double k = mesh_.valence(v);
                pmp::Point p(0, 0, 0);

                for (auto vv : mesh_.vertices(v))
                {
                    p += points_[vv];
                    triplets.emplace_back(v.idx(), vv.idx(), 1.0 / (k * k));
                }

                for (auto f : mesh_.faces(v))
                {
                    p += fpoint[f];
                    double f_val = mesh_.valence(f);
                    for (auto vv : mesh_.vertices(f))
                    {
                        triplets.emplace_back(v.idx(), vv.idx(),
                                              1.0 / (k * k * f_val));
                    }
                }

                p /= (k * k);

                p += ((k - 2.0) / k) * points_[v];

                triplets.emplace_back(v.idx(), v.idx(), (k - 2.0) / k);

                vpoint[v] = p;
            }
        }

        // assign new positions to old vertices
        for (auto v : mesh_.vertices())
        {
            points_[v] = vpoint[v];
        }

        // split edges
        for (auto e : mesh_.edges())
        {
            // normal edge
            auto he = mesh_.insert_vertex(e, epoint[e]);
            v_orig[mesh_.to_vertex(he)] = false;
        }

        // split faces
        for (auto f : mesh_.faces())
        {
            auto h0 = mesh_.halfedge(f);
            mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

            auto h1 = mesh_.next_halfedge(h0);
            auto he = mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);
            v_orig[mesh_.to_vertex(he)] = false;

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }

        // clean-up properties
        mesh_.remove_vertex_property(vpoint);
        mesh_.remove_edge_property(epoint);
        mesh_.remove_face_property(fpoint);

        S_l.setFromTriplets(triplets.begin(), triplets.end());
        if (i == 0)
        {
            A0 = S_l;
        }
        else
        {
            std::cout << "A0: " << A0.rows() << "x" << A0.cols() << std::endl;
            std::cout << "Sl: " << S_l.rows() << "x" << S_l.cols() << std::endl;

            Eigen::SparseMatrix<double> S_l1 = S_l * A0;
            A0 = S_l1;
        }
    }
}

void setup_mod_butterfly_P_matrix(pmp::SurfaceMesh &mesh_,
                                  Eigen::SparseMatrix<double> &A0, int lvl) {

    Eigen::SparseMatrix<double> S_l;
    for (int i = 0; i < lvl; i++)
    {
        size_t nv = mesh_.n_vertices();
        size_t ne = mesh_.n_edges();
        size_t nf = mesh_.n_faces();
        // reserve memory
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.n_vertices() * 4 * 4);
        S_l.resize(nv + ne , nv);

        pmp::VertexProperty<pmp::Point> points_ = mesh_.vertex_property<pmp::Point>("v:point");
        mesh_.reserve(nv + ne, 2 * ne + 3 * nf, 4 * nf);

        // add properties
        auto epoint = mesh_.add_edge_property<Point>("bf:epoint");

        // compute vertex positions
        for (auto v : mesh_.vertices())
        {
            triplets.emplace_back(v.idx(),v.idx(),1.0);
        }

        // compute edge positions
        for (auto e : mesh_.edges())
        {
            // boundary edge?
            if (mesh_.is_boundary(e))
            {
                auto h0 = mesh_.halfedge(e, 0);
                auto h1 = mesh_.halfedge(e, 1);
                Halfedge h_1,h2;
                if(mesh_.is_boundary(h0)){
                    h2=mesh_.next_halfedge(h0);
                    h_1=mesh_.prev_halfedge(h0);
                    if(!mesh_.is_boundary(h2) || !mesh_.is_boundary( h_1)){
                        std::cout << "Wrong boundaries chosen" << std::endl;
                    }
                }else{
                    h2=mesh_.next_halfedge(h1);
                    h_1=mesh_.prev_halfedge(h1);
                    if(!mesh_.is_boundary(h2) || !mesh_.is_boundary( h_1)){
                        std::cout << "Wrong boundaries chosen" << std::endl;
                    }
                }

                Vertex v0 = mesh_.to_vertex(h0);
                Vertex v1 = mesh_.to_vertex(h1);
//
//                triplets.emplace_back(nv+e.idx(),v0.idx(),0.5);
//                triplets.emplace_back(nv+e.idx(), v1.idx(), 0.5);
//                epoint[e] = (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]) *0.5;

                triplets.emplace_back(nv+e.idx(),v0.idx(),9.0/16.0);
                triplets.emplace_back(nv+e.idx(), v1.idx(), 9.0/16.0);

                triplets.emplace_back(nv+e.idx(), mesh_.from_vertex(h_1).idx(), -1.0/16.0);
                triplets.emplace_back(nv+e.idx(), mesh_.to_vertex(h2).idx(), -1.0/16.0);

                epoint[e] = (points_[mesh_.vertex(e, 0)] +
                             points_[mesh_.vertex(e, 1)]) *
                            9.0/16.0-(points_[mesh_.from_vertex(h_1)] +
                             points_[mesh_.to_vertex(h2)])*1.0/16.0;
            }
            // interior edge
            else
            {

                auto h0 = mesh_.halfedge(e, 0);
                auto h1 = mesh_.halfedge(e, 1);
                Vertex v0 = mesh_.to_vertex(h0);
                Vertex v1 = mesh_.to_vertex(h1);
                //to regular vertices => normal butterfly
                if(mesh_.valence(v0) == 6 && mesh_.valence(v1)==6){
                    Point p = points_[v0];
                    p += points_[v1];
                    p *= 8.0;
                    triplets.emplace_back(nv+e.idx(),v0.idx(),0.5);
                    triplets.emplace_back(nv+e.idx(), v1.idx(), 0.5);

                    p += 2.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h0))];
                    p += 2.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h1))];

                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h0)).idx(),0.125);
                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h1)).idx(),0.125);

                    auto h0o0 = mesh_.opposite_halfedge(mesh_.next_halfedge(h0));
                    auto h0o1=mesh_.opposite_halfedge(mesh_.prev_halfedge(h0));

                    p -= points_[mesh_.to_vertex(mesh_.next_halfedge(h0o0))];
                    p -= points_[mesh_.to_vertex(mesh_.next_halfedge(h0o1))];

                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h0o0)).idx(),-0.0625);
                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h0o1)).idx(),-0.0625);

                    auto h1o0 = mesh_.opposite_halfedge(mesh_.next_halfedge(h1));
                    auto h1o1=mesh_.opposite_halfedge(mesh_.prev_halfedge(h1));

                    p -= points_[mesh_.to_vertex(mesh_.next_halfedge(h1o0))];
                    p -= points_[mesh_.to_vertex(mesh_.next_halfedge(h1o1))];

                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h1o0)).idx(),-0.0625);
                    triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h1o1)).idx(),-0.0625);

                    p *= 0.0625;
                    epoint[e] = p;
                }//two extraordinary vertices
                else if(mesh_.valence(v0) != 6 && mesh_.valence(v1) != 6 )
                {
                    Halfedge h_reg, h_ext;
                    Vertex v_reg, v_ext;
                    Point p(0.0);
                    for(int j = 0; j<2;j++){
                        if(j==0){
                            h_reg = h0;
                            h_ext = h1;
                            v_reg = v0;
                            v_ext = v1;
                        }else{
                            h_reg = h1;
                            h_ext = h0;
                            v_reg = v1;
                            v_ext = v0;
                        }
                         if(mesh_.valence(v_ext) ==3){
                            p+=0.5*0.75*points_[v_ext];
                            p+=0.5*5.0/12.0*points_[v_reg] ;
                            triplets.emplace_back(nv+e.idx(),v_reg.idx(),0.5*5.0/12.0);
                            triplets.emplace_back(nv+e.idx(),v_ext.idx(),0.5*0.75);

                            p -= 1.0/12.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h_reg))]*0.5;
                            p -= 1.0/12.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h_ext))]*0.5;

                            triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h_ext)).idx(),-1.0/12.0*0.5);
                            triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h_reg)).idx(),-1.0/12.0*0.5);


                        }else if(mesh_.valence(v_ext)==4){
                            p+=0.5*0.75*points_[v_ext];
                            p+=0.5*3.0/8.0*points_[v_reg] ;
                            triplets.emplace_back(nv+e.idx(),v_reg.idx(),0.5*3.0/8.0);
                            triplets.emplace_back(nv+e.idx(),v_ext.idx(),0.5*0.75);

                            auto h1o0 = mesh_.opposite_halfedge(mesh_.next_halfedge(h_ext));
                            p -= 0.5*1.0/8.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h1o0))];
                            triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h1o0)).idx(),-0.5*1.0/8.0);

                        }else{
                            double k = 0.0;
                            double val = mesh_.valence(v_ext);
                            p+=0.5*0.75*points_[v_ext];
                            triplets.emplace_back(nv+e.idx(),v_ext.idx(),0.5*0.75);


                            //w_j =1/4 + cos(2pij/val) + 1/2 cos(4pij/val))/val
                            p += 0.5*1.75/val*points_[v_reg];
                            triplets.emplace_back(nv+e.idx(),v_reg.idx(),0.5*1.75/val);
                            auto h_next = mesh_.opposite_halfedge(mesh_.prev_halfedge(h_reg));
                            while(h_next != h_reg){
                                k ++;
                                double wj = (0.25 + cos(2.0*M_PI*k/val) + 0.5 *cos(4.0*M_PI*k/val))/val;
                                wj*=0.5;
                                p+=wj*points_[mesh_.to_vertex(h_next)];
                                triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(h_next).idx(),wj);
                                h_next = mesh_.opposite_halfedge(mesh_.prev_halfedge(h_next));
                            }
                        }
                    }
                    epoint[e] = p;
                }// One regular vertex, one extraordinary
                else
                {
                    Halfedge h_reg, h_ext;
                    Vertex v_reg, v_ext;
                    Point p;
                    if(mesh_.valence(v1) != 6){
                        h_reg = h0;
                        h_ext = h1;
                        v_reg = v0;
                        v_ext = v1;
                    }else{
                        h_reg = h1;
                        h_ext = h0;
                        v_reg = v1;
                        v_ext = v0;
                    }
                    p=0.75*points_[v_ext];
                    triplets.emplace_back(nv+e.idx(),v_ext.idx(),0.75);
                    if(mesh_.valence(v_ext) ==3){
                        p+=5.0/12.0*points_[v_reg] ;
                        triplets.emplace_back(nv+e.idx(),v_reg.idx(),5.0/12.0);

                        p -=1.0/12.0* points_[mesh_.to_vertex(mesh_.next_halfedge(h_reg))];
                        p -=1.0/12.0* points_[mesh_.to_vertex(mesh_.next_halfedge(h_ext))];

                        triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h_ext)).idx(),-1.0/12.0);
                        triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h_reg)).idx(),-1.0/12.0);

                    }else if(mesh_.valence(v_ext)==4){
                        p+=3.0/8.0*points_[v_reg] ;
                        triplets.emplace_back(nv+e.idx(),v_reg.idx(),3.0/8.0);

                        auto h1o0 = mesh_.opposite_halfedge(mesh_.next_halfedge(h_ext));
                        p -=1.0/8.0*points_[mesh_.to_vertex(mesh_.next_halfedge(h1o0))];
                        triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(mesh_.next_halfedge(h1o0)).idx(),-1.0/8.0);
                    }else{
                        double j = 0.0;
                        double val = mesh_.valence(v_ext);
                        //w_j =1/4 + cos(2pij/val) + 1/2 cos(4pij/val))/val
                        p += 1.75/val*points_[v_reg];
                        triplets.emplace_back(nv+e.idx(),v_reg.idx(),1.75/val);
                        auto h_next = mesh_.opposite_halfedge(mesh_.prev_halfedge(h_reg));
                        while(h_next != h_reg){
                            j ++;
                            double wj = (0.25 + cos(2.0*M_PI*j/val) + 0.5 *cos(4.0*M_PI*j/val))/val;
                            p+=wj*points_[mesh_.to_vertex(h_next)];
                            triplets.emplace_back(nv+e.idx(),mesh_.to_vertex(h_next).idx(),wj);
                            h_next = mesh_.opposite_halfedge(mesh_.prev_halfedge(h_next));

                        }
                    }
                    epoint[e] = p;
                }
            }
        }

        // insert new vertices on edges
        for (auto e : mesh_.edges())
        {
            mesh_.insert_vertex(e, epoint[e]);
        }

        // split faces
        Halfedge h;
        for (auto f : mesh_.faces())
        {
            h = mesh_.halfedge(f);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
            h = mesh_.next_halfedge(h);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
            h = mesh_.next_halfedge(h);
            mesh_.insert_edge(h, mesh_.next_halfedge(mesh_.next_halfedge(h)));
        }

        S_l.setFromTriplets(triplets.begin(), triplets.end());
        if (i == 0)
        {
            A0 = S_l;
        }
        else
        {
            std::cout << "A0: " << A0.rows() << "x" << A0.cols() << std::endl;
            std::cout << "Sl: " << S_l.rows() << "x" << S_l.cols() << std::endl;

            Eigen::SparseMatrix<double> S_l1 = S_l * A0;
            A0 = S_l1;
        }
        // clean-up properties
        mesh_.remove_edge_property(epoint);
    }
}

void setup_bdry_cc_A0_matrix(pmp::SurfaceMesh &mesh_,
                             Eigen::SparseMatrix<double> &A0, int lvl,
                             std::vector<Eigen::SparseMatrix<double>> &Al)
{
    Eigen::SparseMatrix<double> S_l;
    Al.resize(lvl);
    for (int i = 0; i < lvl; i++)
    {
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.n_vertices() * 4 * 4);

        pmp::VertexProperty<pmp::Point> points_ =
            mesh_.vertex_property<pmp::Point>("v:point");

        // reserve memory
        size_t nv = mesh_.n_vertices();
        size_t ne = mesh_.n_edges();
        size_t nf = mesh_.n_faces();
        mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);
        S_l.resize(nv + ne + nf, nv);

        // get properties
        auto vpoint = mesh_.add_vertex_property<pmp::Point>("catmull:vpoint");
        auto epoint = mesh_.add_edge_property<pmp::Point>("catmull:epoint");
        auto fpoint = mesh_.add_face_property<pmp::Point>("catmull:fpoint");

        // compute face vertices
        for (auto f : mesh_.faces())
        {

            pmp::Point p(0.0, 0.0, 0.0);
            pmp::Scalar c(0.0);
            double val_f = mesh_.valence(f);
            for (auto v : mesh_.vertices(f))
            {
                p += points_[v];
                ++c;
                triplets.emplace_back(nv + ne + f.idx(), v.idx(),
                                      (1.0 / val_f));
            }
            p /= c;
            fpoint[f] = p;
        }

        // compute edge vertices
        for (auto e : mesh_.edges())
        {
            // boundary or feature edge?
            if (mesh_.is_boundary(e) || mesh_.is_boundary(mesh_.face(e, 0)) ||
                mesh_.is_boundary(mesh_.face(e, 1)))
            {
                epoint[e] = 0.5 * (points_[mesh_.vertex(e, 0)] +
                                   points_[mesh_.vertex(e, 1)]);

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.5);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.5);
            }

            // interior edge
            else
            {
                pmp::Point p(0.0, 0.0, 0.0);
                p += points_[mesh_.vertex(e, 0)];
                p += points_[mesh_.vertex(e, 1)];
                p += fpoint[mesh_.face(e, 0)];
                p += fpoint[mesh_.face(e, 1)];
                p *= 0.25;
                epoint[e] = p;

                auto v0 = mesh_.vertex(e, 0);
                auto v1 = mesh_.vertex(e, 1);

                auto f0 = mesh_.face(e, 0);
                auto f1 = mesh_.face(e, 1);
                triplets.emplace_back(nv + e.idx(), v0.idx(), 0.25);
                triplets.emplace_back(nv + e.idx(), v1.idx(), 0.25);

                double val_f0 = mesh_.valence(f0);
                double val_f1 = mesh_.valence(f1);

                for (auto v : mesh_.vertices(f0))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f0) * 0.25);
                }
                for (auto v : mesh_.vertices(f1))
                {
                    triplets.emplace_back(nv + e.idx(), v.idx(),
                                          (1.0 / val_f1) * 0.25);
                }
            }
        }

        // compute new positions for old vertices
        for (auto v : mesh_.vertices())
        {
            // isolated vertex?
            if (mesh_.is_isolated(v))
            {
                vpoint[v] = points_[v];

                triplets.emplace_back(v.idx(), v.idx(), 1.0);
            }

            // boundary vertex?
            else if (mesh_.is_boundary(v))
            {
                auto h1 = mesh_.halfedge(v);
                auto h0 = mesh_.prev_halfedge(h1);

                Point p = points_[v];

                if (mesh_.valence(v) == 2)
                {
                    vpoint[v] = p;
                    triplets.emplace_back(v.idx(), v.idx(), 1.0);
                }
                else if (abs(dot(points_[mesh_.to_vertex(h1)] - p,
                                 points_[mesh_.from_vertex(h0)] - p)) < 0.00001)
                {
                    vpoint[v] = p;
                    triplets.emplace_back(v.idx(), v.idx(), 1.0);
                }
                else
                {

                    auto h1 = mesh_.halfedge(v);
                    auto h0 = mesh_.prev_halfedge(h1);

                    pmp::Point p = points_[v];
                    p *= 6.0;
                    p += points_[mesh_.to_vertex(h1)];
                    p += points_[mesh_.from_vertex(h0)];
                    p *= 0.125;

                    vpoint[v] = p;

                    auto v_to = mesh_.to_vertex(h1);
                    auto v_from = mesh_.from_vertex(h0);

                    triplets.emplace_back(v.idx(), v.idx(), 6.0 * 0.125);
                    triplets.emplace_back(v.idx(), v_to.idx(), 0.125);
                    triplets.emplace_back(v.idx(), v_from.idx(), 0.125);
                }
            }

            // interior vertex
            else
            {
                // weights from SIGGRAPH paper "Subdivision Surfaces in Character Animation"

                const pmp::Scalar k = mesh_.valence(v);
                pmp::Point p(0, 0, 0);

                for (auto vv : mesh_.vertices(v))
                {
                    p += points_[vv];
                    triplets.emplace_back(v.idx(), vv.idx(), 1.0 / (k * k));
                }

                for (auto f : mesh_.faces(v))
                {
                    p += fpoint[f];
                    double f_val = mesh_.valence(f);
                    for (auto vv : mesh_.vertices(f))
                    {
                        triplets.emplace_back(v.idx(), vv.idx(),
                                              1.0 / (k * k * f_val));
                    }
                }

                p /= (k * k);

                p += ((k - 2.0) / k) * points_[v];

                triplets.emplace_back(v.idx(), v.idx(), (k - 2.0) / k);

                vpoint[v] = p;
            }
        }

        // assign new positions to old vertices
        for (auto v : mesh_.vertices())
        {
            points_[v] = vpoint[v];
        }

        // split edges
        for (auto e : mesh_.edges())
        {
            // normal edge
            mesh_.insert_vertex(e, epoint[e]);
        }

        // split faces
        for (auto f : mesh_.faces())
        {
            auto h0 = mesh_.halfedge(f);
            mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

            auto h1 = mesh_.next_halfedge(h0);
            mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }

        // clean-up properties
        mesh_.remove_vertex_property(vpoint);
        mesh_.remove_edge_property(epoint);
        mesh_.remove_face_property(fpoint);

        S_l.setFromTriplets(triplets.begin(), triplets.end());
        //        for(int i = 0 ; i< S_l.rows(); i++){
        //            double sum = S_l.row(i).sum();
        //            if(abs(sum-1.0)> 0.000001){
        //                std::cout << "Weight sum row " << i <<" :  " << sum << std::endl;
        //            }
        //        }
        Al[i] = S_l;

        if (i == 0)
        {
            A0 = S_l;
        }
        else
        {
            std::cout << "A0: " << A0.rows() << "x" << A0.cols() << std::endl;
            std::cout << "Sl: " << S_l.rows() << "x" << S_l.cols() << std::endl;

            Eigen::SparseMatrix<double> S_l1 = S_l * A0;
            A0 = S_l1;
        }
    }
}


void setup_interpolating_sec_A0_matrix(
    pmp::SurfaceMesh &mesh_, Eigen::SparseMatrix<double> &A0, int lvl,
    std::vector<Eigen::SparseMatrix<double>> &Al)
{
    Al.resize(lvl);
    for (int i = 0; i < lvl; i++)
    {

        // reserve memory
        size_t nv = mesh_.n_vertices();
        size_t ne = mesh_.n_edges();
        size_t nf = mesh_.n_faces();
        mesh_.reserve(2 * (nv + ne + nf), 2 * ne + 4 * nf, 4 * nf);

        // get properties
        auto epoint = mesh_.add_edge_property<pmp::Point>("catmull:epoint");
        auto fpoint = mesh_.add_face_property<pmp::Point>("catmull:fpoint");
        auto points_ = mesh_.vertex_property<pmp::Point>("v:point");

        Eigen::SparseMatrix<double> S_l;

        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(mesh_.n_vertices() * 4 * 4);
        S_l.resize(nv + ne + nf, nv);

        for (auto v : mesh_.vertices())
        {
            triplets.emplace_back(v.idx(), v.idx(), 1.0);
        }
        // compute face vertices
        for (auto f : mesh_.faces())
        {
            pmp::Point p(0, 0, 0);
            pmp::Scalar c(0);

            double val_f = mesh_.valence(f);
            for (auto v : mesh_.vertices(f))
            {
                p += points_[v];
                ++c;
                triplets.emplace_back(nv + ne + f.idx(), v.idx(),
                                      (1.0 / val_f));
            }
            p /= c;
            fpoint[f] = p;
        }

        // compute edge vertices
        for (auto e : mesh_.edges())
        {
            epoint[e] = 0.5 * (points_[mesh_.vertex(e, 0)] +
                               points_[mesh_.vertex(e, 1)]);
            triplets.emplace_back(nv + e.idx(), mesh_.vertex(e, 0).idx(), 0.5);
            triplets.emplace_back(nv + e.idx(), mesh_.vertex(e, 1).idx(), 0.5);
        }

        // split edges
        for (auto e : mesh_.edges())
        {
            mesh_.insert_vertex(e, epoint[e]);
        }

        // split faces
        for (auto f : mesh_.faces())
        {
            auto h0 = mesh_.halfedge(f);
            mesh_.insert_edge(h0, mesh_.next_halfedge(mesh_.next_halfedge(h0)));

            auto h1 = mesh_.next_halfedge(h0);
            mesh_.insert_vertex(mesh_.edge(h1), fpoint[f]);

            auto h = mesh_.next_halfedge(
                mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            while (h != h0)
            {
                mesh_.insert_edge(h1, h);
                h = mesh_.next_halfedge(
                    mesh_.next_halfedge(mesh_.next_halfedge(h1)));
            }
        }
        // clean-up properties
        mesh_.remove_edge_property(epoint);
        mesh_.remove_face_property(fpoint);
        //        mesh_.remove_vertex_property(points_);
        S_l.setFromTriplets(triplets.begin(), triplets.end());
        Al[i] = S_l;
        if (i == 0)
        {
            A0 = S_l;
        }
        else
        {
            Eigen::SparseMatrix<double> S_l1 = S_l * A0;
            A0 = S_l1;
        }
    }
}

void normalize_sec_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
                             const Eigen::VectorXd &h)
{
    Eigen::SparseMatrix<double> A_0;
    SurfaceMesh subdiv_mesh = mesh;
    std::vector<Eigen::SparseMatrix<double>> Al;
    setup_cc_P_matrix(subdiv_mesh, A_0, sec_laplace_lvl, Al);
    normalize_poly_gradients(subdiv_mesh, g, A_0 * h);
}

void setup_sec_divergence_operator(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &D)
{
    SurfaceMesh subdiv_mesh = mesh;
    Eigen::SparseMatrix<double> Div_subdiv, A_0;
    std::vector<Eigen::SparseMatrix<double>> Al;
    setup_cc_P_matrix(subdiv_mesh, A_0, sec_laplace_lvl, Al);
    poly_laplace_lambda_ = sec_laplace_lambda_;
    setup_poly_divergence_operator(subdiv_mesh, Div_subdiv);
    D = A_0.transpose() * Div_subdiv;
}

void solve_poisson_non_dirichlet(pmp::SurfaceMesh &mesh,int function){

    //setup Prolongation

    //setup Stiffness and Mass

    //

}
