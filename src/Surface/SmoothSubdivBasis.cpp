#include "pmp/algorithms/SurfaceSubdivision.h"
#include "SmoothSubdivBasis.h"
#include "DisneyLaplace.h"
#include "LaplaceConstruction.h"
#include "Diamond_2D.h"
#include "diffgeo.h"
#include "Poisson_System.h"
#include "igl/slice.h"

//=============================================================================

using SurfaceMesh = pmp::SurfaceMesh;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

int subdiv_lvl = 1;
//=============================================================================
enum Basis
{
    CC_Subdiv_CC_P = 0,
    Lin_Subdiv_CC_P = 1,
    Lin_Subdiv_Lin_P = 2,
    SEC = 3
};

enum Function
{
    LinearPrecision = 0,
    Franke = 1,
    SphericalHarmonics = 2
};
void setup_subdiv_stiffness_matrix(pmp::SurfaceMesh &mesh,
                                   Eigen::SparseMatrix<double> &S, int kind)
{
    SurfaceMesh basis_mesh = mesh;
    SurfaceMesh lin_interpolation_mesh = mesh;

    Eigen::SparseMatrix<double> S_subdiv, P;
    for (int i = 0; i < subdiv_lvl; i++)
    {
        linear_interpolation_catmull_clark(lin_interpolation_mesh);
    }
    std::cout << subdiv_lvl << std::endl;
    setup_cc_prolongation(basis_mesh, P, subdiv_lvl);
    setup_sandwich_stiffness_matrix(lin_interpolation_mesh, S_subdiv, 0);
    S = P.transpose() * S_subdiv * P;
}

void setup_subdiv_mass_matrix(pmp::SurfaceMesh &mesh,
                              Eigen::SparseMatrix<double> &M, int kind)
{

    SurfaceMesh basis_mesh = mesh;
    SurfaceMesh lin_interpolation_mesh = mesh;

    Eigen::SparseMatrix<double> M_subdiv, P;
    for (int i = 0; i < subdiv_lvl; i++)
    {
        linear_interpolation_catmull_clark(lin_interpolation_mesh);
    }
    std::cout << subdiv_lvl << std::endl;
    setup_cc_prolongation(basis_mesh, P, subdiv_lvl);
    setup_sandwich_mass_matrix(lin_interpolation_mesh, M_subdiv, 0);
    M = P.transpose() * M_subdiv * P;
}

void setup_subdiv_system_matrices(pmp::SurfaceMesh &mesh,
                                  Eigen::SparseMatrix<double> &S,
                                  Eigen::SparseMatrix<double> &M,
                                  Eigen::SparseMatrix<double> &P, int kind)
{
    SurfaceMesh basis_mesh = mesh;
    SurfaceMesh lin_interpolation_mesh = mesh;
    SurfaceMesh cc_mesh = mesh;

    Eigen::SparseMatrix<double> M_subdiv, S_subdiv, Pt;
    SurfaceSubdivision divider = SurfaceSubdivision(cc_mesh);
    Eigen::MatrixXd Pos(mesh.n_vertices(), 3), PosCC;
    for (auto v : mesh.vertices())
    {
        for (int i = 0; i < 3; i++)
        {
            Pos(v.idx(), i) = mesh.position(v)[i];
        }
    }
    for (int i = 0; i < subdiv_lvl; i++)
    {
        linear_interpolation_catmull_clark(lin_interpolation_mesh);
        divider.catmull_clark();
    }

    if (kind == CC_Subdiv_CC_P)
    {
        setup_cc_prolongation(basis_mesh, P, subdiv_lvl);
        std::cout << "CC Subdiv, CC Prolongation" << std::endl;
        setup_sandwich_mass_matrix(cc_mesh, M_subdiv, 0);
        setup_sandwich_stiffness_matrix(cc_mesh, S_subdiv, 0);
    }
    else if (kind == Lin_Subdiv_CC_P)
    {
        std::vector<Eigen::SparseMatrix<double>> Al;
        setup_cc_prolongation(basis_mesh, P, subdiv_lvl);
        std::cout << "Linear Subdiv, CC Prolongation" << std::endl;
        setup_sandwich_mass_matrix(lin_interpolation_mesh, M_subdiv, 0);
        setup_sandwich_stiffness_matrix(lin_interpolation_mesh, S_subdiv, 0);
    }

    if (kind == Lin_Subdiv_Lin_P)
    {
        setup_interpolating_prolongation(basis_mesh, P, subdiv_lvl);
        setup_sandwich_mass_matrix(lin_interpolation_mesh, M_subdiv, 0);
        setup_sandwich_stiffness_matrix(lin_interpolation_mesh, S_subdiv, 0);
        std::cout << "Linear Subdiv, Linear Prolongation" << std::endl;
    }
    if (kind == SEC)
    {
        setup_cc_prolongation(basis_mesh, P, subdiv_lvl);
        setup_poly_mass_matrix(lin_interpolation_mesh, M_subdiv);
        setup_poly_Laplace_matrix(lin_interpolation_mesh, S_subdiv);
        S_subdiv *= 0.5;

        std::cout << "SEC" << std::endl;
    }

    Pt = P.transpose();
    S = Pt * S_subdiv * P;
    M = Pt * M_subdiv * P;

    double area = 0.0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(M, k); it; ++it)
        {
            area += it.value();
        }
    }
    std::cout << "Surface area: " << area << std::endl;
}

double solve_poisson_non_dirichlet(pmp::SurfaceMesh &mesh, int function,
                                   int kind)
{
    Eigen::SparseMatrix<double> S, M, P;
    setup_subdiv_system_matrices(mesh, S, M, P, kind);
    int nv = mesh.n_vertices();
    int l=2, m=4;
    Eigen::VectorXi orig(mesh.n_vertices()), bdry_v, in;
    int nb = 0;
    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
            nb++;
        orig(v.idx()) = v.idx();
    }

    Eigen::MatrixXd X, B;
    Eigen::VectorXd b;
    if (function == LinearPrecision)
    {
        B.resize(nv + nb, 2);
    }
    else if (function == Franke)
    {
        B.resize(nv + nb, 1);
        b.resize(nv);
    }else if (function == SphericalHarmonics){
        B.resize(nv, 1);
        b.resize(nv);
    }
    else
    {
//        int k = nv;
//        for (auto v : mesh.vertices())
//        {
//            analytic_solution(v.idx()) =
//                sphericalHarmonic(mesh.position(v), l, m);
//        }
//        if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
//            laplace == quadratic_Triangle_Laplace)
//        {
//            k = nv + ne;
//            for (auto e : mesh.edges())
//            {
//                pmp::Point p0 = mesh.position(mesh.vertex(e, 0));
//                Point p1 = mesh.position(mesh.vertex(e, 1));
//                Point p = 0.5 * (p0 + p1);
//                p.normalize();
//                analytic_solution(nv + e.idx()) = sphericalHarmonic(p, l, m);
//            }
//        }
        solver.analyzePattern(M);
        solver.factorize(M);
        analytic_solution.normalize();
        X = solver.solve(S * analytic_solution);
        if (solver.info() != Eigen::Success)
        {
            std::cout << "Issue: " << solver.info() << std::endl;
            std::cerr << "Could not solve linear system\n";
        }
        //        std::cout << "sum M*X: " << (M * X).sum() << std::endl;
        //        std::cout << "sum X : " << X.sum() << std::endl;
        double eval = -l * (l + 1);
        error = (analytic_solution - 1.0 / eval * X).transpose() * M *
                (analytic_solution - 1.0 / eval * X);
        error = sqrt(error / double(k));
        std::cout << "error inner product : " << error << std::endl;
        return error;
    }

    B.setZero();
    double error = 0.0;

    bdry_v.resize(nb);
    in.resize(nv - nb);

    int c_bdry = 0;
    int c_in = 0;
    // Set the constraints at the locked vertices to the evluation of the Franke function
    for (auto v : mesh.vertices())
    {
        if (mesh.is_boundary(v))
        {
            // right-hand side: fix boundary values with franke function of the vertices
            if (function == LinearPrecision)
            {

                for (int i = 0; i < 2; i++)
                {
                    B(nv + c_bdry, i) = mesh.position(v)[i];
                    //                      B(nv + c_bdry, i) = 1.0;
                }
            }
            else if (function == Franke)
            {

                B(nv + c_bdry, 0) =
                    franke_function(mesh.position(v)[0], mesh.position(v)[1]);
            }
            bdry_v(c_bdry) = v.idx();
            c_bdry++;
        }
        else
        {
            if (c_in == 0)
            {
                bdry_v(c_bdry) = v.idx();
            }
            in(c_in) = v.idx();
            c_in++;
        }
        if (function == Franke)
        {
            b(v.idx()) = -laplace_franke_function(mesh.position(v)[0],
                                                  mesh.position(v)[1]);
        }else if (function == SphericalHarmonics){
            b(v.idx()) = sphericalHarmonic(mesh.position(v), l, m);
        }
    }

    std::vector<Eigen::Triplet<double>> tripletsA;
    tripletsA.reserve((nv + nb) * (8));

    // Slice boundary function values of boundary vertices for all original nodes out of the Prolongation matrix
    Eigen::SparseMatrix<double> P_boundary, P_basis, A_sparse;
    igl::slice(P, bdry_v, orig, P_boundary);
    igl::slice(P, orig, orig, P_basis);
    if (function == Franke)
    {
        //        lump_matrix(M);
        b = M * b;
        B.block(0, 0, nv, 1) = b;
    }
    //std::cout << P_boundary << std::endl;
    // is always one, no issue here
    //    for(int i =0; i< P_boundary.rows();i++){
    //        double row_sum = P_boundary.row(i).sum();
    //        std::cout << row_sum << std::endl;
    //    }

    for (int k = 0; k < S.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(S, k); it; ++it)
        {
            tripletsA.emplace_back(it.row(), it.col(), -it.value());
        }
    }

    for (int k = 0; k < P_boundary.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(P_boundary, k); it;
             ++it)
        {
            tripletsA.emplace_back(nv + it.row(), it.col(), it.value());
            tripletsA.emplace_back(it.col(), nv + it.row(), it.value());
        }
    }

    A_sparse.resize(nv + nb, nv + nb);
    A_sparse.setFromTriplets(tripletsA.begin(), tripletsA.end());

    double difference = 0.0;
    for (int k = 0; k < S.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(S, k); it; ++it)
        {
            difference += sqrt(
                pow(-A_sparse.coeffRef(it.row(), it.col()) - it.value(), 2.0));
        }
    }
    //    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    //    Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>
    //        solver;
    solver.compute(A_sparse);
    Eigen::MatrixXd x = solver.solve(B);
    Eigen::MatrixXd sol;
    if (function == Franke)
    {
        sol = P_basis * x.block(0, 0, nv, 1);
    }
    else if (function == LinearPrecision)
    {
        sol = P_basis * x.block(0, 0, nv, 2);
    }
    for (auto v : mesh.vertices())
    {
        Eigen::Vector2d p, p1;
        p << mesh.position(v)[0], mesh.position(v)[1];

        p1 << sol(v.idx(), 0), sol(v.idx(), 1);
        if (function == LinearPrecision)
        {


            std::cout << "Point: (" << mesh.position(v)[0] << ","
                      << mesh.position(v)[1] << ") "
                      << " computed pos: (" << sol(v.idx(), 0) << ","
                      << sol(v.idx(), 1) << ")" << std::endl;
            //            }
            error += pow((p1 - p).norm(), 2.);
        }
        else if (function == Franke)
        {
            //            if (!mesh.is_boundary(v))
            //            {

            //                std::cout << "Franke:"
            //                          << franke_function(mesh.position(v)[0],
            //                                             mesh.position(v)[1])
            //                          << " sol: " << sol(v.idx(), 0) << std::endl;
            //            }
            error +=
                pow(franke_function(mesh.position(v)[0], mesh.position(v)[1]) -
                        sol(v.idx(), 0),
                    2.);
        }
    }
    std::cout << "Difference between A_sparse and S: " << difference
              << std::endl;
    std::cout << "Mesh vertices: " << mesh.n_vertices()
              << " sol dim: " << sol.rows() << std::endl;
    //    std::cout << "Rank System: "<< solver.rank() << " Dimension: " << A_sparse.rows() << "x" << A_sparse.cols()<<std::endl;
    std::cout << "Subdiv Lvl: " << subdiv_lvl << std::endl;
    std::cout << " RMSE:" << sqrt(error / (double)mesh.n_vertices())
              << std::endl;

    //====================================================================
    //Visualize Results
    std::vector<Scalar> values;
    values.reserve(mesh.n_vertices());
    // generate 1D texture coordiantes
    auto tex = mesh.vertex_property<TexCoord>("v:tex");
    auto correct_tex = mesh.vertex_property<TexCoord>("v:correct_tex");

    double sum = 0.0;
    double neg_offset = 0.0;
    if (function == LinearPrecision)
    {
        for (auto v : mesh.vertices())
        {
            if (sol(v.idx(), 0) < neg_offset)
            {
                neg_offset = sol(v.idx(), 0);
            }
        }
        for (auto v : mesh.vertices())
        {
            tex[v] = TexCoord(sol(v.idx(), 0) - neg_offset, neg_offset);
        }
    }
    else if (function == Franke)
    {
        for (auto v : mesh.vertices())
        {
            tex[v] = TexCoord(sol(v.idx(), 0), 0.0);
        }
    }
    // remove per-halfedge texture coordinates
    auto htex = mesh.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh.remove_halfedge_property(htex);

    //====================================================================
    return sqrt(error / (double)mesh.n_vertices());
}
double curvature_non_dirichlet(pmp::SurfaceMesh &mesh, int kind)
{

    // properties
    auto points = mesh.vertex_property<Point>("v:point");
    auto curvatures = mesh.add_vertex_property<Scalar>("v:curv");

    const unsigned int nv = mesh.n_vertices();

    unsigned k = 0;

    Eigen::SparseMatrix<double> M, S, P, P_basis;
    Eigen::MatrixXd B(nv, 3);
    Eigen::VectorXd H(nv), test(3);
    Eigen::VectorXi orig(mesh.n_vertices());
    for (auto v : mesh.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
        orig(v.idx()) = v.idx();
    }

    setup_subdiv_system_matrices(mesh, S, M, P, kind);
    igl::slice(P, orig, orig, P_basis);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    //        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;

    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::MatrixXd X = solver.solve(S * B);
    B = P_basis * X;
    //     compute mean curvature
    for (unsigned int i = 0; i < nv; i++)
    {
        H(i) = B.row(i).norm();
    }
    double rms = 0.0;
    for (auto v : mesh.vertices())
    {
        curvatures[v] = fabs(0.5 * H(v.idx()));
        double c = 0.5 * H(v.idx());

                    std::cout << c << std::endl;
        rms += (c - 1.0) * (c - 1.0);
    }
    rms /= (double)nv;
    rms = sqrt(rms);

    std::cout << "RMSE unit sphere: " << rms << std::endl;

    //Curvature to texture
    // sort curvature values
    std::vector<Scalar> values;
    values.reserve(mesh.n_vertices());

    for (auto v : mesh.vertices())
    {
        values.push_back(curvatures[v]);
    }
    std::sort(values.begin(), values.end());
    unsigned int n = values.size() - 1;

    // clamp upper/lower 5%
    unsigned int i = n / 20;
    Scalar kmin = values[i];
    Scalar kmax = values[n - 1 - i];

    // generate 1D texture coordiantes
    auto tex = mesh.vertex_property<TexCoord>("v:tex");
    if (kmin < 0.0) // signed
    {
        kmax = std::max(fabs(kmin), fabs(kmax));
        for (auto v : mesh.vertices())
        {
            tex[v] = TexCoord((0.5f * curvatures[v] / kmax) + 0.5f, 0.0);
        }
    }
    else // unsigned
    {
        for (auto v : mesh.vertices())
        {
            tex[v] = TexCoord((curvatures[v] - kmin) / (kmax - kmin), 0.0);
        }
    }

    // remove per-halfedge texture coordinates
    auto htex = mesh.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh.remove_halfedge_property(htex);
    mesh.remove_vertex_property<Scalar>(curvatures);
}

void setup_cc_prolongation(pmp::SurfaceMesh &mesh_,
                           Eigen::SparseMatrix<double> &P, int lvl)
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
        if (i == 0)
        {
            P = S_l;
        }
        else
        {
            std::cout << "P: " << P.rows() << "x" << P.cols() << std::endl;
            std::cout << "Sl: " << S_l.rows() << "x" << S_l.cols() << std::endl;

            Eigen::SparseMatrix<double> S_l1 = S_l * P;
            P = S_l1;
        }
    }
}

void linear_interpolation_catmull_clark(pmp::SurfaceMesh &mesh_)
{
    // reserve memory
    size_t nv = mesh_.n_vertices();
    size_t ne = mesh_.n_edges();
    size_t nf = mesh_.n_faces();
    mesh_.reserve(nv + ne + nf, 2 * ne + 4 * nf, 4 * nf);

    // get properties
    auto epoint = mesh_.add_edge_property<pmp::Point>("lin:epoint");
    auto fpoint = mesh_.add_face_property<pmp::Point>("lin:fpoint");
    auto points_ = mesh_.vertex_property<pmp::Point>("v:point");

    // compute face vertices
    for (auto f : mesh_.faces())
    {
        pmp::Point p(0, 0, 0);
        double c(0);
        for (auto v : mesh_.vertices(f))
        {
            p += points_[v];
            ++c;
        }
        p /= c;
        fpoint[f] = p;
    }

    // compute edge vertices
    for (auto e : mesh_.edges())
    {
        epoint[e] =
            0.5f * (points_[mesh_.vertex(e, 0)] + points_[mesh_.vertex(e, 1)]);
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

        auto h =
            mesh_.next_halfedge(mesh_.next_halfedge(mesh_.next_halfedge(h1)));
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
}

void setup_interpolating_prolongation(pmp::SurfaceMesh &mesh_,
                                      Eigen::SparseMatrix<double> &P, int lvl)
{
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
        if (i == 0)
        {
            P = S_l;
        }
        else
        {
            Eigen::SparseMatrix<double> S_l1 = S_l * P;
            P = S_l1;
        }
    }
}

void setup_smooth_basis_matrix(pmp::SurfaceMesh &mesh, Eigen::MatrixXd &B,
                               int kind)
{
    const unsigned int nv = mesh.n_vertices();
    Eigen::SparseMatrix<double> S_subdiv, M_subdiv;
    SurfaceMesh linear_interpolation_mesh = mesh;
    SurfaceMesh catmull_clark_mesh = mesh;
    Eigen::SparseMatrix<double> P, P_cc, P_init, P_lin, L_subdiv;

    if (kind == CC_Subdiv_CC_P || kind == Lin_Subdiv_CC_P)
    {
        setup_cc_prolongation(catmull_clark_mesh, P, subdiv_lvl);
        std::cout << "Catmull Clark Prolongation" << std::endl;
    }

    else if (kind == Lin_Subdiv_Lin_P)
    {
        setup_interpolating_prolongation(catmull_clark_mesh, P, subdiv_lvl);
        std::cout << "interpolating Prolongation" << std::endl;
    }
    Eigen::MatrixXd E_i(nv, nv);
    //    Eigen::MatrixXd E_ii(catmull_clark_mesh.n_vertices(), catmull_clark_mesh.n_vertices());

    E_i.setIdentity();
    //    E_ii.setIdentity();

    B = P * E_i;
    //    std::cout << B << std::endl;
    //    std::cout << "-----------------" <<std::endl;
    //    std::cout << P.transpose()*E_ii << std::endl;
}
