
#include "LaplaceConstruction.h"
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "Diamond_2D.h"
#include "DisneyLaplace.h"
#include "SECLaplace.h"
#include "AQAPoly_Laplacian.h"
#include "unsupported/Eigen/SparseExtra"

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

using namespace std;

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    Disney = 5,
    SEC = 6,
    AQAPoly_Laplace = 7,
    quadratic_Triangle_Laplace = 8

};

enum InsertedPoint
{
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2,
    Triangle_Circumcenter = 3
};

//----------------------------------------------------------------------------------

void setup_sandwich_stiffness_matrix(SurfaceMesh &mesh,
                                     Eigen::SparseMatrix<double> &S,
                                     int minpoint)
{

    const int nv = mesh.n_vertices();

    Eigen::MatrixXd Si;
    Eigen::Vector3d min;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f : mesh.faces())
    {
        const int n = mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh.vertices(f))
        {
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        if (minpoint == Centroid)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else if (minpoint == AbsAreaMinimizer)
        {
            Eigen::Vector3d point;
            optimizeAbsoluteTriangleArea(poly, point);
            find_weights_for_point(poly, point, w);
        }
        else if (minpoint == Triangle_Circumcenter)
        {
            Eigen::Vector3d point;
            point = triangle_circumcenter(poly);
            barycentric_weights_triangle(poly, point, w);
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min;

        min = poly.transpose() * w;
        localCotanMatrix(poly, min, w, Si);

        int j = 0;
        int k;
        for (Vertex v : mesh.vertices(f))
        {
            k = 0;
            for (Vertex vv : mesh.vertices(f))
            {

                trip.emplace_back(vv.idx(), v.idx(), Si(k, j));
                k++;
            }
            j++;
        }
    }

    S.resize(nv, nv);
    S.setFromTriplets(trip.begin(), trip.end());
    S *= -1.0;
}
//----------------------------------------------------------------------------------

void setup_sandwich_mass_matrix(SurfaceMesh &mesh,
                                Eigen::SparseMatrix<double> &M, int minpoint)
{
    const int nv = mesh.n_vertices();
    const int nf = mesh.n_faces();

    Eigen::MatrixXd Mi;
    Eigen::Vector3d min;
    Eigen::VectorXd w;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f : mesh.faces())
    {
        const int n = mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh.vertices(f))
        {
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // setup polygon weights
        if (minpoint == Centroid)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else if (minpoint == AbsAreaMinimizer)
        {
            Eigen::Vector3d point;
            optimizeAbsoluteTriangleArea(poly, point);
            find_weights_for_point(poly, point, w);
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }

        Eigen::Vector3d min = poly.transpose() * w;
        localMassMatrix(poly, min, w, Mi);

        int j = 0;
        int k;
        for (Vertex v : mesh.vertices(f))
        {
            k = 0;
            for (Vertex vv : mesh.vertices(f))
            {

                trip.emplace_back(vv.idx(), v.idx(), Mi(k, j));
                k++;
            }
            j++;
        }
    }
    M.resize(nv, nv);
    M.setFromTriplets(trip.begin(), trip.end());
}

//----------------------------------------------------------------------------------

void setup_sandwich_gradient_matrix(SurfaceMesh &mesh,
                                    Eigen::SparseMatrix<double> &G,
                                    int minpoint)
{

    const int nv = mesh.n_vertices();
    const int nf = mesh.n_faces();

    Eigen::MatrixXd Gi;
    Eigen::Vector3d min;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;
    int nr_triangles = 0;
    int s = 0;
    for (Face f : mesh.faces())
    {
        const int n = mesh.valence(f);
        nr_triangles += n;
        poly.resize(n, 3);
        int i = 0;

        for (auto h : mesh.halfedges(f))
        {
            Vertex v = mesh.from_vertex(h);
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }
        // compute weights for the polygon
        if (minpoint == Centroid)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else if (minpoint == AbsAreaMinimizer)
        {
            Eigen::Vector3d point;
            optimizeAbsoluteTriangleArea(poly, point);
            find_weights_for_point(poly, point, w);
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min;

        min = poly.transpose() * w;
        localGradientMatrix(poly, min, w, Gi);

        // Sandwichen
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; k++)
                    Gi(3 * i + k, j) += w(j) * Gi(3 * i + k, n);

        int j = 0;
        int k;
        for (auto vv : mesh.vertices(f))
        {
            k = 0;
            for (auto h : mesh.halfedges(f))
            {
                Vertex v = mesh.from_vertex(h);
                for (int i = 0; i < 3; i++)
                {
                    trip.emplace_back(3 * s + i, v.idx(), Gi(3 * j + i, k));
                }
                k++;
            }
            j++;
            s++;
        }
    }

    G.resize(3 * nr_triangles, nv);
    G.setFromTriplets(trip.begin(), trip.end());
}
//----------------------------------------------------------------------------------

void setup_sandwich_divergence_matrix(SurfaceMesh &mesh,
                                      Eigen::SparseMatrix<double> &D,
                                      int minpoint)
{
    const int nv = mesh.n_vertices();
    const int nf = mesh.n_faces();

    Eigen::MatrixXd Gi, Di;
    Eigen::Vector3d min;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;
    int nr_triangles = 0;
    int s = 0;
    for (Face f : mesh.faces())
    {

        const int n = mesh.valence(f);
        nr_triangles += n;
        poly.resize(n, 3);
        int i = 0;

        for (auto h : mesh.halfedges(f))
        {
            Vertex v = mesh.from_vertex(h);
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        if (minpoint == Centroid)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else if (minpoint == AbsAreaMinimizer)
        {
            Eigen::Vector3d point;
            optimizeAbsoluteTriangleArea(poly, point);
            find_weights_for_point(poly, point, w);
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }

        Eigen::Vector3d min = poly.transpose() * w;
        localGradientMatrix(poly, min, w, Gi);
        Di = -Gi.transpose();

        //Triangle Area Diagonal matrix
        Eigen::MatrixXd Ai;
        Ai.resize(3 * n, 3 * n);
        Ai.setZero();
        for (int i = 0; i < n; ++i)
        {
            const int i1 = (i + 1) % n;

            Eigen::Vector3d p0 = poly.row(i);
            Eigen::Vector3d p1 = poly.row(i1);

            double area = 0.5 * ((p0 - min).cross(p1 - min)).norm();
            for (int k = 0; k < 3; k++)
            {
                Ai(3 * i + k, 3 * i + k) = area;
            }
        }
        Di *= Ai;

        // Sandwichen

        for (int j = 0; j < n; ++j)
        {
            for (int i = 0; i < n; ++i)
            {
                for (int k = 0; k < 3; k++)
                    Di(i, 3 * j + k) += w(i) * Di(n, 3 * j + k);
            }
        }

        int j = 0;
        int k;
        for (auto vv : mesh.vertices(f))
        {
            k = 0;
            for (auto h : mesh.halfedges(f))
            {
                Vertex v = mesh.from_vertex(h);
                for (int i = 0; i < 3; i++)
                {
                    trip.emplace_back(v.idx(), 3 * s + i, Di(k, 3 * j + i));
                }
                k++;
            }
            j++;
            s++;
        }
    }

    D.resize(nv, 3 * nr_triangles);
    D.setFromTriplets(trip.begin(), trip.end());
}
//----------------------------------------------------------------------------------

void localCotanMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min,
                      Eigen::VectorXd &w, Eigen::MatrixXd &L)
{
    const int n = (int)poly.rows();

    L.resize(n, n);
    L.setZero();

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i)
    {
        const int i1 = (i + 1) % n;

        l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
        l2[0] = (poly.row(i1) - min.transpose()).squaredNorm();
        l2[1] = (poly.row(i) - min.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.5 * sqrt(arg);
        if (area > 1e-9)
        {
            l[0] = 0.25 * (l2[1] + l2[2] - l2[0]) / area;
            l[1] = 0.25 * (l2[2] + l2[0] - l2[1]) / area;
            l[2] = 0.25 * (l2[0] + l2[1] - l2[2]) / area;

            L(i1, i1) += l[0];
            L(i, i) += l[1];
            L(i1, i) -= l[2];
            L(i, i1) -= l[2];
            L(i, i) += l[2];
            L(i1, i1) += l[2];

            ln(i1) -= l[0];
            ln(i) -= l[1];
            ln(n) += l[0] + l[1];
        }
    }

    // Sandwiching
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)

            L(i, j) += w(i) * ln(j) + w(j) * ln(i) + w(i) * w(j) * ln(n);
}

//----------------------------------------------------------------------------------

void localMassMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min,
                     Eigen::VectorXd &w, Eigen::MatrixXd &M)
{
    const int n = (int)poly.rows();

    M.resize(n, n);

    M.setZero();

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i)
    {

        const int i1 = (i + 1) % n;

        l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
        l2[0] = (poly.row(i1) - min.transpose()).squaredNorm();
        l2[1] = (poly.row(i) - min.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.25 * sqrt(arg);

        l[0] = 1.0 / 6.0 * area;
        l[1] = 1.0 / 12.0 * area;

        M(i1, i1) += 1.0 / 6.0 * area;
        M(i, i) += 1.0 / 6.0 * area;
        M(i1, i) += 1.0 / 12.0 * area;
        M(i, i1) += 1.0 / 12.0 * area;

        ln(i1) += l[1];
        ln(i) += l[1];
        ln(n) += l[0];
    }
    // Sandwiching
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            M(i, j) += w(i) * ln(j) + w(j) * ln(i) + w(i) * w(j) * ln(n);
}

//----------------------------------------------------------------------------------

void localGradientMatrix(const Eigen::MatrixXd &poly,
                         const Eigen::Vector3d &min, Eigen::VectorXd &w,
                         Eigen::MatrixXd &G)
{

    const int n = (int)poly.rows();

    G.resize(3 * n, n + 1);

    G.setZero();

    Eigen::Vector3d gradient_p, gradient_p0, gradient_p1, p, p0, p1;

    p = min;
    for (int i = 0; i < n; ++i)
    {

        const int i1 = (i + 1) % n;

        p0 = poly.row(i);
        p1 = poly.row(i1);

        gradient_p = gradient_hat_function(p, p0, p1);
        gradient_p0 = gradient_hat_function(p0, p1, p);
        gradient_p1 = gradient_hat_function(p1, p, p0);
        for (int j = 0; j < 3; j++)
        {
            G(3 * i + j, n) = gradient_p(j);
            G(3 * i + j, i) = gradient_p0(j);
            G(3 * i + j, i1) = gradient_p1(j);
        }
    }

    //    std::cout << "Local Gradient pre Sandwich: " << G << std::endl;
}
//----------------------------------------------------------------------------------

void setup_stiffness_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S,
                              int Laplace, int minpoint, int degree,
                              CoarseDimension coarseningType)
{

    if (Laplace == AlexaLaplace)
    {
        setup_poly_Laplace_matrix(mesh, S);
        S *= 0.5;
    }
    else if (Laplace == CotanLaplace)
    {
        setup_triangle_Laplace_matrix(mesh, S);
        std::cout << S << std::endl;
    }
    else if (Laplace == SandwichLaplace)
    {
        setup_sandwich_stiffness_matrix(mesh, S, minpoint);
    }
    else if (Laplace == Diamond)
    {
        FaceProperty<Point> area_points;
        FaceProperty<Eigen::VectorXd> area_weights;
        if (!mesh.has_face_property("f:point"))
        {
            area_points = mesh.add_face_property<Point>("f:point");
            area_weights = mesh.add_face_property<Eigen::VectorXd>("f:weights");
        }
        setup_face_point_properties(mesh, minpoint);
        Eigen::SparseMatrix<double> G, D, G2, D2, Gra, Div, P;
        setup_prolongation_matrix(mesh, P);
        compute_primal_points(mesh, minpoint);
        setup_diamond_gradient_divergence_intrinsic(mesh, G, D);
        Gra = G * P;
        Div = P.transpose() * D;
        S = Div * Gra;
    }
    else if (Laplace == Disney)
    {
        setup_disney_laplace_operator(mesh, S);
    }
    else if (Laplace == SEC)
    {
        Eigen::SparseMatrix<double> M;
        setup_smooth_laplace_matrices(mesh, S, M);
//        setup_sec_Laplace_matrix(mesh, S);
        //        S *= 0.5;
    }
    else if (Laplace == AQAPoly_Laplace)
    {
        setup_AQAPoly_matrices(mesh, S, Stiffness, coarseningType, minpoint,
                               degree);
        if (mesh.is_triangle_mesh() && coarseningType == Edges)
        {
            Eigen::SparseMatrix<double> S_;
            setup_triangle_FEM_stiffness_matrix(mesh, S_, false);
            std::cout << "norm stiffness matrices on triangle mesh: "
                      << (S - S_).norm() << std::endl;
        }
        else if (degree == 1)
        {
            Eigen::SparseMatrix<double> S_;
            setup_sandwich_stiffness_matrix(mesh, S_, minpoint);
            std::cout << "norm stiffness matrices  : " << (S - S_).norm()
                      << std::endl;
        }
        //        }else if(degree == 2 && coarseningType == Vertices){
        //            Eigen::SparseMatrix<double> S_;
        //            setup_AQAPoly_matrices(mesh, S_, Stiffness, coarseningType, minpoint,
        //                                   1);
        //            std::cout << "norm stiffness matrices degree Vertices vs linear vertices : " << (S - S_).norm()
        //                      << std::endl;
        //        }
    }
    else if (Laplace == quadratic_Triangle_Laplace)
    {
        setup_triangle_FEM_stiffness_matrix(mesh, S, false);
    }
    //    std::cout << S << std::endl;
}

//----------------------------------------------------------------------------------

void setup_mass_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M,
                         int Laplace, int minpoint, int degree,
                         CoarseDimension coarseningType, bool lumped)
{

    if (Laplace == AlexaLaplace)
    {
        setup_poly_mass_matrix(mesh, M);
    }
    else if (Laplace == CotanLaplace)
    {
        setup_triangle_mass_matrix(mesh, M);
    }
    else if (Laplace == SandwichLaplace)
    {
        setup_sandwich_mass_matrix(mesh, M, minpoint);
        if (lumped)
        {
            lump_matrix(M);
        }
    }
    else if (Laplace == Diamond)
    {
        Eigen::SparseMatrix<double> M_, P;
        setup_prolongation_matrix(mesh, P);
        setup_diamond_mass_matrix(mesh, M_);
        M = P.transpose() * M_ * P;
    }
    else if (Laplace == Disney)
    {
        setup_disney_mass_matrix(mesh, M);
    }
    else if (Laplace == AQAPoly_Laplace)
    {
        setup_AQAPoly_matrices(mesh, M, Mass, coarseningType, minpoint, degree);
        if (degree == 1)
        {
            Eigen::SparseMatrix<double> M_;
            setup_sandwich_mass_matrix(mesh, M_, minpoint);
            if (lumped)
            {
                lump_matrix(M);
            }
            std::cout << "norm mass matrices  : " << (M_ - M).norm()
                      << std::endl;
        }
    }
    else if (Laplace == SEC)
    {

//        setup_sandwich_mass_matrix(mesh, M, minpoint);
        Eigen::SparseMatrix<double> S;
        setup_smooth_laplace_matrices(mesh, S, M);
//        setup_sec_mass_matrix(mesh, M);
        //        lump_matrix(M);
    }
    else if (Laplace == quadratic_Triangle_Laplace)
    {
        setup_triangle_FEM_mass_matrix(mesh, M);
    }
    double area = 0.0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(M, k); it; ++it)
        {
            area += it.value();
        }
    }
    //    std::cout << "------------" <<std::endl;
    //        std::cout << M << std::endl;
    std::cout << "Surface area: " << area << std::endl;
}
//----------------------------------------------------------------------------------

void testOperator(SurfaceMesh &mesh, unsigned int laplace,
                  unsigned int min_point)
{
    Eigen::SparseMatrix<double> S, M;
    setup_stiffness_matrices(mesh, S, laplace, min_point);
    setup_mass_matrices(mesh, M, laplace, min_point);
    invert_mass_matrix(M);

    if (laplace == SandwichLaplace)
    {
        std::cout << "Sandwich Laplace: " << std::endl;
    }
    else if (laplace == AlexaLaplace)
    {
        std::cout << "Alexa l=" << poly_laplace_lambda_ << " :" << std::endl;
    }
    else if (laplace == CotanLaplace)
    {
        std::cout << "Cotan Laplace" << std::endl;
    }
    else if (laplace == Diamond)
    {
        std::cout << "Diamond Laplace" << std::endl;
    }
    else if (laplace == Disney)
    {
        std::cout << "Disney Laplace" << std::endl;
    }

    std::cout << "Negative Entries Laplace: "
              << round(((double)num_negative(S) / (double)num_nonzero(S) *
                        100.0) *
                       1000) /
                     1000
              << "%" << std::endl;

    std::cout << "Nonzero Entries Laplace "
              << round(((double)sparsity(S) / (S.cols() * S.rows()) * 100.0) *
                       1000) /
                     1000

              << "%" << std::endl;

    std::cout << "nr vertices: " << mesh.n_vertices() << std::endl;
}

//------------------------------naive prolongation----------------------------
void setup_simple_2D_prolongation(SurfaceMesh &mesh,
                                  Eigen::SparseMatrix<double> &P)
{
    int nv = mesh.n_vertices();
    int ne = mesh.n_edges();
    int nf = mesh.n_faces();
    auto v_idx = mesh.add_vertex_property<Eigen::Vector2d>("v:idx");
    auto e_idx = mesh.add_edge_property<Eigen::Vector2d>("e:idx");

    Eigen::SparseMatrix<double> Pf, Pe, Pm;
    std::vector<Eigen::Triplet<double>> trip_pf;
    std::vector<Eigen::Triplet<double>> trip_pe;
    std::vector<Eigen::Triplet<double>> trip_pm;

    for (auto v : mesh.vertices())
    {
        trip_pf.emplace_back(v.idx(), v.idx(), 1.0);
        trip_pe.emplace_back(v.idx(), v.idx(), 1.0);
        trip_pm.emplace_back(v.idx(), v.idx(), 1.0);

        v_idx[v] = Eigen::Vector2d(v.idx(), v.idx());
    }
    for (auto e : mesh.edges())
    {
        trip_pf.emplace_back(nv + e.idx(), nv + e.idx(), 1.0);
        trip_pe.emplace_back(nv + e.idx(), nv + e.idx(), 1.0);
        trip_pm.emplace_back(nv + e.idx(), nv + e.idx(), 1.0);

        Vertex v0 = mesh.vertex(e, 0);
        Vertex v1 = mesh.vertex(e, 1);
        if (v0.idx() > v1.idx())
            e_idx[e] = Eigen::Vector2d(v0.idx(), v1.idx());
        else
            Eigen::Vector2d(v1.idx(), v0.idx());
    }

    int sum_val = 0;
    int ve_idx = 0;
    for (auto f : mesh.faces())
    {
        int val = mesh.valence(f);
        sum_val += val;
        Eigen::MatrixXd poly(2 * val, 3);
        Eigen::VectorXd weights;
        int i = 0;
        for (auto hf : mesh.halfedges(f))
        {
            Vertex v0 = mesh.from_vertex(hf);

            trip_pe.emplace_back(nv + ne + nf + ve_idx, v0.idx(), 0.5);
            trip_pe.emplace_back(nv + ne + nf + ve_idx, nv + ne + f.idx(), 0.5);
            ve_idx++;

            Vertex v1 = mesh.to_vertex(hf);
            Point ve = 0.5 * (mesh.position(v0) + mesh.position(v1));
            Eigen::Vector3d p0(mesh.position(v0)[0], mesh.position(v0)[1],
                               mesh.position(v0)[2]),
                pe(ve[0], ve[1], ve[2]);
            poly.row(i) = p0;
            i++;
            poly.row(i) = pe;
            i++;
        }
        find_area_minimizer_weights(poly, weights);
        int k = 0;
        i = 0;
        for (auto hf : mesh.halfedges(f))
        {
            Vertex v0 = mesh.from_vertex(hf);
            trip_pf.emplace_back(nv + ne + f.idx(), v0.idx(), weights(i));
            i++;
            trip_pf.emplace_back(nv + ne + f.idx(), nv + mesh.edge(hf).idx(),
                                 weights(i));
            i++;

            //Mishas experiment

            if (v0.idx() == val - 1)
            {
                for(int l = 0; l< nv;l++){
                    trip_pm.emplace_back(  nv + ne+l, 2 * v0.idx()+1, 0.5*weights(k));
                    trip_pm.emplace_back(  nv + ne+l, 2 * v0.idx()-1, 0.5*weights(k+1));
                }
                trip_pm.emplace_back(2 * nv + ne, 2 * v0.idx()+1, weights(k));
                k++;
                trip_pm.emplace_back(2 * nv + ne, 2 * v0.idx()-1, weights(k));
                k++;

//                trip_pm.emplace_back(2 * nv + ne, k+2, weights(k+1));
                trip_pm.emplace_back(2 * nv + ne-1, 2 * v0.idx()+1, 0.5);
            }
            else if (v0.idx() == val - 2)
            {
                for(int l = 0; l< nv;l++){
                    trip_pm.emplace_back(  nv + ne+l, 2*v0.idx(), 0.5*weights(k));
                    trip_pm.emplace_back(  nv + ne+l, 2*v0.idx()+2, 0.5*weights(k+1));
                }
                trip_pm.emplace_back(2 * nv + ne, 2*v0.idx(), weights(k));
                k++;
                trip_pm.emplace_back(2 * nv + ne, 2*v0.idx()+2, weights(k));
                k++;
                trip_pm.emplace_back( nv+ne+v0.idx(), 2*v0.idx(), 0.5);
            }
            else
            {
                for(int l = 0; l< nv;l++){
                    trip_pm.emplace_back(  nv + ne+l, 2*v0.idx(), 0.5*weights(k));
                    trip_pm.emplace_back(  nv + ne+l, 2*v0.idx()+1, 0.5*weights(k+1));

                }
                trip_pm.emplace_back(2 * nv + ne, 2*v0.idx(), weights(k));
                k++;
                trip_pm.emplace_back(2 * nv + ne, 2*v0.idx()+1, weights(k));
                k++;
                trip_pm.emplace_back( nv+ne+v0.idx(), 2*v0.idx(), 0.5);
            }
        }
        trip_pe.emplace_back(nv + ne + f.idx(), nv + ne + f.idx(), 1.0);
    }

    Pf.resize(nv + ne + nf, nv + ne);
    Pe.resize(nv + ne + nf + sum_val, nv + ne + nf);
    Pm.resize(2*nv + ne + nf , nv + ne  );

    Pf.setFromTriplets(trip_pf.begin(), trip_pf.end());
    Pe.setFromTriplets(trip_pe.begin(), trip_pe.end());
    Pm.setFromTriplets(trip_pm.begin(), trip_pm.end());
    P = Pe * Pf;

    Eigen::saveMarket(Pm, "P.txt");

    std::cout << "--------------Pm------------" << std::endl;
    std::cout << Pm << std::endl;

    std::cout << "--------------P------------" << std::endl;
    std::cout << P << std::endl;
}