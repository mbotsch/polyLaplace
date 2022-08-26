//=============================================================================

#include "Curvature.h"
#include "LaplaceConstruction.h"
#include "DisneyLaplace.h"
#include "SECLaplace.h"
#include <pmp/algorithms/SurfaceNormals.h>

//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    IntrinsicDelaunay = 4,
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

void Curvature::visualize_curvature(unsigned int laplace,
                                    unsigned int min_point, bool lumped,
                                    int degree, CoarseDimension coarseningType)
{

    if (!mesh_.n_vertices())
        return;

    // properties

    auto points = mesh_.vertex_property<Point>("v:point");
    auto curvatures = mesh_.add_vertex_property<Scalar>("v:curv");

    const unsigned int nv = mesh_.n_vertices();
    const unsigned int ne = mesh_.n_edges();

    unsigned k = 0;

    Eigen::SparseMatrix<double> M, S;
    Eigen::MatrixXd B(nv, 3);
    Eigen::VectorXd H(nv), test(3);
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        B.resize(nv + ne, 3);
        H.resize(nv + ne);
    }

    for (auto v : mesh_.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        for (auto e : mesh_.edges())
        {
            pmp::Point e_point = 0.5 * (points[(mesh_.vertex(e, 0))] +
                                        points[(mesh_.vertex(e, 1))]);
            e_point.normalize();
            B(nv + e.idx(), 0) = e_point[0];
            B(nv + e.idx(), 1) = e_point[1];
            B(nv + e.idx(), 2) = e_point[2];
        }
    }
    setup_stiffness_matrices(mesh_, S, laplace, min_point, degree,
                             coarseningType);
    setup_mass_matrices(mesh_, M, laplace, min_point, degree, coarseningType,
                        lumped);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::MatrixXd X = solver.solve(S * B);
    B = X;
    //     compute mean curvature
    for (unsigned int i = 0; i < nv; i++)
    {
        H(i) = B.row(i).norm();
    }
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        for (unsigned int i = 0; i < ne; i++)
        {
            H(nv + i) = B.row(nv + i).norm();
        }
    }
    double rms = 0.0;
    for (auto v : mesh_.vertices())
    {
        curvatures[v] = fabs(0.5 * H(k));
        if (compare_to_sphere)
        {
            double c = 0.5 * H(k);
//            std::cout << c << std::endl;
            rms += (c - 1.0) * (c - 1.0);
        }
        k++;
    }
    if ((laplace == AQAPoly_Laplace && degree == 2 &&
         coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    //    if (laplace == AQAPoly_Laplace && coarseningType == Edges)
    {
        if (compare_to_sphere)
        {
            for (unsigned int i = 0; i < ne; i++)
            {
                double c = 0.5 * H(k);
//                                std::cout << c << std::endl;
                rms += (c - 1.0) * (c - 1.0);
                k++;
            }
        }
    }

    if (compare_to_sphere)
    {
        if ((laplace == AQAPoly_Laplace && coarseningType == Vertices) ||
            (laplace != quadratic_Triangle_Laplace && degree != 2))
        {
            std::cout << "linear Dof" << std::endl;
            rms /= (double)nv;
        }
        else
        {
            std::cout << "quadratic Dof" << std::endl;
            rms /= (double)(nv + ne);
        }

        rms = sqrt(rms);

        if (laplace == AlexaLaplace)
        {
            std::cout << "Curvature deviation (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace == SEC)
        {
            std::cout << "Curvature deviation (SEC, l=" << sec_laplace_lambda_
                      << ", lvl= " << sec_laplace_lvl << "): " << rms
                      << std::endl;
        }
        else if (laplace == IntrinsicDelaunay)
        {
            std::cout << "Curvature deviation intrinsic Delaunay: " << rms
                      << std::endl;
        }
        else if (laplace == CotanLaplace)
        {
            std::cout << "Curvature deviation Cotan Laplacian: " << rms
                      << std::endl;
        }
        else if (laplace == Disney)
        {
            std::cout << "Curvature deviation (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace == quadratic_Triangle_Laplace)
        {
            std::cout << "Curvature deviation quadratic Triangle Laplacian: "
                      << rms << std::endl;
        }
        else
        {
            if (laplace == Diamond)
            {
                std::cout << "Diamond Laplace ";
            }
            else if (laplace == AQAPoly_Laplace)
            {
                std::cout << "AQAPoly ";
                if (degree == 2)
                {
                    std::cout << "quadratic: ";
                }
                else
                {
                    std::cout << "linear: ";
                }
            }
            else
            {
                std::cout << "Sandwich Laplace ";
            }

            if (min_point == AbsAreaMinimizer)
            {
                std::cout << "Curvature deviation (abs area): " << rms
                          << std::endl;
            }
            else if (min_point == AreaMinimizer)
            {
                std::cout << "Curvature deviation (our Point): " << rms
                          << std::endl;
            }
            else if (min_point == Triangle_Circumcenter)
            {
                std::cout << "Curvature deviation (Triangle Circumcenter): "
                          << rms << std::endl;
            }
            else
            {

                std::cout << "Curvature deviation (centroid): " << rms
                          << std::endl;
            }
        }
    }

    curvature_to_texture_coordinates();
    mesh_.remove_vertex_property<Scalar>(curvatures);
}

//-----------------------------------------------------------------------------

double Curvature::compute_curvature_error(unsigned int laplace,
                                          unsigned int min_point, bool lumped,
                                          int degree,
                                          CoarseDimension coarseningType)
{
    if (!mesh_.n_vertices())
        return -10000;

    auto points = mesh_.vertex_property<Point>("v:point");
    const unsigned int nv = mesh_.n_vertices();
    const unsigned int ne = mesh_.n_edges();

    unsigned k = 0;

    Eigen::SparseMatrix<double> M, S;
    Eigen::MatrixXd B(nv, 3);
    Eigen::VectorXd H(nv), test(3);
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        B.resize(nv + ne, 3);
        H.resize(nv + ne);
    }

    for (auto v : mesh_.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        for (auto e : mesh_.edges())
        {
            pmp::Point e_point = 0.5 * (points[(mesh_.vertex(e, 0))] +
                                        points[(mesh_.vertex(e, 1))]);
            e_point.normalize();
            B(nv + e.idx(), 0) = e_point[0];
            B(nv + e.idx(), 1) = e_point[1];
            B(nv + e.idx(), 2) = e_point[2];
        }
    }
    setup_stiffness_matrices(mesh_, S, laplace, min_point, degree,
                             coarseningType);
    setup_mass_matrices(mesh_, M, laplace, min_point, degree, coarseningType,
                        lumped);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::MatrixXd X = solver.solve(S * B);
    B = X;

    //     compute mean curvature
    for (unsigned int i = 0; i < nv; i++)
    {
        H(i) = B.row(i).norm();
    }
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        for (unsigned int i = 0; i < ne; i++)
        {
            H(nv + i) = B.row(nv + i).norm();
        }
    }
    double rms = 0.0;
    for (auto v : mesh_.vertices())
    {
        double c = 0.5 * H(k);
        rms += (c - 1.0) * (c - 1.0);
//        std::cout << c << std::endl;

        k++;
    }
    if ((laplace == AQAPoly_Laplace && coarseningType == Edges) ||
        laplace == quadratic_Triangle_Laplace)
    {
        for (unsigned int i = 0; i < ne; i++)
        {
            double c = 0.5 * H(k);
            rms += (c - 1.0) * (c - 1.0);
            k++;
        }
    }

    if ((laplace == AQAPoly_Laplace &&
         (degree == 1 || coarseningType == Vertices)) ||
        laplace != quadratic_Triangle_Laplace)
    {
        std::cout << "linear Dof" << std::endl;
        rms /= (double)nv;
    }
    else
    {
        std::cout << "quadratic Dof" << std::endl;
        rms /= (double)(nv + ne);
    }

    rms = sqrt(rms);

    if (laplace == AlexaLaplace)
    {
        std::cout << "Curvature deviation (Alexa, l=" << poly_laplace_lambda_
                  << "): " << rms << std::endl;
    }
    else if (laplace == SEC)
    {
        std::cout << "Curvature deviation (SEC, l=" << sec_laplace_lambda_
                  << ", lvl= " << sec_laplace_lvl << "): " << rms << std::endl;
    }
    else if (laplace == IntrinsicDelaunay)
    {
        std::cout << "Curvature deviation intrinsic Delaunay: " << rms
                  << std::endl;
    }
    else if (laplace == CotanLaplace)
    {
        std::cout << "Curvature deviation Cotan Laplacian: " << rms
                  << std::endl;
    }
    else if (laplace == Disney)
    {
        std::cout << "Curvature deviation (Disney, l=" << disney_laplace_lambda_
                  << "): " << rms << std::endl;
    }
    else
    {
        if (laplace == Diamond)
        {
            std::cout << "Diamond Laplace ";
        }
        else if (laplace == AQAPoly_Laplace)
        {
            std::cout << "AQAPoly ";
            if (degree == 2)
            {
                std::cout << "quadratic: ";
            }
            else
            {
                std::cout << "linear: ";
            }
        }
        else
        {
            std::cout << "Sandwich Laplace ";
        }

        if (min_point == AbsAreaMinimizer)
        {
            std::cout << "Curvature deviation (abs area): " << rms << std::endl;
        }
        else if (min_point == AreaMinimizer)
        {
            std::cout << "Curvature deviation (our Point): " << rms
                      << std::endl;
        }
        else if (min_point == Triangle_Circumcenter)
        {
            std::cout << "Curvature deviation (Triangle Circumcenter): " << rms
                      << std::endl;
        }
        else
        {

            std::cout << "Curvature deviation (centroid): " << rms << std::endl;
        }
    }

    return rms;
}
double Curvature::compute_quad_Tri_Laplace_curvature_error(
    SurfaceMesh& initial_mesh)
{
    auto points = mesh_.vertex_property<Point>("v:point");
    const unsigned int nv = mesh_.n_vertices();
    const unsigned int ne = mesh_.n_edges();

    const unsigned int nv_orig = initial_mesh.n_vertices();
    const unsigned int ne_orig = initial_mesh.n_edges();

    Eigen::SparseMatrix<double> M, S;
    Eigen::MatrixXd B(nv + ne, 3);
    Eigen::VectorXd H(nv + ne);

    for (auto v : mesh_.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }
    for (auto e : mesh_.edges())
    {
        pmp::Point e_point =
            0.5 * (points[(mesh_.vertex(e, 0))] + points[(mesh_.vertex(e, 1))]);
        e_point.normalize();
        B(nv + e.idx(), 0) = e_point[0];
        B(nv + e.idx(), 1) = e_point[1];
        B(nv + e.idx(), 2) = e_point[2];
    }

    setup_stiffness_matrices(mesh_, S, quadratic_Triangle_Laplace);
    setup_mass_matrices(mesh_, M, quadratic_Triangle_Laplace);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    Eigen::MatrixXd X = solver.solve(S * B);
    B = X;

    //     compute mean curvature
    for (unsigned int i = 0; i < nv; i++)
    {
        H(i) = B.row(i).norm();
    }
    for (unsigned int i = 0; i < ne; i++)
    {
        H(nv + i) = B.row(nv + i).norm();
    }

    double rms = 0.0;
    for (auto v : initial_mesh.vertices())
    {
        double c = 0.5 * H(v.idx());
        std::cout << c << std::endl;
        rms += (c - 1.0) * (c - 1.0);
    }

    for (auto e : initial_mesh.edges())
    {
        double c = 0.5 * H(nv+e.idx());
        rms += (c - 1.0) * (c - 1.0);
    }

    rms /= (double)(nv_orig + ne_orig);
    rms = sqrt(rms);

    std::cout << "quad Tri Laplace curvature: " << rms << std::endl;
    return rms;
}

void Curvature::curvature_to_texture_coordinates() const
{
    auto curvatures = mesh_.get_vertex_property<Scalar>("v:curv");
    assert(curvatures);

    // sort curvature values
    std::vector<Scalar> values;
    values.reserve(mesh_.n_vertices());

    for (auto v : mesh_.vertices())
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
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    if (kmin < 0.0) // signed
    {
        kmax = std::max(fabs(kmin), fabs(kmax));
        for (auto v : mesh_.vertices())
        {
            tex[v] = TexCoord((0.5f * curvatures[v] / kmax) + 0.5f, 0.0);
        }
    }
    else // unsigned
    {
        for (auto v : mesh_.vertices())
        {
            tex[v] = TexCoord((curvatures[v] - kmin) / (kmax - kmin), 0.0);
        }
    }

    // remove per-halfedge texture coordinates
    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh_.remove_halfedge_property(htex);
}

//=============================================================================
