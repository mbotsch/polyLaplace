//=============================================================================

#include "Curvature.h"


//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

enum LaplaceMethods
{
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    deGoesLaplace = 4
};

enum InsertedPoint
{
    Centroid = 0,
    AreaMinimizer = 2
};

void Curvature::visualize_curvature(unsigned int laplace, unsigned int min_point,bool lumped)
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

    for (auto v : mesh_.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }

    setup_stiffness_matrices(mesh_, S, laplace, min_point);
    setup_mass_matrices(mesh_, M, laplace, min_point, lumped);
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

    if (compare_to_sphere)
    {

        rms /= (double)nv;
        rms = sqrt(rms);

        if (laplace == AlexaWardetzkyLaplace)
        {
            std::cout << "Curvature deviation (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace == CotanLaplace)
        {
            std::cout << "Curvature deviation Cotan Laplacian: " << rms
                      << std::endl;
        }
        else if (laplace == deGoesLaplace)
        {
            std::cout << "Curvature deviation (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else
        {
            if (laplace == Diamond)
            {
                std::cout << "Diamond Laplace ";
            }
            else
            {
                std::cout << "Sandwich Laplace ";
            }
           if (min_point == AreaMinimizer)
            {
                std::cout << "Curvature deviation (our Point): " << rms
                          << std::endl;
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
                                          unsigned int min_point, bool lumped)
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

    for (auto v : mesh_.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }


    setup_stiffness_matrices(mesh_, S, laplace, min_point);
    setup_mass_matrices(mesh_, M, laplace, min_point,lumped);

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
    double rms = 0.0;
    for (auto v : mesh_.vertices())
    {
        double c = 0.5 * H(k);
        rms += (c - 1.0) * (c - 1.0);
//        std::cout << c << std::endl;

        k++;
    }

        rms /= (double)nv;

    rms = sqrt(rms);

    if (laplace == AlexaWardetzkyLaplace)
    {
        std::cout << "Curvature deviation (Alexa, l=" << poly_laplace_lambda_
                  << "): " << rms << std::endl;
    }
    else if (laplace == CotanLaplace)
    {
        std::cout << "Curvature deviation Cotan Laplacian: " << rms
                  << std::endl;
    }
    else if (laplace == deGoesLaplace)
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

        else
        {
            std::cout << "Sandwich Laplace ";
        }

      if (min_point == AreaMinimizer)
        {
            std::cout << "Curvature deviation (our Point): " << rms
                      << std::endl;
        }
        else
        {

            std::cout << "Curvature deviation (centroid): " << rms << std::endl;
        }
    }

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
