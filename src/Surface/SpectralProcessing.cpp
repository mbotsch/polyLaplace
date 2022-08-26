
//=============================================================================

#include "SpectralProcessing.h"
#include "LaplaceConstruction.h"
#include "[dGBD20]Laplace.h"
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/SymGEigsSolver.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/Util/GEigsMode.h>
#include <iostream>
#include <fstream>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/Util/CompInfo.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include "HarmonicBasis2D.h"

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
    AreaMinimizer = 2,
};

double solve_eigenvalue_problem(SurfaceMesh &mesh, int laplace, int face_point,
                                int degree, CoarseDimension coarseningType,
                                std::string meshname)
{
    std::string filename;
    if (laplace == Diamond)
    {

       filename = "eigenvalues_[BBA21]_" + meshname + ".csv";

    }
    else if (laplace == AlexaWardetzkyLaplace)
    {
        filename = "eigenvalues_[AW11]_" + meshname + ".csv";
    }
    else if (laplace == deGoesLaplace)
    {
        filename = "eigenvalues_[dGBD20]_" + meshname + ".csv";
    }
    else if (laplace == PolySimpleLaplace)
    {
        filename = "eigenvalues_[BHKB20]_" + meshname + ".csv";
    }
    std::ofstream ev_file(filename);

    ev_file << "computed,analytic,offset" << std::endl;
    Eigen::SparseMatrix<double> M, S;

    setup_stiffness_matrices(mesh, S, laplace, face_point);
    setup_mass_matrices(mesh, M, laplace, face_point);

    // Construct matrix operation object using the wrapper class SparseGenMatProd

    int num_eval = 49;
    int converge_speed = 5 * num_eval;

    // S and M are sparse
    using OpType =
        Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    OpType op(S, M);
    BOpType Bop(M);

    // Construct generalized eigen solver object, seeking three generalized
    // eigenvalues that are closest to zero. This is equivalent to specifying
    // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule

    Spectra::SymGEigsShiftSolver<OpType, BOpType,
                                 Spectra::GEigsMode::ShiftInvert>
        geigs(op, Bop, num_eval, converge_speed, 1e-8);
    geigs.init();
    geigs.compute(Spectra::SortRule::LargestMagn);
    std::cout << "compute & init" << std::endl;

    Eigen::VectorXd analytic, evalues;
    Eigen::MatrixXd evecs;
    // Retrieve results
    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        evalues.resize(num_eval);
        analytic.resize(num_eval);
        evalues = geigs.eigenvalues();
    }
    else
    {
        std::cout << "Eigenvalue computation failed!\n" << std::endl;
    }

    analytic_eigenvalues_unitsphere(analytic, num_eval);
    double error = 0.0;

    for (int i = 1; i < evalues.size(); i++)
    {
        std::cout << evalues(i) << "," << analytic(i) << ","
                  << evalues(i) - analytic(i) << std::endl;
        ev_file << evalues(i) << "," << analytic(i) << ","
                << evalues(i) - analytic(i) << std::endl;
        error += pow(evalues(i) - analytic(i), 2);
    }

    error = sqrt(error / (double)evalues.size());
    std::cout << "Root mean squared error: " << error << std::endl;
    ev_file.close();
    return error;
}

void analytic_eigenvalues_unitsphere(Eigen::VectorXd &eval, int n)
{
    eval.resize(n);
    int i = 1;
    int band = 1;
    eval(0) = 0.0;
    for (int k = 0; k < 10; k++)
    {
        for (int j = 0; j <= 2 * band; j++)
        {
            eval(i) = -(double)band * ((double)band + 1.0);
            i++;
            if (i == n)
            {
                break;
            }
        }
        if (i == n)
        {
            break;
        }
        band++;
    }
}
//----------------------------------------------------------------------------

double factorial(int n)
{
    if (n == 0)
        return 1.0;
    return (double)n * factorial(n - 1);
}

//----------------------------------------------------------------------------

double scale(int l, int m)
{
    double temp = ((2.0 * (double)l + 1.0) * factorial(l - m)) /
                  (4.0 * M_PI * factorial(l + m));
    return sqrt(temp);
}

//----------------------------------------------------------------------------

double legendre_Polynomial(int l, int m, double x)
{
    // evaluate an Associated Legendre Polynomial P(l,m,x) at x
    double pmm = 1.0;
    if (m > 0)
    {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= m; i++)
        {
            pmm *= (-fact) * somx2;
            fact += 2.0;
        }
    }
    if (l == m)
        return pmm;
    double pmmp1 = x * (2.0 * (double)m + 1.0) * pmm;
    if (l == m + 1)
        return pmmp1;
    double pll = 0.0;
    for (int ll = m + 2; ll <= l; ++ll)
    {
        pll = ((2.0 * (double)ll - 1.0) * x * pmmp1 -
               ((double)ll + (double)m - 1.0) * pmm) /
              ((double)ll - (double)m);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}

//----------------------------------------------------------------------------

double sphericalHarmonic(pmp::Point p, int l, int m)
{

    // l is the band, range [0..n]
    // m in the range [-l..l]
    // transform cartesian to spherical coordinates, assuming r = 1

    double phi = atan2(p[0], p[2]) + M_PI;
    double cos_theta = p[1] / norm(p);
    const double sqrt2 = sqrt(2.0);
    if (m == 0)
        return scale(l, 0) * legendre_Polynomial(l, m, cos_theta);
    else if (m > 0)
        return sqrt2 * scale(l, m) * cos((double)m * phi) *
               legendre_Polynomial(l, m, cos_theta);
    else
        return sqrt2 * scale(l, -m) * sin(-(double)m * phi) *
               legendre_Polynomial(l, -m, cos_theta);
}
//----------------------------------------------------------------------------

double rmse_sh(SurfaceMesh &mesh, unsigned int laplace, unsigned int min_point_,
               bool lumped, int degree, CoarseDimension coarseningType)
{
    auto points = mesh.vertex_property<Point>("v:point");

    // comparing eigenvectors up to the 8th Band of legendre polynomials
    int band = 4;

    double error;
    double sum = 0.0;
    int nv = mesh.n_vertices();
    Eigen::VectorXd y(mesh.n_vertices()), y2(mesh.n_vertices());


    Eigen::SparseMatrix<double> S, M;
    setup_stiffness_matrices(mesh, S, laplace, min_point_);
    setup_mass_matrices(mesh, M, laplace, min_point_, lumped);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);

    int n = 0;
    for (int l = 1; l <= band; l++)
    {
        double eval = -l * (l + 1);
        for (int m = -l; m <= l; m++)
        {
            n++;
            for (auto v : mesh.vertices())
            {
                y(v.idx()) = sphericalHarmonic(points[v], l, m);
                if (l == 4 && m == 2)
                {
                    double x, y_, z;
                    x = points[v][0];
                    y_ = points[v][1];
                    z = points[v][2];

                    y2(v.idx()) = (3.0 / 8.0 * sqrt(5.0 / M_PI)) *
                                  (pow(x, 2.0) - pow(y_, 2.0)) *
                                  (7.0 * pow(z, 2.0) - 1.0);
                }
            }
            y.normalize();
            Eigen::MatrixXd X = solver.solve(S * y);
            Eigen::MatrixXd X2;
            if (l == 4 && m == 2)
            {
                y2.normalize();
                X2 = solver.solve(S * y2);
            }
            error = (y - 1.0 / eval * X).transpose() * M * (y - 1.0 / eval * X);
            if (l == 4 && m == 2)
            {
                std::cout << "error band: " << l << "," << m << ": " << error
                          << std::endl;
                std::cout << "error band explicit: " << l << "," << m << ": "
                          << (y2 - 1.0 / eval * X2).transpose() * M *
                                 (y2 - 1.0 / eval * X2)
                          << std::endl;
            }
            sum += error;
        }
    }
    if (laplace == AlexaWardetzkyLaplace)
    {
        std::cout << "Error SH band recreation  (Alexa, l="
                  << poly_laplace_lambda_ << "): " << sum << std::endl;
    }
    else if (laplace == deGoesLaplace)
    {
        std::cout << "Error SH band recreation  (Disney, l="
                  << disney_laplace_lambda_ << "): " << sum << std::endl;
    }
    else
    {
        if (laplace == Diamond)
        {
            std::cout << "Diamond Laplace: ";
        }
        else if (laplace == PolySimpleLaplace)
        {
            std::cout << "Polysimple Laplace: ";
        }
        if (min_point_ == Centroid)
        {
            std::cout << "Error SH band recreation (centroid): " << sum
                      << std::endl;
        }
        else if (min_point_ == AreaMinimizer)
        {
            std::cout << "Error SH band recreation (area Minimizer): " << sum
                      << std::endl;
        }
    }

    return sum;
}
//=============================================================================
