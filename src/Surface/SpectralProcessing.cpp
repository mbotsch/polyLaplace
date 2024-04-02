//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "../common_util.h"
#include "SpectralProcessing.h"
#include "LaplaceConstruction.h"
#include "[dGBD20]Laplace.h"
#include "HarmonicBasis2D.h"
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/MatOp/SparseSymShiftSolve.h>
#include <Spectra/MatOp/SparseCholesky.h>
#include <Spectra/Util/GEigsMode.h>
#include <iostream>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/Util/CompInfo.h>
#include <Spectra/SymGEigsShiftSolver.h>
#include "Spectra/SymGEigsSolver.h"
#include <Spectra/SymEigsShiftSolver.h>
#include "Spectra/SymEigsSolver.h"
#include <iomanip>

//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

double solve_eigenvalue_problem(SurfaceMesh& mesh, int laplace, int face_point,
                                const std::string& meshname)
{
    std::string filename;
    if (laplace == Diamond)
    {
        filename = "eigenvalues_[BBA21]_" + meshname + ".csv";
    }
    else if (laplace == AlexaWardetzkyLaplace)
    {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << poly_laplace_lambda_;
        std::string s = stream.str();
        filename = "eigenvalues_[AW11]_l=" + s + "_" + meshname + ".csv";
    }
    else if (laplace == deGoesLaplace)
    {
        std::stringstream stream;
        stream << std::fixed << std::setprecision(2) << deGoes_laplace_lambda_;
        std::string s = stream.str();
        filename = "eigenvalues_[dGBD20]_l=" + s + "_" + meshname + ".csv";
    }
    else if (laplace == PolySimpleLaplace)
    {
        filename = "eigenvalues_[BHKB20]_" + meshname + ".csv";
    }
    else if (laplace == Harmonic)
    {
        filename = "eigenvalues_[MKB08]_" + meshname + ".csv";
    }
    std::ofstream ev_file(filename);

    ev_file << "computed,analytic,offset" << std::endl;
    Eigen::SparseMatrix<double> M, S;
    double error = 0.0;
    if (laplace == Harmonic)
    {
        buildStiffnessAndMass2d(mesh, S, M);
        lump_matrix(M);
        S *= -1.0;
    }
    else
    {
        setup_stiffness_matrices(mesh, S, laplace, face_point);
        setup_mass_matrices(mesh, M, laplace, face_point);
    }
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

    for (int i = 1; i < evalues.size(); i++)
    {
        ev_file << evalues(i) << "," << analytic(i) << ","
                << evalues(i) - analytic(i) << std::endl;
        error += pow(evalues(i) - analytic(i), 2);
    }

    error = sqrt(error / (double)evalues.size());
    std::cout << "Root mean squared error: " << error << std::endl;
    ev_file.close();

    return error;
}

void analytic_eigenvalues_unitsphere(Eigen::VectorXd& eval, int n)
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
                  (4.0 * std::numbers::pi * factorial(l + m));
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

    double phi = atan2(p[0], p[2]) + std::numbers::pi;
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

double rmse_sh(SurfaceMesh& mesh, int laplace, int min_point_, bool lumped)
{
    auto points = mesh.vertex_property<Point>("v:point");

    // comparing eigenvectors up to the 8th Band of legendre polynomials
    int band = 3;

    double error;
    double sum = 0.0;
    Eigen::VectorXd y(mesh.n_vertices());

    Eigen::SparseMatrix<double> S, M;
    setup_stiffness_matrices(mesh, S, laplace, min_point_);
    setup_mass_matrices(mesh, M, laplace, min_point_, lumped);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.analyzePattern(M);
    solver.factorize(M);
    for (int l = 1; l <= band; l++)
    {
        double eval = -l * (l + 1);
        for (int m = -l; m <= l; m++)
        {
            for (auto v : mesh.vertices())
            {
                y(v.idx()) = sphericalHarmonic(points[v], l, m);
            }
            Eigen::MatrixXd X = solver.solve(S * y);
            error = (y - 1.0 / eval * X).transpose() * M * (y - 1.0 / eval * X);
            error = sqrt(error / double(mesh.n_vertices()));
            sum += error;
        }
    }

    if (laplace == AlexaWardetzkyLaplace)
    {
        std::cout << "Error SH band recreation  (AlexaWardetzky Laplace, l="
                  << poly_laplace_lambda_ << "): " << sum << std::endl;
    }
    else if (laplace == deGoesLaplace)
    {
        std::cout << "Error SH band recreation  (deGoes Laplace, l="
                  << deGoes_laplace_lambda_ << "): " << sum << std::endl;
    }
    else if (laplace == Harmonic)
    {
        std::cout << "Error SH band recreation  (Harmonic Laplace: " << sum
                  << std::endl;
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
        if (min_point_ == Centroid_)
        {
            std::cout << "Error SH band recreation (centroid): " << sum
                      << std::endl;
        }
        else if (min_point_ == AreaMinimizer)
        {
            std::cout << "Error SH band recreation (area Minimizer): " << sum
                      << std::endl;
        }
        else
        {
            std::cout << "Error SH band recreation (trace minimizer): " << sum
                      << std::endl;
        }
    }
    return sum;
}

double condition_number(pmp::SurfaceMesh& mesh, int laplace, int minpoint,
                        Eigen::Vector3d& values, bool generalized)
{
    const int numVerts = mesh.n_vertices();
    Eigen::SparseMatrix<double> S, M;
    if (laplace == Harmonic)
    {
        buildStiffnessAndMass2d(mesh, S, M);
        lump_matrix(M);
        S *= -1.0;
    }
    else
    {
        setup_stiffness_matrices(mesh, S, laplace, minpoint);
        setup_mass_matrices(mesh, M, laplace, minpoint, true);
    }

    //slice matrices so that only rows and cols for inner vertices remain
    std::vector<int> innerVertIdxs;
    for (auto v : mesh.vertices())
    {
        if (!mesh.is_boundary(v))
        {
            innerVertIdxs.push_back(v.idx());
        }
    }
    int nInnerVertIdxs = innerVertIdxs.size();

    Eigen::SparseMatrix<double> S_in_in(nInnerVertIdxs, nInnerVertIdxs);
    Eigen::SparseMatrix<double> M_in_in(nInnerVertIdxs, nInnerVertIdxs);
    if (nInnerVertIdxs == numVerts)
    {
        S_in_in = S;
        M_in_in = M;
    }
    else
    {
        Eigen::SparseMatrix<double> S_columns(S.rows(), nInnerVertIdxs);
        Eigen::SparseMatrix<double> M_columns(M.rows(), nInnerVertIdxs);
        Eigen::SparseMatrix<double, Eigen::RowMajor> S_rows(nInnerVertIdxs,
                                                            nInnerVertIdxs);
        Eigen::SparseMatrix<double, Eigen::RowMajor> M_rows(nInnerVertIdxs,
                                                            nInnerVertIdxs);

        // process rows and columns separately for linear runtime
        for (int i = 0; i < nInnerVertIdxs; i++)
        {
            S_columns.col(i) = S.col(innerVertIdxs[i]);
            M_columns.col(i) = M.col(innerVertIdxs[i]);
        }
        for (int i = 0; i < nInnerVertIdxs; i++)
        {
            S_rows.row(i) = S_columns.row(innerVertIdxs[i]);
            M_rows.row(i) = M_columns.row(innerVertIdxs[i]);
        }
        S_in_in = S_rows;
        M_in_in = M_rows;
    }

    int numEigValues = 3;
    int convergenceSpeed = std::min(40 * numEigValues, (int)S_in_in.rows());

    Eigen::VectorXd eigValsMax;
    Eigen::VectorXd eigValsMin;
    if (generalized)
    {
        // Construct generalized eigen solver object, requesting the largest generalized eigenvalue
        Spectra::SparseSymMatProd<double> sOpMax(-S_in_in);
        Spectra::SparseCholesky<double> sBOpMax(M_in_in);
        Spectra::SymGEigsSolver<Spectra::SparseSymMatProd<double>,
                                Spectra::SparseCholesky<double>,
                                Spectra::GEigsMode::Cholesky>
            eigSolverMax(sOpMax, sBOpMax, 1, convergenceSpeed);
        eigSolverMax.init();
        eigSolverMax.compute(Spectra::SortRule::LargestAlge);
        eigValsMax = eigSolverMax.eigenvalues();

        switch (eigSolverMax.info())
        {
            case Spectra::CompInfo::NotComputed:
                std::cout << "Max Eig: Not Computed" << std::endl;
                break;
            case Spectra::CompInfo::NotConverging:
                std::cout << "Max Eig: Not Converging" << std::endl;
                break;
            case Spectra::CompInfo::NumericalIssue:
                std::cout << "Max Eig: Numerical Issue" << std::endl;
                break;
            default:
                break;
        }

        // Construct generalized eigen solver object, seeking three generalized
        // eigenvalues that are closest to zero. This is equivalent to specifying
        // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
        Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse> sOpMin(
            -S_in_in, M_in_in);
        Spectra::SparseSymMatProd<double> sBOpMin(M_in_in);
        Spectra::SymGEigsShiftSolver<
            Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>,
            Spectra::SparseSymMatProd<double>, Spectra::GEigsMode::ShiftInvert>
            eigSolverMin(sOpMin, sBOpMin, numEigValues, convergenceSpeed, -0.1);
        eigSolverMin.init();
        eigSolverMin.compute(Spectra::SortRule::LargestMagn);
        eigValsMin = eigSolverMin.eigenvalues();

        switch (eigSolverMin.info())
        {
            case Spectra::CompInfo::NotComputed:
                std::cout << "Min Eig: Not Computed" << std::endl;
                break;
            case Spectra::CompInfo::NotConverging:
                std::cout << "Min Eig: Not Converging" << std::endl;
                break;
            case Spectra::CompInfo::NumericalIssue:
                std::cout << "Min Eig: Numerical Issue" << std::endl;
                break;
            default:
                break;
        }
    }
    else
    {
        // Max Eigenvalue solver
        Spectra::SparseSymMatProd<double> sOpMax(-S_in_in);
        Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigSolverMax(
            sOpMax, 1, convergenceSpeed);
        eigSolverMax.init();
        eigSolverMax.compute(Spectra::SortRule::LargestAlge);
        eigValsMax = eigSolverMax.eigenvalues();

        switch (eigSolverMax.info())
        {
            case Spectra::CompInfo::NotComputed:
                std::cout << "Not Computed" << std::endl;
                break;
            case Spectra::CompInfo::NotConverging:
                std::cout << "Not Converging" << std::endl;
                break;
            case Spectra::CompInfo::NumericalIssue:
                std::cout << "Numerical Issue" << std::endl;
                break;
            default:
                break;
        }

        // Min Eigenvalue solver
        Spectra::SparseSymShiftSolve<double> sOpMin(-S_in_in);
        Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>>
            eigSolverMin(sOpMin, numEigValues, convergenceSpeed, -0.1);
        eigSolverMin.init();
        eigSolverMin.compute(Spectra::SortRule::LargestMagn);
        eigValsMin = eigSolverMin.eigenvalues();

        switch (eigSolverMin.info())
        {
            case Spectra::CompInfo::NotComputed:
                std::cout << "Not Computed" << std::endl;
                break;
            case Spectra::CompInfo::NotConverging:
                std::cout << "Not Converging" << std::endl;
                break;
            case Spectra::CompInfo::NumericalIssue:
                std::cout << "Numerical Issue" << std::endl;
                break;
            default:
                break;
        }
    }

    values(0) = eigValsMax.coeff(0);
    values(1) =
        eigValsMin.coeff(numEigValues - 1 - (numVerts == innerVertIdxs.size()));
    values(2) = values(0) / values(1);
    return values(2);
}
//=============================================================================

double get_condition_number(const Eigen::SparseMatrix<double>& M,
                            bool firstEigZero)
{
    int numEigValues = 3;
    int convergenceSpeed = std::min(40 * numEigValues, (int)M.rows());

    Eigen::VectorXd eigValsMax;
    Eigen::VectorXd eigValsMin;

    Spectra::SparseSymMatProd<double> sOpMax(M);
    Spectra::SymEigsSolver<Spectra::SparseSymMatProd<double>> eigSolverMax(
        sOpMax, 1, convergenceSpeed);
    eigSolverMax.init();
    eigSolverMax.compute(Spectra::SortRule::LargestAlge);
    eigValsMax = eigSolverMax.eigenvalues();

    switch (eigSolverMax.info())
    {
        case Spectra::CompInfo::NotComputed:
            std::cout << "Not Computed" << std::endl;
            break;
        case Spectra::CompInfo::NotConverging:
            std::cout << "Not Converging" << std::endl;
            break;
        case Spectra::CompInfo::NumericalIssue:
            std::cout << "Numerical Issue" << std::endl;
            break;
        default:
            break;
    }

    // Min Eigenvalue solver
    Spectra::SparseSymShiftSolve<double> sOpMin(M);
    Spectra::SymEigsShiftSolver<Spectra::SparseSymShiftSolve<double>>
        eigSolverMin(sOpMin, numEigValues, convergenceSpeed, -0.1);
    eigSolverMin.init();
    eigSolverMin.compute(Spectra::SortRule::LargestMagn);
    eigValsMin = eigSolverMin.eigenvalues();

    switch (eigSolverMin.info())
    {
        case Spectra::CompInfo::NotComputed:
            std::cout << "Not Computed" << std::endl;
            break;
        case Spectra::CompInfo::NotConverging:
            std::cout << "Not Converging" << std::endl;
            break;
        case Spectra::CompInfo::NumericalIssue:
            std::cout << "Numerical Issue" << std::endl;
            break;
        default:
            break;
    }

    //std::cout << eigValsMin << std::endl << eigValsMax << std::endl;
    std::cout << "Condition number of application matrix: "
              << eigValsMax.coeff(0) /
                     eigValsMin.coeff(numEigValues - 1 - (firstEigZero))
              << std::endl;
    return eigValsMax.coeff(0) /
           eigValsMin.coeff(numEigValues - 1 - (firstEigZero));
}