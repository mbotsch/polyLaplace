//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch, Philipp Herholz.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#include "HarmonicBasis.h"
#include "unsupported/Eigen/SparseExtra"
#include "Spectra/MatOp/SymShiftInvert.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Spectra/SymGEigsShiftSolver.h"
#include "Eigenmodes.h"
#include "pmp/Timer.h"
#include "../Surface/diffgeo.h"

double Franke(const Eigen::Vector3d& p)
{
    double x = p(0), y = p(1), z = p(2);

    double cx2 = (9. * x - 2.) * (9. * x - 2.);
    double cy2 = (9. * y - 2.) * (9. * y - 2.);
    double cz2 = (9. * z - 2.) * (9. * z - 2.);

    double cx1 = (9. * x + 1.) * (9. * x + 1.);
    double cx7 = (9. * x - 7.) * (9. * x - 7.);

    double cy3 = (9. * y - 3.) * (9. * y - 3.);
    double cx4 = (9. * x - 4.) * (9. * x - 4.);

    double cy7 = (9. * y - 7.) * (9. * y - 7.);
    double cz5 = (9. * z - 5.) * (9. * z - 5.);

    return (3. / 4.) *
               exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2 - (1. / 4.) * cz2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10. -
                           (9. / 10.) * z - 1. / 10.) +
           (1. / 2.) *
               exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3 - (1. / 4.) * cz5) -
           (1. / 5.) * exp(-cx4 - cy7 - cz5);
}

double FrankeLaplacian(Eigen::Vector3d& p)
{
    double x = p(0), y = p(1), z = p(2);

    return (243. * (-2299. + 1800. * x * (2. + 9. * x))) /
               (480200. *
                exp((9. * (12. + 10. * x * (2. + 9. * x) + 49. * y + 49. * z)) /
                    490.)) -
           (486. *
            exp(-pow(4. - 9. * x, 2) - pow(7. - 9. * y, 2) -
                pow(5. - 9. * z, 2)) *
            (59. + 6. * x * (-8. + 9. * x) + 6. * y * (-14. + 9. * y) +
             6. * z * (-10 + 9 * z))) /
               5. +
           (81. *
            exp((-pow(7. - 9. * x, 2) - 9 * pow(1. - 3. * y, 2) -
                 pow(5. - 9. * z, 2)) /
                4.) *
            (77. + 9. * x * (-14. + 9. * x) + 27. * y * (-2. + 3. * y) +
             9 * z * (-10. + 9. * z))) /
               8. +
           (729. * (2. + 3. * x * (-4. + 9. * x) + 3. * y * (-4. + 9. * y) +
                    3. * z * (-4. + 9. * z))) /
               (16. *
                exp((3. * (4. + 3. * x * (-4. + 9. * x) +
                           3. * y * (-4. + 9. * y) + 3. * z * (-4. + 9. * z))) /
                    4.0));
}

void setup_3D_harmonic_matrices(PolyhedralMesh& mesh,
                                Eigen::SparseMatrix<double>& K,
                                Eigen::SparseMatrix<double>& M)
{
    using namespace std;

    const int nc = mesh.numCells();

    std::vector<Eigen::MatrixXd> allK(nc);
    std::vector<Eigen::MatrixXd> allM(nc);
    std::vector<std::vector<int>> allCellVertices(nc);

    std::vector<HarmonicPolyhedron> harmonicPolyhedra(nc);

    threadHelper<6>(
        [&](const int i) {
            vector<Eigen::Vector3d> vertices;
            vector<vector<int>> faces;
            allCellVertices[i] = mesh.getCellGeometry(i, vertices, faces);
            harmonicPolyhedra[i] = HarmonicPolyhedron(vertices, faces);

            if (i % 10 == 0)
                std::cout << i << " / " << nc << std::endl;
        },
        nc);

    // cannot parallelize quadrature initialization because tetgen seems to be not thread safe
    for (int i = 0; i < nc; ++i)
        harmonicPolyhedra[i].initQuadrature();

    threadHelper<6>(
        [&](const int i) {
            harmonicPolyhedra[i].stiffnessMatrix(allK[i]);
            harmonicPolyhedra[i].massMatrix(allM[i]);

            if (i % 10 == 0)
                std::cout << i << " / " << nc << std::endl;
        },
        nc);

    vector<Eigen::Triplet<double>> tripK, tripM;

    for (int i = 0; i < mesh.numCells(); ++i)
    {
        auto& cellVertices = allCellVertices[i];
        Eigen::MatrixXd& Ki = allK[i];

        for (int j = 0; j < (int)cellVertices.size(); ++j)
        {
            for (int k = 0; k < (int)cellVertices.size(); ++k)
            {
                tripK.emplace_back(cellVertices[j], cellVertices[k], Ki(j, k));
            }
        }

        Eigen::MatrixXd& Mi = allM[i];

        for (int j = 0; j < (int)cellVertices.size(); ++j)
        {
            for (int k = 0; k < (int)cellVertices.size(); ++k)
            {
                tripM.emplace_back(cellVertices[j], cellVertices[k], Mi(j, k));
            }
        }
    }

    K.resize(mesh.numVertices(), mesh.numVertices());
    K.setFromTriplets(tripK.begin(), tripK.end());

    M.resize(mesh.numVertices(), mesh.numVertices());
    M.setFromTriplets(tripM.begin(), tripM.end());
}

void setup_3D_harmonic_stiffness_matrix(PolyhedralMesh& mesh,
                                        Eigen::SparseMatrix<double>& K)
{
    using namespace std;

    const int nc = mesh.numCells();

    std::vector<Eigen::MatrixXd> allK(nc);
    std::vector<Eigen::MatrixXd> allM(nc);
    std::vector<std::vector<int>> allCellVertices(nc);

    std::vector<HarmonicPolyhedron> harmonicPolyhedra(nc);

    threadHelper<6>(
        [&](const int i) {
            vector<Eigen::Vector3d> vertices;
            vector<vector<int>> faces;
            allCellVertices[i] = mesh.getCellGeometry(i, vertices, faces);
            harmonicPolyhedra[i] = HarmonicPolyhedron(vertices, faces);
        },
        nc);

    // cannot parallelize quadrature initialization because tetgen seems to be not thread safe
    for (int i = 0; i < nc; ++i)
        harmonicPolyhedra[i].initQuadrature();

    threadHelper<6>(
        [&](const int i) { harmonicPolyhedra[i].stiffnessMatrix(allK[i]); },
        nc);

    vector<Eigen::Triplet<double>> tripK;

    for (int i = 0; i < mesh.numCells(); ++i)
    {
        auto& cellVertices = allCellVertices[i];
        Eigen::MatrixXd& Ki = allK[i];

        for (int j = 0; j < (int)cellVertices.size(); ++j)
        {
            for (int k = 0; k < (int)cellVertices.size(); ++k)
            {
                tripK.emplace_back(cellVertices[j], cellVertices[k], Ki(j, k));
            }
        }
    }

    K.resize(mesh.numVertices(), mesh.numVertices());
    K.setFromTriplets(tripK.begin(), tripK.end());
}

double solve_3D_Franke_harmonic(const std::string& meshname)
{
    PolyhedralMesh mesh(meshname);
    Eigen::SparseMatrix<double> S, M;
    setup_3D_harmonic_matrices(mesh, S, M);
    lump_matrix(M);
    double count = 0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (Eigen::SparseMatrix<double>::InnerIterator it(M, k); it; ++it)
        {
            count += it.value();
        }
    }
    std::cout << "Volume mass matrix: " << count << std::endl;
    Eigen::VectorXd b(mesh.points.size());

    for (int i = 0; i < (int)mesh.points.size(); i++)
    {
        Eigen::Vector3d p(mesh.points[i][0], mesh.points[i][1],
                          mesh.points[i][2]);
        b(i) = -FrankeLaplacian(p);
    }
    b = M * b;
    for (int i = 0; i < (int)mesh.points.size(); i++)
    {
        if (mesh.isBoundaryVertex(i))
        {
            // right-hand side: fix boundary values with franke function of the vertices
            Eigen::Vector3d p(mesh.points[i][0], mesh.points[i][1],
                              mesh.points[i][2]);
            b(i) = Franke(p);
        }
    }
    // Adjust the right-hand-side to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter;
             ++iter)
        {
            int row = (int)iter.row();
            int col = (int)iter.col();
            if (!mesh.isBoundaryVertex(row) && mesh.isBoundaryVertex(col))
            {
                b[iter.row()] -= b[iter.col()] * iter.value();
            }
        }
    // Adjust the system matrix to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter;
             ++iter)
        {
            int row = (int)iter.row();
            int col = (int)iter.col();
            if (mesh.isBoundaryVertex(row))
                iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
            else if (mesh.isBoundaryVertex(col))
                iter.valueRef() = 0;
        }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;

    solver.compute(S);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "Size x :" << x.size() << std::endl;
    double error = 0.0;
    for (int i = 0; i < (int)mesh.points.size(); i++)
    {
        Eigen::Vector3d p(mesh.points[i][0], mesh.points[i][1],
                          mesh.points[i][2]);
        error += pow(x[i] - Franke(p), 2.);
    }

    std::cout << "DoF " << mesh.points.size() << std::endl;
    std::cout << "Franke RMSE error inner vertices: "
              << sqrt(error / (double)mesh.points.size()) << std::endl;
    return sqrt(error / (double)mesh.points.size());
}

double setup_3D_harmonic_stiffness_matrix(PolyhedralMesh& mesh,
                                          Eigen::SparseMatrix<double>& S,
                                          std::ofstream& file)
{
    using namespace std;
    pmp::Timer t;
    double time = 0.0;

    const int nc = mesh.numCells();

    std::vector<Eigen::MatrixXd> allK(nc);
    std::vector<std::vector<int>> allCellVertices(nc);

    std::vector<HarmonicPolyhedron> harmonicPolyhedra(nc);
    t.start();
    threadHelper<8>(
        [&](const int i) {
            vector<Eigen::Vector3d> vertices;
            vector<vector<int>> faces;
            allCellVertices[i] = mesh.getCellGeometry(i, vertices, faces);
            harmonicPolyhedra[i] = HarmonicPolyhedron(vertices, faces);

            if (i % 10 == 0)
                std::cout << i << " / " << nc << std::endl;
        },
        nc);
    t.stop();
    time += t.elapsed();
    file << t.elapsed() << ",";
    t.start();
    // cannot parallelize quadrature initialization because tetgen seems to be not thread safe
    for (int i = 0; i < nc; ++i)
        harmonicPolyhedra[i].initQuadrature();
    t.stop();
    time += t.elapsed();
    file << t.elapsed() << ",";
    t.start();
    threadHelper<8>(
        [&](const int i) { harmonicPolyhedra[i].stiffnessMatrix(allK[i]); },
        nc);

    vector<Eigen::Triplet<double>> tripK;

    for (int i = 0; i < mesh.numCells(); ++i)
    {
        auto& cellVertices = allCellVertices[i];
        Eigen::MatrixXd& Ki = allK[i];

        for (int j = 0; j < (int)cellVertices.size(); ++j)
        {
            for (int k = 0; k < (int)cellVertices.size(); ++k)
            {
                tripK.emplace_back(cellVertices[j], cellVertices[k], Ki(j, k));
            }
        }
    }

    S.resize(mesh.numVertices(), mesh.numVertices());
    S.setFromTriplets(tripK.begin(), tripK.end());
    t.stop();
    time += t.elapsed();

    file << t.elapsed() << ",";
    file << time << std::endl;
    std::cout << "Stiffness matrix construction: " << time << std::endl;
    return time;
}

double solve_harmonic_eigenvalue_problem(std::string& meshname,
                                         Eigen::VectorXd& evalues,
                                         std::string& meshname_file)
{
    PolyhedralMesh mesh(meshname);
    Eigen::SparseMatrix<double> S, M;
    setup_3D_harmonic_matrices(mesh, S, M);
    lump_matrix(M);
    S *= -1.0;
    std::string filename = "eigenmodes_[MKB08]" + meshname_file + ".csv";
    std::ofstream ev_file(filename);
    ev_file << "computed,analytic,offset" << std::endl;

    Eigen::VectorXd b(mesh.numVertices());

    // collect indices of inner vertices
    std::vector<int> indices;
    for (int i = 0; i < (int)mesh.points.size(); i++)
    {
        if (!mesh.isBoundaryVertex(i))
        {
            indices.push_back(i);
        }
    }

    Eigen::VectorXi in(indices.size());
    std::cout << "inner indices: " << indices.size() << std::endl;
    //Rewrite indices to Eigen::Vector
    for (int i = 0; i < (int)indices.size(); ++i)
    {
        in(i) = indices[i];
    }

    Eigen::SparseMatrix<double> S_in_in, M_in_in;

    //slice matrices so that only rows and cols for inner vertices remain

    Eigen::SparseMatrix<double> column_subset_S(S.rows(), (int)indices.size());
    Eigen::SparseMatrix<double> column_subset_M(M.rows(), (int)indices.size());

    Eigen::SparseMatrix<double, Eigen::RowMajor> row_subset_S(
        (int)indices.size(), (int)indices.size());
    Eigen::SparseMatrix<double, Eigen::RowMajor> row_subset_M(
        (int)indices.size(), (int)indices.size());

    for (int j = 0; j != column_subset_S.cols(); ++j)
    {
        column_subset_S.col(j) = S.col(indices[j]);
        column_subset_M.col(j) = M.col(indices[j]);
    }
    for (int j = 0; j != row_subset_S.rows(); ++j)
    {
        row_subset_S.row(j) = column_subset_S.row(indices[j]);
        row_subset_M.row(j) = column_subset_M.row(indices[j]);
    }
    S_in_in = row_subset_S;
    M_in_in = row_subset_M;

    //--------------

    int num_eval = 34;
    int converge_speed = 5 * num_eval;
    // S and M are sparse
    using OpType =
        Spectra::SymShiftInvert<double, Eigen::Sparse, Eigen::Sparse>;
    using BOpType = Spectra::SparseSymMatProd<double>;
    OpType op(S_in_in, M_in_in);
    BOpType Bop(M_in_in);

    // Construct generalized eigen solver object, seeking three generalized
    // eigenvalues that are closest to zero. This is equivalent to specifying
    // a shift sigma = 0.0 combined with the SortRule::LargestMagn selection rule
    Spectra::SymGEigsShiftSolver<OpType, BOpType,
                                 Spectra::GEigsMode::ShiftInvert>
        geigs(op, Bop, num_eval, converge_speed, 0.0);
    geigs.init();
    geigs.compute(Spectra::SortRule::LargestMagn);

    Eigen::VectorXd evectors, analytic;
    // Retrieve results
    if (geigs.info() == Spectra::CompInfo::Successful)
    {
        evalues = geigs.eigenvalues();
    }
    analytic_eigenvalues_unitBall(analytic, num_eval);
    double error = 0.0;
    for (int i = 0; i < evalues.size(); i++)
    {
        ev_file << -evalues(i) << "," << analytic(i) << ","
                << -evalues(i) - analytic(i) << std::endl;
        std::cout << "Computed evalue: " << -evalues(i)
                  << " analytical Bessel: " << analytic(i) << std::endl;

        error += pow(-evalues(i) - analytic(i), 2);
    }
    error = sqrt(error / (double)evalues.size());
    std::cout << "Root mean squared error: " << error << std::endl;
    return error;
}
