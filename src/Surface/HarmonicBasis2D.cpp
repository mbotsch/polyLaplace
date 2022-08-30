#include "HarmonicBasis2D.h"
#include "unsupported/Eigen/SparseExtra"
#include "Spectra/MatOp/SymShiftInvert.h"
#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Spectra/SymGEigsShiftSolver.h"
#include "Poisson_System.h"
#include "Eigen/CholmodSupport"
#include "SpectralProcessing.h"
#include "../HarmonicBase3D2D/HarmonicPolygon.hpp"

void buildStiffnessAndMass2d(pmp::SurfaceMesh &mesh, Eigen::SparseMatrix<double>& S, Eigen::SparseMatrix<double>& M)
{
    Eigen::MatrixXd V(mesh.n_vertices(),3);
    std::vector<std::vector<int>> poly;

    for(auto v: mesh.vertices()){
        Eigen::Vector3d p (mesh.position(v)[0],mesh.position(v)[1],mesh.position(v)[2]);
        V.row(v.idx()) = p;
    }
    for(auto f : mesh.faces()){
        std::vector<int> face;
        for (auto fv: mesh.vertices(f)){
            face.emplace_back(fv.idx());
        }
        poly.emplace_back(face);
    }

    const int nv = (int)V.rows();
    std::vector<Eigen::Triplet<double>> tripS, tripM;

    for(auto& p : poly) {
        Eigen::MatrixXd pts(p.size(), 3);
        for(unsigned int i = 0; i < p.size(); ++i) {
            pts(i, 0) = V(p[i], 0);
            pts(i, 1) = V(p[i], 1);
            pts(i, 2) = V(p[i], 2);
        }

        HarmonicPolygon hp(pts);

        Eigen::MatrixXd Si, Mi;

        hp.stiffnessMatrix(Si);
        hp.massMatrix(Mi);

        for(unsigned int i = 0; i < p.size(); ++i) {
            for(unsigned int j = 0; j < p.size(); ++j) {
                tripS.emplace_back(p[i], p[j], Si(i, j));
                tripM.emplace_back(p[i], p[j], Mi(i, j));
            }
        }
    }

    S.resize(nv, nv);
    S.setFromTriplets(tripS.begin(), tripS.end());

    M.resize(nv, nv);
    M.setFromTriplets(tripM.begin(), tripM.end());
}

double solve_2D_Franke_harmonic(pmp::SurfaceMesh &mesh)
{
    Eigen::SparseMatrix<double> S, M;
    buildStiffnessAndMass2d(mesh, S, M);

    Eigen::VectorXd b(mesh.n_vertices());

    for (auto v: mesh.vertices())
    {
        b(v.idx()) = -laplace_franke_function(mesh.position(v)[0], mesh.position(v)[1]);
    }
    b = M * b;
    for (auto v: mesh.vertices())
    {
        if(mesh.is_boundary(v)){
            // right-hand side: fix boundary values with franke function of the vertices
            b(v.idx()) = franke_function(mesh.position(v)[0], mesh.position(v)[1]);
        }
    }
    // Adjust the right-hand-side to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter;
             ++iter)
        {
            pmp::Vertex row = pmp::Vertex(iter.row());
            pmp::Vertex col = pmp::Vertex(iter.col());
            if (!mesh.is_boundary(row) && mesh.is_boundary(col))
            {
                b[iter.row()] -= b[iter.col()] * iter.value();
            }
        }
    // Adjust the system matrix to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter;
             ++iter)
        {
            pmp::Vertex row = pmp::Vertex(iter.row());
            pmp::Vertex col = pmp::Vertex(iter.col());
            if (mesh.is_boundary(row)) iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
            else if (mesh.is_boundary(col)) iter.valueRef() = 0;
        }

    //    Eigen::SimplicialLLT< Eigen::SparseMatrix< double > > solver;
    Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver;

    solver.compute(S);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "Size x :" << x.size() << std::endl;
    double error = 0.0;
    int k = 0;
    if (solver.info() != Eigen::Success) {
        std::cerr << "harmonic(): Could not solve linear system\n";
    } else {
        // copy solution
        for (auto v: mesh.vertices()) {
            if (!mesh.is_boundary(v)) {
                error += pow(x(v.idx()) - franke_function(mesh.position(v)[0],mesh.position(v)[1]), 2.0);
                k++;
            }
        }
    }
//    for (auto v: mesh.vertices())
//    {
////        std::cout << "x: " << x[v.idx()] << " Franke : " << franke_function(mesh.position(v)[0], mesh.position(v)[1]) << std::endl;
//        error += pow(x[v.idx()] - franke_function(mesh.position(v)[0], mesh.position(v)[1]), 2.);
//    }

    std::cout << "DoF " << mesh.n_vertices() << std::endl;
    std::cout << "Franke RMSE error inner vertices: "
              << sqrt(error / (double)k) << std::endl;
    return sqrt(error / (double)k);
}

double solve_2D_harmonic_eigenvalue_problem(pmp::SurfaceMesh &mesh, Eigen::VectorXd &evalues, std::string& meshname_file)
{

    Eigen::SparseMatrix<double> S, M;
    buildStiffnessAndMass2d(mesh, S, M);
    S*=-1.0;
    std::string filename = "eigenvalues_Harmonic_" + meshname_file + ".csv";
    std::ofstream ev_file(filename);
    ev_file << "computed,analytic,offset" << std::endl;

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

    Spectra::SymGEigsShiftSolver<OpType, BOpType,Spectra::GEigsMode::ShiftInvert>geigs(op, Bop, num_eval, converge_speed, 1e-8);
    geigs.init();
    geigs.compute(Spectra::SortRule::LargestMagn);
    std::cout << "compute & init" << std::endl;

    Eigen::VectorXd analytic;
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
