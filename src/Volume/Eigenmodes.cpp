
#include "Eigenmodes.h"
#include <cmath>
#include "LaplaceConstruction_3D.h"
#include "unsupported/Eigen/SparseExtra"
#include <igl/slice.h>
#include <Spectra/MatOp/SparseGenMatProd.h>
#include <Spectra/MatOp/SparseSymMatProd.h>
#include <Spectra/Util/GEigsMode.h>
#include <fstream>
#include <Spectra/Util/SelectionRule.h>
#include <Spectra/MatOp/SymShiftInvert.h>
#include <Spectra/SymGEigsShiftSolver.h>

enum LaplaceMethods {
    Diamond = 0,
    PolySimpleLaplace = 2,
};

enum VolumePoints {
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

double solve_eigenvalue_problem(VolumeMesh &mesh, Eigen::VectorXd &evalues,
                                int Laplace, int face_point, int cell_point,
                                const std::string& meshname) {
    std::string filename;
    if (Laplace == Diamond) {
        filename = "eigenmodes_[BBA21]" + meshname + ".csv";
    } else if (Laplace == PolySimpleLaplace) {

        filename = "eigenmodes_[BHBK20]" + meshname + ".csv";

    }

    std::ofstream ev_file(filename);
    ev_file << "computed,analytic,offset" << std::endl;

    Eigen::SparseMatrix<double> M, S;
    Eigen::VectorXd b(mesh.n_vertices());

    setup_3D_stiffness_matrix(mesh, S, Laplace, face_point, cell_point);
    setup_3D_mass_matrix(mesh, M, Laplace, face_point, cell_point);
    // collect indices of inner vertices
    std::vector<int> indices;
    for (auto v: mesh.vertices()) {
        if (!mesh.is_boundary(v)) {
            indices.push_back(v.idx());
        }
    }
    Eigen::VectorXi in(indices.size());
    std::cout << "inner indices: " << indices.size() << std::endl;
    //Rewrite indices to Eigen::Vector
    for (int i = 0; i < (int) indices.size(); ++i) {
        in(i) = indices[i];
    }

    Eigen::SparseMatrix<double> S_in_in, M_in_in;

    //slice matrices so that only rows and cols for inner vertices remain

    igl::slice(S, in, in, S_in_in);
    igl::slice(M, in, in, M_in_in);

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
    if (geigs.info() == Spectra::CompInfo::Successful) {
        evalues = geigs.eigenvalues();
    }
    analytic_eigenvalues_unitBall(analytic, num_eval);
    double error = 0.0;
    for (int i = 0; i < evalues.size(); i++) {
        ev_file << -evalues(i) << "," << analytic(i) << ","
                << -evalues(i) - analytic(i) << std::endl;
        std::cout << "Computed evalue: " << -evalues(i)
                  << " analytical Bessel: " << analytic(i) << std::endl;

        error += pow(-evalues(i) - analytic(i), 2);
    }
    error = sqrt(error / (double) evalues.size());
    std::cout << "Root mean squared error: " << error << std::endl;
    return error;
}

void analytic_eigenvalues_unitBall(Eigen::VectorXd &eval, int n) {

    //n so far only <= 34
    if (n > 34) {
        std::cout << "n must be lower than 34! \n";
    }
    eval.resize(n);
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            //Bessel j zero (1/2,1)
            eval(i) = 9.869604401089358;
        } else if (i > 0 && i < 4) {
            //Bessel j zero (3/2,1)
            eval(i) = 20.190728556426629;
        } else if (i >= 4 && i < 9) {
            //Bessel j zero (5/2,1)
            eval(i) = 33.217461914268368;
        } else if (i == 9) {
            //Bessel j zero (1/2,2)

            eval(i) = 39.47841760435743447;
        } else if (i > 9 && i < 17) {
            //Bessel j zero (7/2,1)
            eval(i) = 48.831193643619198876;
        } else if (i >= 17 && i < 20) {
            //Bessel j zero (3/2,2)
            eval(i) = 59.6795159441094188805;
        } else if (i >= 20 && i < 29) {
            //Bessel j zero (9/2,1)
            eval(i) = 66.95431192510480532587;
        } else if (i >= 29 && i < 34) {
            //Bessel j zero (5/2,2)
            eval(i) = 82.719231101493279988217;
        }
    }
}
