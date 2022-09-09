#include <iomanip>
#include <igl/slice.h>
#include "Poisson_System.h"
#include "SpectralProcessing.h"

using namespace std;

enum LaplaceMethods {
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
};

enum InsertedPoint {
    Centroid = 0,
    AreaMinimizer = 2
};

enum Function {
    poisson_SH = 0,
    Franke2d = 2,
};

//-----------------------------------------------------------------------------

double solve_poisson_system(pmp::SurfaceMesh &mesh, int laplace, int minpoint,
                            int function, int l, int m) {
#define PMP_SCALAR_TYPE_64

    Eigen::SparseMatrix<double> S, M;
    int nv = (int) mesh.n_vertices();
    Eigen::VectorXd b(nv), analytic_solution(nv);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::MatrixXd X;

    setup_stiffness_matrices(mesh, S, laplace, minpoint);
    setup_mass_matrices(mesh, M, laplace, minpoint, true);
    double error = 0.0;
    if (function == Franke2d) {

        for (auto v: mesh.vertices()) {
            b(v.idx()) =
                    laplace_of_poisson_function(mesh.position(v), function);
        }

        b = M * b;

        // Set the constraints at the locked vertices to the evluation of the Franke function
        for (auto v: mesh.vertices()) {
            if (mesh.is_boundary(v)) {
                // right-hand side: fix boundary values with franke function of the vertices
                b(v.idx()) = poisson_function(mesh.position(v), function);
            }
        }

        // Adjust the right-hand-side to account for the locked nodes
        for (unsigned int i = 0; i < S.outerSize(); i++)
            for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i);
                 iter; ++iter) {
                Vertex row = pmp::Vertex(iter.row());
                Vertex col = pmp::Vertex(iter.col());
                if (!mesh.is_boundary(row) && mesh.is_boundary(col)) {
                    b[iter.row()] -= b[iter.col()] * iter.value();
                }
            }

        // Adjust the system matrix to account for the locked nodes
        for (unsigned int i = 0; i < S.outerSize(); i++)
            for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i);
                 iter; ++iter) {
                Vertex row = pmp::Vertex(iter.row());
                Vertex col = pmp::Vertex(iter.col());
                if (mesh.is_boundary(row))
                    iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
                else if (mesh.is_boundary(col))
                    iter.valueRef() = 0;
            }

        solver.compute(S);
        Eigen::VectorXd x = solver.solve(b);
        for (auto v: mesh.vertices()) {
            error += pow(
                    x[v.idx()] - poisson_function(mesh.position(v), function),
                    2.);
        }

        std::cout << "Franke RMSE error inner vertices: "
                  << sqrt(error / (double) mesh.n_vertices()) << std::endl;
        return sqrt(error / (double) mesh.n_vertices());
    } else if (function == poisson_SH) {
        for (auto v: mesh.vertices()) {
            analytic_solution(v.idx()) =
                    sphericalHarmonic(mesh.position(v), l, m);
        }
        solver.analyzePattern(M);
        solver.factorize(M);
// Astrid
//        analytic_solution.normalize();
// Fernando
        analytic_solution /= sqrt(analytic_solution.transpose() * M * analytic_solution);

        X = solver.solve(S * analytic_solution);
        if (solver.info() != Eigen::Success) {
            std::cout << "Issue: " << solver.info() << std::endl;
            std::cerr << "Could not solve linear system\n";
        }
        double eval = -l * (l + 1);
        error = (analytic_solution - 1.0 / eval * X).transpose() * M *
                (analytic_solution - 1.0 / eval * X);
        error = sqrt(error / double(nv));
        std::cout << "error spherical harmonics "<<"Y_"<<l<<"^"<<m << ": " << error << std::endl;
        return error;
    } else {
        std::cerr << "Function not implemented!" << std::endl;
        return -50.0;
    }
}

//-----------------------------------------------------------------------------

double poisson_function(pmp::Point &p, int function) {

    if (function == Franke2d) {
        return franke_function(p[0], p[1]);
    } else {
        std::cerr << "Function not implemented" << std::endl;
        return -50.0;
    }
}

//-----------------------------------------------------------------------------

double laplace_of_poisson_function(Point &p, int function) {
    if (function == Franke2d) {
        return laplace_franke_function(p[0], p[1]);
    } else {
        std::cerr << "Function not implemented" << std::endl;
        return -50.0;
    }
}

//-----------------------------------------------------------------------------

double franke_function(double x, double y) {
    double cx2 = (9. * x - 2.) * (9. * x - 2.);
    double cy2 = (9. * y - 2.) * (9. * y - 2.);

    double cx1 = (9. * x + 1.) * (9. * x + 1.);
    double cx7 = (9. * x - 7.) * (9. * x - 7.);

    double cy3 = (9. * y - 3.) * (9. * y - 3.);
    double cx4 = (9. * x - 4.) * (9. * x - 4.);

    double cy7 = (9. * y - 7.) * (9. * y - 7.);

    return (3. / 4.) * exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10.) +
           (1. / 2.) * exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3) -
           (1. / 5.) * exp(-cx4 - cy7);
}

//-----------------------------------------------------------------------------

double laplace_franke_function(double x, double y) {

    double mathematica =
            64.8 * exp(-pow(-4. + 9. * x, 2.0) - pow(-7. + 9. * y, 2.0)) -
            40.5 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) -
            60.75 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) -
            1.8720918367346937 * exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
                                     0.1 * (1. + 9. * y)) +
            10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
            pow(-7. + 9. * x, 2) -
            64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
            pow(-4. + 9. * x, 2) +
            15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
            pow(-2. + 9. * x, 2) +
            0.1012078300708038 *
            exp(-0.02040816326530612 * pow(1. + 9. * x, 2) -
                0.1 * (1. + 9. * y)) *
            pow(1. + 9. * x, 2) -
            64.8 * exp(-pow(-4. + 9. * x, 2) - pow(-7. + 9. * y, 2)) *
            pow(-7. + 9. * y, 2) +
            10.125 * exp(0.25 * (-pow(-7. + 9. * x, 2) - pow(-3. + 9. * y, 2))) *
            pow(-3. + 9. * y, 2) +
            15.1875 * exp(0.25 * (-pow(-2. + 9. * x, 2) - pow(-2. + 9. * y, 2))) *
            pow(-2. + 9. * y, 2);
    return mathematica;
}

//-----------------------------------------------------------------------------

