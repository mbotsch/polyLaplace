#include <iomanip>
#include <igl/slice.h>
#include "Poisson_System.h"
#include "SpectralProcessing.h"
#include "HarmonicBasis2D.h"

using namespace std;

enum LaplaceMethods
{
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
    Harmonic = 4
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
    int nv = (int)mesh.n_vertices();
    Eigen::VectorXd b(nv), analytic_solution(nv);

    if (laplace == Harmonic) {
        return solve_2D_Franke_harmonic(mesh);
    }
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    Eigen::MatrixXd X;

    setup_stiffness_matrices(mesh, S, laplace, minpoint);
    setup_mass_matrices(mesh, M, laplace, minpoint, true);
    double error = 0.0;
    if (function != poisson_SH) {

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
    } else {
        int k = nv;
        for (auto v: mesh.vertices()) {
            analytic_solution(v.idx()) =
                    sphericalHarmonic(mesh.position(v), l, m);
        }
        solver.analyzePattern(M);
        solver.factorize(M);
        analytic_solution.normalize();
        X = solver.solve(S * analytic_solution);
        if (solver.info() != Eigen::Success) {
            std::cout << "Issue: " << solver.info() << std::endl;
            std::cerr << "Could not solve linear system\n";
        }
        double eval = -l * (l + 1);
        error = (analytic_solution - 1.0 / eval * X).transpose() * M *
                (analytic_solution - 1.0 / eval * X);
        error = sqrt(error / double(k));
        std::cout << "error inner product : " << error << std::endl;
        return error;
    }
 }

//-----------------------------------------------------------------------------

double poisson_function(pmp::Point &p, int function) {

    if (function == Franke2d) {
        return franke_function(p[0], p[1]);
    } else if (function == poisson_SH) {
        return spherical_harmonic_function_scaled(p[0], p[1], p[2]);
    }else{
        std::cerr << "Function not implemented" << std::endl;
    }
}

//-----------------------------------------------------------------------------

double laplace_of_poisson_function(Point &p, int function) {
    if (function == Franke2d) {
        return laplace_franke_function(p[0], p[1]);
    } else if (function == poisson_SH) {
        return spherical_harmonic_function(p[0], p[1], p[2]);
    }else{
        std::cerr << "Function not implemented" << std::endl;
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

double spherical_harmonic_function(double x, double y, double z) {

    //    y_4,-3 (evalue 20)
    //    return -(3*pow(x, 2.0) - pow(y, 2.0))*y*z; //* (3.0/8.0*sqrt(5.0/M_PI));#

    //y_4,2 (evalue 20)
    //    return -(pow(x, 2.0) - pow(y, 2.0)) * (7.0 * pow(z, 2.0) - 1.0); //* (3.0/8.0*sqrt(5.0/M_PI));#

    // y_3,-1 (evalue 12)
    //    return -y * (4.0 * pow(z, 2.0) - pow(x, 2.0) - pow(y, 2.0));

    // y_2,0 (evalue 6)
    //   return -2.0*pow(z, 2.0)+pow(x, 2.0)+pow(y, 2.0);

    // y_3,0 (evalue 12)
    return -z * (2.0 * pow(z, 2.0) - 3.0 * pow(x, 2.0) - 3.0 * pow(y, 2.0));
}

//-----------------------------------------------------------------------------

double spherical_harmonic_function_scaled(double x, double y, double z) {
    //    double l = 4.0;
    double l = 3.0;
    //    double l = 2.0;

    double evalue = -l * (l + 1.0);

    return spherical_harmonic_function(x, y, z) / evalue;
}

//-----------------------------------------------------------------------------

void solve_laplace_equation(SurfaceMesh &mesh, int laplace, int face_point) {
    Eigen::SparseMatrix<double> M, S, S_f;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero((int)mesh.n_vertices(), 3);

    setup_stiffness_matrices(mesh, S, laplace, face_point);
    int nb = 0;
    for (auto v: mesh.vertices()) {
        //count nr outer vertices
        if (mesh.is_boundary(v)) {
            nb++;
        }
    }

    Eigen::SparseMatrix<double> G, V, Div, P;

    unsigned ins = 0;
    unsigned out = 0;

    Eigen::VectorXi in(mesh.n_vertices() - nb), bound(nb), x(3);
    x << 0, 1, 2;
    for (auto v: mesh.vertices()) {
        // save indices of inner and outer vertices
        if (!mesh.is_boundary(v)) {
            in(ins) = (int)v.idx();
            ins++;
        } else {
            bound(out) = (int)v.idx();
            out++;
            // right-hand side: fix x coordinate of boundary vertices for the righthandsite
            B(v.idx(), 0) = mesh.position(v)[0];
            B(v.idx(), 1) = mesh.position(v)[1];
            B(v.idx(), 2) = mesh.position(v)[2];
        }
    }

    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out;
    // slice S and b and bring boundary values on the righthandsite to solve only for inner vertices
    igl::slice(S, in, in, L_in_in);
    igl::slice(S, in, bound, L_in_b);
    igl::slice(B, in, x, b_in);
    igl::slice(B, bound, x, b_out);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_in_in);
    Eigen::MatrixXd X = solver.solve(b_in - L_in_b * b_out);

    double error = 0.0;
    int k = 0;
    if (solver.info() != Eigen::Success) {
        std::cerr << "harmonic(): Could not solve linear system\n";
    } else {
        // copy solution
        for (auto v: mesh.vertices()) {
            if (!mesh.is_boundary(v)) {
                error += pow(X(k, 0) - mesh.position(v)[0], 2) +
                         pow(X(k, 1) - mesh.position(v)[1], 2) +
                         pow(X(k, 2) - mesh.position(v)[2], 2);
                k++;
            }
        }
    }
    std::cout << "RMSE inner vertex positions : " << sqrt(error / (double) k)
              << std::endl;
}
