//=============================================================================

#include "diffgeo.h"
#include <Eigen/Dense>
#include "[AW11]Laplace.h"
#include <unsupported/Eigen/NonLinearOptimization>

const double eps = 1e-10;

#include <pmp/Timer.h>

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

enum InsertedPoint {
    Centroid = 0,
    AreaMinimizer = 2,
};

//=================== Setup P matrix ==========================================================

void setup_prolongation_matrix(pmp::SurfaceMesh &mesh, SparseMatrix &P) {

    auto area_weights = mesh.get_face_property<Eigen::VectorXd>("f:weights");

    const unsigned int nv = mesh.n_vertices();
    const unsigned int nf = mesh.n_faces();
    Eigen::VectorXd w;

    std::vector<Triplet> tripletsA;
    pmp::Vertex v;
    pmp::Face f;
    for (auto v : mesh.vertices()) {
        tripletsA.emplace_back(v.idx(), v.idx(), 1.0);
    }

    unsigned int j = 0;
    for (auto f : mesh.faces()) {
        w = area_weights[f];
        unsigned int i = 0;
        for (auto v : mesh.vertices(f)) {
            tripletsA.emplace_back(nv + j, v.idx(), w(i));
            i++;
        }
        j++;
    }

    // build sparse matrix from triplets
    P.resize(nv + nf, nv);
    P.setFromTriplets(tripletsA.begin(), tripletsA.end());
}

//-----------------------------------------------------------------------------

void lump_matrix(SparseMatrix &D) {

    std::vector<Triplet> triplets;
    triplets.reserve(D.rows() * 20);
    for (int k = 0; k < D.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(D, k); it; ++it) {
            triplets.emplace_back(it.row(), it.row(), it.value());
        }
    }
    D.setFromTriplets(triplets.begin(), triplets.end());
}

//===================Gradient matrix Computation ==========================================================

Eigen::Vector3d gradient_hat_function(pmp::Point i,pmp::Point j, pmp::Point k) {
    pmp::Point base, side, grad;
    Eigen::Vector3d gradient;
    double area;
    area = triangle_area(i, j, k);
    side = i - j;
    base = k - j;
    grad = side - (dot(side, base) / norm(base)) * base / norm(base);
    if (area < eps) {
        gradient = Eigen::Vector3d(0, 0, 0);
    } else {
        grad = norm(base) * grad / norm(grad);
        gradient = Eigen::Vector3d(grad[0], grad[1], grad[2]) / (2.0 * area);
    }

    return gradient;
}

//-----------------------------------------------------------------------------
Eigen::Vector3d gradient_hat_function(Eigen::Vector3d i, Eigen::Vector3d j, Eigen::Vector3d k) {
    Eigen::Vector3d gradient, side, base, grad;
    double area = 0.5 * ((j - i).cross(k - i)).norm();
    side = i - j;
    base = k - j;
    grad = side - (side.dot(base) / base.norm()) * base / base.norm();
    if (area < eps) {
        gradient = Eigen::Vector3d(0, 0, 0);
    } else {
        grad = base.norm() * grad / grad.norm();
        gradient = grad / (2.0 * area);
    }
    return gradient;
}

//===================Minimization for Squared Area Point through vertex weights derivatives=============================

void setup_face_point_properties(pmp::SurfaceMesh &mesh, unsigned int min_point) {

    auto area_points = mesh.get_face_property<pmp::Point>("f:point");
    auto area_weights = mesh.get_face_property<Eigen::VectorXd>("f:weights");

    Eigen::MatrixXd T, P, PP;
    std::vector<Eigen::VectorXd> weights, testWeights, testWeightsPhil;
    Eigen::VectorXd w;
    Eigen::MatrixXd poly;
    Eigen::Vector3d p;
    std::vector<Eigen::Triplet<double>> trip;

    for (pmp::Face f: mesh.faces()) {
        const int n = mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (pmp::Vertex v : mesh.vertices(f)) {
            for (int h = 0; h < 3; h++) {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }
        if (min_point == Centroid) {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double) val;
        } else {
            // All other methods not relevant atm, so squared triangle area minimization is default
            find_area_minimizer_weights(poly, w);
        }

        Eigen::Vector3d min = poly.transpose() * w;
        pmp::Point point = pmp::Point(min(0), min(1), min(2));
        area_points[f] = point;
        area_weights[f] = w;
    }
}

//--------------------------------------------------------------------------------

void find_area_minimizer_weights(const Eigen::MatrixXd &poly,
                       Eigen::VectorXd &weights) {

    int val = poly.rows();
    Eigen::MatrixXd J(val, val);
    Eigen::VectorXd b(val);
    weights.resize(val);

    for (int i = 0; i < val; i++) {
        Eigen::Vector3d pk = poly.row(i);

        double Bk1_d2 = 0.0;
        double Bk1_d1 = 0.0;

        double Bk2_d0 = 0.0;
        double Bk2_d2 = 0.0;

        double Bk3_d0 = 0.0;
        double Bk3_d1 = 0.0;

        double CBk = 0.0;
        Eigen::Vector3d d = Eigen::MatrixXd::Zero(3, 1);

        for (int j = 0; j < val; j++) {
            Eigen::Vector3d pi = poly.row(j);
            Eigen::Vector3d pj = poly.row((j + 1) % val);
            d = pi - pj;


            double Bik1 = d(1) * pk(2) - d(2) * pk(1);
            double Bik2 = d(2) * pk(0) - d(0) * pk(2);
            double Bik3 = d(0) * pk(1) - d(1) * pk(0);

            double Ci1 = d(1) * pi(2) - d(2) * pi(1);
            double Ci2 = d(2) * pi(0) - d(0) * pi(2);
            double Ci3 = d(0) * pi(1) - d(1) * pi(0);

            Bk1_d1 += d(1) * Bik1;
            Bk1_d2 += d(2) * Bik1;

            Bk2_d0 += d(0) * Bik2;
            Bk2_d2 += d(2) * Bik2;

            Bk3_d0 += d(0) * Bik3;
            Bk3_d1 += d(1) * Bik3;

            CBk += Ci1 * Bik1 + Ci2 * Bik2 + Ci3 * Bik3;
        }
        for (int k = 0; k < val; k++) {
            Eigen::Vector3d xj = poly.row(k);
            J(i, k) = 0.5 * (xj(2) * Bk1_d1 - xj(1) * Bk1_d2 + xj(0) * Bk2_d2 -
                             xj(2) * Bk2_d0 + xj(1) * Bk3_d0 - xj(0) * Bk3_d1);
        }
        b(i) = 0.5 * CBk;

    }



//====================================================================

    Eigen::MatrixXd M(val + 1, val);
    M.block(0, 0, val, val) = 4 * J;
    M.block(val, 0, 1, val).setOnes();

//====================================================================

    Eigen::VectorXd b_(val + 1);
    b_.block(0, 0, val, 1) = 4 * b;
    b_(val) = 1.;
    weights = M.completeOrthogonalDecomposition().solve(b_).topRows(val);
}

//----------------Abs point ----------------------------------------------------------------

void absAreaJacobian(const Eigen::Matrix<double, -1, 3> &poly, const Eigen::Vector3d &a, Eigen::Vector3d &jac) {
    const int n = (int) poly.rows();
    jac.setZero();

    for (int i = 0; i < n; ++i) {
        Eigen::RowVector3d b = poly.row(i);
        Eigen::RowVector3d c = poly.row((i + 1) % n);

        const double expr0 = 2 *
                             pow(pow(b[0] * a[1] - c[0] * a[1] - a[0] * b[1] + c[0] * b[1] + a[0] * c[1] -
                                     b[0] * c[1],
                                     2) +
                                 pow(b[1] * a[2] - c[1] * a[2] - a[1] * b[2] + c[1] * b[2] + a[1] * c[2] -
                                     b[1] * c[2],
                                     2) +
                                 pow(c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2]), 2),
                                 0.5);

        assert(std::abs(expr0) > 1e-5);


        jac[0] += (2 * (-b[1] + c[1]) * (c[0] * (-a[1] + b[1]) + b[0] * (a[1] - c[1]) + a[0] * (-b[1] + c[1])) +
                   2 * (-b[2] + c[2]) * (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2]))) /
                  expr0;

        jac[1] += (2 * (b[0] - c[0]) * (c[0] * (-a[1] + b[1]) + b[0] * (a[1] - c[1]) + a[0] * (-b[1] + c[1])) +
                   2 * (-b[2] + c[2]) * (c[1] * (-a[2] + b[2]) + b[1] * (a[2] - c[2]) + a[1] * (-b[2] + c[2]))) /
                  expr0;

        jac[2] += (2 * (b[0] - c[0]) * (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2])) +
                   2 * (b[1] - c[1]) * (c[1] * (-a[2] + b[2]) + b[1] * (a[2] - c[2]) + a[1] * (-b[2] + c[2]))) /
                  expr0;
    }
}

void absAreaHessian(const Eigen::Matrix<double, -1, 3> &poly, const Eigen::Vector3d &a, Eigen::Matrix3d &hess) {
    using std::pow;

    const int n = (int) poly.rows();
    hess.setZero();
    for (int i = 0; i < n; ++i) {
        Eigen::RowVector3d b = poly.row(i);
        Eigen::RowVector3d c = poly.row((i + 1) % n);

        const double expr0 = pow(
                pow(b[0] * a[1] - c[0] * a[1] - a[0] * b[1] + c[0] * b[1] + a[0] * c[1] - b[0] * c[1], 2) +
                pow(b[1] * a[2] - c[1] * a[2] - a[1] * b[2] + c[1] * b[2] + a[1] * c[2] - b[1] * c[2], 2) +
                pow(c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2]), 2), 1.5);

        const double expr1 = pow(b[0] - c[0], 2) + pow(b[1] - c[1], 2) + pow(b[2] - c[2], 2);

        const double f = expr1 / expr0;

        hess(0, 0) += f * pow(b[1] * a[2] - c[1] * a[2] - a[1] * b[2] + c[1] * b[2] + a[1] * c[2] - b[1] * c[2], 2);

        hess(1, 0) += f * ((c[1] * (a[2] - b[2]) + a[1] * (b[2] - c[2]) + b[1] * (-a[2] + c[2])) *
                           (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2])));

        hess(2, 0) += f * ((c[0] * (a[1] - b[1]) + a[0] * (b[1] - c[1]) + b[0] * (-a[1] + c[1])) *
                           (c[1] * (a[2] - b[2]) + a[1] * (b[2] - c[2]) + b[1] * (-a[2] + c[2])));


        hess(0, 1) += f * ((c[1] * (a[2] - b[2]) + a[1] * (b[2] - c[2]) + b[1] * (-a[2] + c[2])) *
                           (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2])));

        hess(1, 1) += f * pow(c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2]), 2);

        hess(2, 1) += f * ((c[0] * (a[1] - b[1]) + a[0] * (b[1] - c[1]) + b[0] * (-a[1] + c[1])) *
                           (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2])));


        hess(0, 2) += f * ((c[0] * (a[1] - b[1]) + a[0] * (b[1] - c[1]) + b[0] * (-a[1] + c[1])) *
                           (c[1] * (a[2] - b[2]) + a[1] * (b[2] - c[2]) + b[1] * (-a[2] + c[2])));

        hess(1, 2) += f * ((c[0] * (a[1] - b[1]) + a[0] * (b[1] - c[1]) + b[0] * (-a[1] + c[1])) *
                           (c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2])));

        hess(2, 2) += f * pow(c[0] * (-a[1] + b[1]) + b[0] * (a[1] - c[1]) + a[0] * (-b[1] + c[1]), 2);
    }
}

double absTriangleArea(const Eigen::MatrixXd &poly, const Eigen::VectorXd &a) {
    using std::pow;

    const int n = (int) poly.rows();
    double area = 0;

    for (int i = 0; i < n; ++i) {
        Eigen::RowVector3d b = poly.row(i);
        Eigen::RowVector3d c = poly.row((i + 1) % n);

        area += pow(pow(b[0] * a[1] - c[0] * a[1] - a[0] * b[1] + c[0] * b[1] + a[0] * c[1] - b[0] * c[1], 2) +
                    pow(b[1] * a[2] - c[1] * a[2] - a[1] * b[2] + c[1] * b[2] + a[1] * c[2] - b[1] * c[2], 2) +
                    pow(c[0] * (-a[2] + b[2]) + b[0] * (a[2] - c[2]) + a[0] * (-b[2] + c[2]), 2), 0.5);
    }

    return area;
}


void optimizeAbsoluteTriangleArea(const Eigen::MatrixXd &poly, Eigen::Vector3d &x,
                                  const int maxIter, const int maxLineSearchIter) {
    // initial guess: mean of polygon
    x = poly.colwise().mean();
    double currArea = absTriangleArea(poly, x);

    // newton iterations
    for (int i = 0; i < maxIter; ++i) {
        Eigen::Matrix3d hess;
        Eigen::Vector3d jac;

        absAreaJacobian(poly, x, jac);
        absAreaHessian(poly, x, hess);

        const Eigen::Vector3d dir = hess.ldlt().solve(jac);

        // line search: half the step size at every iteration.
        double step = 1.;
        for (int j = 0; j < maxLineSearchIter; ++j) {
            Eigen::Vector3d xj = x - step * dir;
            double ar = absTriangleArea(poly, xj);

            if (ar < currArea) {
                currArea = ar;
                x = xj;
                step = 0;
                break;
            }

            step *= 0.5;
        }

        // could not improve the objective in search direction. Quit.
        if (step > 0) return;
    }
}

void find_weights_for_point(const Eigen::MatrixXd &poly,
                            const Eigen::Vector3d &point,
                            Eigen::VectorXd &weights) {
    const int n = (int) poly.rows();
    Eigen::MatrixXd M(4, n);
    M.block(0, 0, 3, n) = poly.transpose();
    M.row(3).setOnes();

    Eigen::VectorXd b(4);
    b.setZero();
    b.topRows(3) = point;
    b(3) = 1.;

    weights = M.completeOrthogonalDecomposition().solve(b);
}

Eigen::Vector3d triangle_circumcenter(const Eigen::MatrixXd &poly) {
//    Vector3f a,b,c // are the 3 pts of the tri
    Eigen::RowVector3d a = poly.row(0);
    Eigen::RowVector3d b = poly.row(1);
    Eigen::RowVector3d c = poly.row(2);

    Eigen::Vector3d ac = c - a;
    Eigen::Vector3d ab = b - a;
    Eigen::Vector3d bc = c - b;

    Eigen::Vector3d abXac = ab.cross(ac);

// this is the vector from a TO the circumsphere center
    Eigen::Vector3d toCircumsphereCenter =
            (abXac.cross(ab) * ac.squaredNorm() + ac.cross(abXac) * ab.squaredNorm()) / (2.0 * abXac.squaredNorm());

    double circumsphereRadius = toCircumsphereCenter.norm();


    // The 3 space coords of the circumsphere center then:
    Eigen::Vector3d ccs = a.transpose() + toCircumsphereCenter; // now this is the actual 3space location

    //Test: Distance ccs - point = radius
    if (abs((a.transpose() - ccs).norm() - circumsphereRadius) > 0.0001) {
        std::cout << "Point a not on circumsphere!" << std::endl;
    }
    if (abs((b.transpose() - ccs).norm() - circumsphereRadius) > 0.0001) {
        std::cout << "Point b not on circumsphere!" << std::endl;
    }
    if (abs((c.transpose() - ccs).norm() - circumsphereRadius) > 0.0001) {
        std::cout << "Point c not on circumsphere!" << std::endl;
    }
    //Test: bicenter orthogonal
    if (abs((0.5 * (a + b).transpose() - ccs).dot(ab)) > 0.00001)std::cout << "not orthogonal to ab" << std::endl;
    if (abs((0.5 * (a + c).transpose() - ccs).dot(ac)) > 0.00001)std::cout << "not orthogonal to ac" << std::endl;
    if (abs((0.5 * (c + b).transpose() - ccs).dot(bc)) > 0.00001)std::cout << "not orthogonal to bc" << std::endl;

    return ccs;
}

void barycentric_weights_triangle(const Eigen::MatrixXd &poly,
                                  const Eigen::Vector3d &point,
                                  Eigen::VectorXd &weights) {

    Eigen::RowVector3d a = poly.row(0);
    Eigen::RowVector3d b = poly.row(1);
    Eigen::RowVector3d c = poly.row(2);

    Eigen::Vector3d v0 = b - a, v1 = c - a, v2 = point - a.transpose();
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    weights.resize(3);
    weights(1) = (d11 * d20 - d01 * d21) / denom;
    weights(2) = (d00 * d21 - d01 * d20) / denom;
    weights(0) = 1.0f - weights(1) - weights(2);
}
