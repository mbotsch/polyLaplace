//=============================================================================

#include "diffgeo.h"
#include <Eigen/Dense>
#include "[AW11]Laplace.h"
#include <unsupported/Eigen/NonLinearOptimization>



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
        const int n = (int)mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (pmp::Vertex v : mesh.vertices(f)) {
            for (int h = 0; h < 3; h++) {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }
        if (min_point == Centroid) {
            int val = (int)poly.rows();
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

    int val = (int)poly.rows();
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



