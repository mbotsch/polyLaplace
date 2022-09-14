//=============================================================================
#include "[AW11]Laplace.h"
//=============================================================================

using SurfaceMesh = pmp::SurfaceMesh;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

float poly_laplace_lambda_ = 2.0;

//=============================================================================

void setup_E_and_B_perFace(pmp::SurfaceMesh &mesh, pmp::Face f,
                           Eigen::MatrixXd &E, Eigen::MatrixXd &B) {
    const int n = (int) mesh.valence(f);
    Eigen::MatrixXd X(n, 3);

    int i = 0;
    for (auto v: mesh.vertices(f)) {
        const auto p = mesh.position(v);
        X(i, 0) = p[0];
        X(i, 1) = p[1];
        X(i, 2) = p[2];
        ++i;
    }
    // compute pre-laplacian for face i
    E.resizeLike(X);
    B.resizeLike(X);

    for (i = 0; i < n; ++i) {
        E.row(i) = X.row((i + 1) % n) - X.row(i);
        B.row(i) = 0.5 * (X.row((i + 1) % n) + X.row(i));
    }
}
//-----------------------------------------------------------------------------

void setup_poly_gradient_operator(pmp::SurfaceMesh &mesh,
                                  Eigen::SparseMatrix<double> &G) {
    std::vector<Eigen::Triplet<double>> triplets;
    int cnt = 0;

    for (auto f: mesh.faces()) {
        const unsigned int n = mesh.valence(f);
        Eigen::VectorXi F(n);

        int idx = 0;
        for (auto v: mesh.vertices(f)) {
            F(idx) = (int)v.idx();
            ++idx;
        }

        for (unsigned int i = 0; i < n; ++i) {
            triplets.emplace_back(cnt, F(i), -1);
            triplets.emplace_back(cnt, F((i + 1) % n), 1);
            ++cnt;
        }
    }

    G.resize(cnt, (int)mesh.n_vertices());
    G.setFromTriplets(triplets.begin(), triplets.end());


}
//-----------------------------------------------------------------------------

void setup_poly_divergence_operator(pmp::SurfaceMesh &mesh,
                                    Eigen::SparseMatrix<double> &D) {
    double lambda = poly_laplace_lambda_;

    std::vector<Eigen::Triplet<double>> triplets;
    int colCnt = 0;

    for (auto f: mesh.faces()) {
        const int n = (int)mesh.valence(f);
        Eigen::VectorXi F(n);
        Eigen::MatrixXd E, B;
        setup_E_and_B_perFace(mesh, f, E, B);
        int idx = 0;
        for (auto v: mesh.vertices(f)) {
            F(idx) = (int)v.idx();
            ++idx;
        }
        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));

        const double area = af.norm();
        af /= area;

        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();

        for ( int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }

        Eigen::MatrixXd Lf;

        if (lambda > 0.0) {
            // fill kernel
            E -= E * af * af.transpose();

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(E.transpose(),
                                                  Eigen::ComputeFullV);
            const Eigen::MatrixXd C = svd.matrixV().rightCols(n - 2);

            assert((E.transpose() * C).norm() < 1e-10);

            // assemble face divergence
            Lf = d.transpose() *
                 ((B * B.transpose()) / area + lambda * C * C.transpose());
        } else
            Lf = d.transpose() * ((B * B.transpose()) / area);

        // add local laplacian to global matrix entries
        for ( int k = 0; k < n; ++k) {
            for ( int l = 0; l < n; ++l) {
                triplets.emplace_back(F(k), colCnt + l, -2.0 * Lf(k, l));
            }
        }

        colCnt += n;
    }

    D.resize((int)mesh.n_vertices(), colCnt);
    D.setFromTriplets(triplets.begin(), triplets.end());
}
//-----------------------------------------------------------------------------

void normalize_poly_gradients(pmp::SurfaceMesh &mesh, Eigen::VectorXd &g,
                              const Eigen::VectorXd &h) {
    double lambda = poly_laplace_lambda_;
    assert(h.rows() == (int)mesh.n_vertices());
    //  assert(mesh.n_halfedges() == g.rows());

    int hedgeCnt = 0;

    for (auto f: mesh.faces()) {
        const unsigned int n = mesh.valence(f);
        Eigen::VectorXi F(n);
        Eigen::MatrixXd E, B;
        setup_E_and_B_perFace(mesh, f, E, B);
        int idx = 0;
        for (auto v: mesh.vertices(f)) {
            F(idx) = (int)v.idx();
            ++idx;
        }
        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));

        const double area = af.norm();
        af /= area;

        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();

        for (unsigned int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }

        Eigen::MatrixXd Lf;

        if (lambda >0.0) {
            // fill kernel
            E -= E * af * af.transpose();

            Eigen::JacobiSVD<Eigen::MatrixXd> svd(E.transpose(),
                                                  Eigen::ComputeFullV);
            const Eigen::MatrixXd C = svd.matrixV().rightCols(n - 2);

            assert((E.transpose() * C).norm() < 1e-10);

            // assemble face laplacian
            Lf = d.transpose() *
                 ((B * B.transpose()) / area + lambda * C * C.transpose()) * d;
        } else
            Lf = d.transpose() * ((B * B.transpose()) / area) * d;

        // get heat values for face f

        Eigen::VectorXd hi(n);
        for (unsigned int i = 0; i < n; ++i)
            hi(i) = h(F(i));
        const double factor = sqrt(hi.transpose().dot(Lf * hi) / area);

        for (unsigned int i = 0; i < n; ++i)
            g(hedgeCnt++) /= factor;
    }

    assert(g.rows() == hedgeCnt);
}
//-----------------------------------------------------------------------------

void setup_poly_Laplace_matrix(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &L) {

    double lambda = poly_laplace_lambda_;
    std::vector<Eigen::Triplet<double>> triplets;
    //    Eigen::MatrixXd E, B;
    for (auto f: mesh.faces()) {
        const unsigned int n = mesh.valence(f);
        Eigen::VectorXi F(n);
        Eigen::MatrixXd E, B;
        setup_E_and_B_perFace(mesh, f, E, B);
        int idx = 0;
        for (auto v: mesh.vertices(f)) {
            F(idx) = (int)v.idx();
            ++idx;
        }
        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));

        const double area = af.norm();
        af /= area;

        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();

        for (unsigned int i = 0; i < n; ++i) {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }

        Eigen::MatrixXd Lf;

        if (lambda >0) {
#if 0 // Philipp's verion to flatten matrix E
            E -= E * af * af.transpose();
#else // Marc's version
            Eigen::MatrixXd X(3, n);
            int i = 0;
            for (auto v: mesh.vertices(f)) {
                X(0, i) = mesh.position(v)[0];
                X(1, i) = mesh.position(v)[1];
                X(2, i) = mesh.position(v)[2];
                ++i;
            }

            Eigen::MatrixXd flatX = X - af * af.transpose() * X;

            Eigen::MatrixXd flatE(3, n);
            for (unsigned int j = 0; j < n; ++j)
                flatE.col(j) = flatX.col((j + 1) % n) - flatX.col(j);

            E = flatE.transpose();
#endif

#if 0 // Philipp's version to compute the kernel of E^T
            Eigen::JacobiSVD<Eigen::MatrixXd> svd(E.transpose(),
                                                  Eigen::ComputeFullV);
            const Eigen::MatrixXd C = svd.matrixV().rightCols(n - 2);
#else // Marc's version
            Eigen::MatrixXd CLU =
                    Eigen::FullPivLU<Eigen::MatrixXd>(E.transpose()).kernel();
            const Eigen::MatrixXd C =
                    Eigen::JacobiSVD<Eigen::MatrixXd>(
                            CLU, Eigen::ComputeThinU | Eigen::ComputeThinV)
                            .matrixU();

            // uncomment to see that Marc's flattening leads to round-off errors
            // that in turn lead to a wrong estimation of the kernel dimension
            //if (CLU.cols() != n-2)
            //std::cout << "kernel dimension " << CLU.cols() << " != " << (n-2) << std::endl;
#endif

            if ((E.transpose() * C).norm() > 1e-10)
                std::cerr << "Should not happen\n";

            // assemble face laplacian
            Lf = d.transpose() *
                 ((B * B.transpose()) / area + lambda * C * C.transpose()) * d;
        } else
            Lf = d.transpose() * ((B * B.transpose()) / area) * d;

        // add local laplacian to global matrix entries
        for (unsigned int k = 0; k < n; ++k) {
            for (unsigned int l = 0; l < n; ++l) {
                triplets.emplace_back(F(k), F(l), -2.0 * Lf(k, l));
            }
        }
    }

    L.resize((int)mesh.n_vertices(), (int)mesh.n_vertices());
    L.setFromTriplets(triplets.begin(), triplets.end());
}

//-----------------------------------------------------------------------------

void setup_poly_mass_matrix(pmp::SurfaceMesh &mesh,
                            Eigen::SparseMatrix<double> &M) {
    M.resize((int)mesh.n_vertices(), (int)mesh.n_vertices());

    std::vector<Triplet> tripletsM;
    double sum = 0.0;
    Eigen::MatrixXd E, B;
    for (auto v: mesh.vertices()) {
        for (auto f: mesh.faces(v)) {
            const unsigned int n = mesh.valence(f);
            setup_E_and_B_perFace(mesh, f, E, B);
            // compute vector area
            const Eigen::Matrix3d Af = E.transpose() * B;
            Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
            const double area = af.norm();
            sum += area / n;
        }
        tripletsM.emplace_back(v.idx(), v.idx(), sum);
        sum = 0.0;
    }
    // build sparse matrix from triplets
    M.setFromTriplets(tripletsM.begin(), tripletsM.end());
}


