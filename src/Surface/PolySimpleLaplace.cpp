
#include "PolySimpleLaplace.h"
#include "diffgeo.h"

//----------------------------------------------------------------------------------
enum InsertedPoint {
    Centroid = 0,
    AreaMinimizer = 2,
};

void setup_stiffness_matrix(SurfaceMesh &mesh,
                            Eigen::SparseMatrix<double> &S,
                            int minpoint) {

    const int nv = (int)mesh.n_vertices();

    Eigen::MatrixXd Si;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f: mesh.faces()) {
        int n = (int)mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v: mesh.vertices(f)) {
            for (int h = 0; h < 3; h++) {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        if (minpoint == Centroid) {
            int val = (int)poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double) val;
        } else {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min;

        min = poly.transpose() * w;
        localCotanMatrix(poly, min, w, Si);

        int j = 0;
        int k;
        for (Vertex v: mesh.vertices(f)) {
            k = 0;
            for (Vertex vv: mesh.vertices(f)) {

                trip.emplace_back(vv.idx(), v.idx(), Si(k, j));
                k++;
            }
            j++;
        }
    }

    S.resize(nv, nv);
    S.setFromTriplets(trip.begin(), trip.end());
    S *= -1.0;
}
//----------------------------------------------------------------------------------

void setup_mass_matrix(SurfaceMesh &mesh,
                       Eigen::SparseMatrix<double> &M, int minpoint) {
    const int nv = (int)mesh.n_vertices();

    Eigen::MatrixXd Mi;
    Eigen::VectorXd w;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;

    for (Face f: mesh.faces()) {
        const int n = (int)mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v: mesh.vertices(f)) {
            for (int h = 0; h < 3; h++) {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // setup polygon weights
        if (minpoint == Centroid) {
            int val = (int)poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double) val;
        } else {
            find_area_minimizer_weights(poly, w);
        }

        Eigen::Vector3d min = poly.transpose() * w;
        localMassMatrix(poly, min, w, Mi);

        int j = 0;
        int k;
        for (Vertex v: mesh.vertices(f)) {
            k = 0;
            for (Vertex vv: mesh.vertices(f)) {

                trip.emplace_back(vv.idx(), v.idx(), Mi(k, j));
                k++;
            }
            j++;
        }
    }
    M.resize(nv, nv);
    M.setFromTriplets(trip.begin(), trip.end());
}

//----------------------------------------------------------------------------------

void setup_gradient_matrix(SurfaceMesh &mesh,
                           Eigen::SparseMatrix<double> &G,
                           int minpoint) {

    const int nv = (int)mesh.n_vertices();

    Eigen::MatrixXd Gi;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;
    int nr_triangles = 0;
    int s = 0;
    for (Face f: mesh.faces()) {
        const int n = (int)mesh.valence(f);
        nr_triangles += n;
        poly.resize(n, 3);
        int row = 0;
        for (auto h: mesh.halfedges(f)) {
            Vertex v = mesh.from_vertex(h);
            for (int j = 0; j < 3; j++) {
                poly.row(row)(j) = mesh.position(v)[j];
            }
            row++;
        }
        // compute weights for the polygon
        if (minpoint == Centroid) {
            int val = (int)poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double) val;
        } else {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min;

        min = poly.transpose() * w;
        localGradientMatrix(poly, min, Gi);

        // Prolongation
        for (int j = 0; j < n; ++j)
            for (int i = 0; i < n; ++i)
                for (int k = 0; k < 3; k++)
                    Gi(3 * i + k, j) += w(j) * Gi(3 * i + k, n);

        int j = 0;
        int k;
        for (auto vv: mesh.vertices(f)) {
            k = 0;
            for (auto h: mesh.halfedges(f)) {
                Vertex v = mesh.from_vertex(h);
                for (int i = 0; i < 3; i++) {
                    trip.emplace_back(3 * s + i, v.idx(), Gi(3 * j + i, k));
                }
                k++;
            }
            j++;
            s++;
        }
    }

    G.resize(3 * nr_triangles, nv);
    G.setFromTriplets(trip.begin(), trip.end());
}
//----------------------------------------------------------------------------------

void setup_divergence_matrix(SurfaceMesh &mesh,
                             Eigen::SparseMatrix<double> &D,
                             int minpoint) {
    const int nv = (int)mesh.n_vertices();

    Eigen::MatrixXd Gi, Di;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;

    std::vector<Eigen::Triplet<double>> trip;
    int nr_triangles = 0;
    int s = 0;
    for (Face f: mesh.faces()) {

        const int n = (int)mesh.valence(f);
        nr_triangles += n;
        poly.resize(n, 3);
        int row = 0;

        for (auto h: mesh.halfedges(f)) {
            Vertex v = mesh.from_vertex(h);
            for (int j = 0; j < 3; j++) {
                poly.row(row)(j) = mesh.position(v)[j];
            }
            row++;
        }

        // compute weights for the polygon
        if (minpoint == Centroid) {
            int val = (int)poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double) val;
        } else {
            find_area_minimizer_weights(poly, w);
        }

        Eigen::Vector3d min = poly.transpose() * w;
        localGradientMatrix(poly, min, Gi);
        Di = -Gi.transpose();

        //Triangle Area Diagonal matrix
        Eigen::MatrixXd Ai;
        Ai.resize(3 * n, 3 * n);
        Ai.setZero();
        for (int i = 0; i < n; ++i) {
            const int i1 = (i + 1) % n;

            Eigen::Vector3d p0 = poly.row(i);
            Eigen::Vector3d p1 = poly.row(i1);

            double area = 0.5 * ((p0 - min).cross(p1 - min)).norm();
            for (int k = 0; k < 3; k++) {
                Ai(3 * i + k, 3 * i + k) = area;
            }
        }
        Di *= Ai;

        // prolongation

        for (int j = 0; j < n; ++j) {
            for (int i = 0; i < n; ++i) {
                for (int k = 0; k < 3; k++)
                    Di(i, 3 * j + k) += w(i) * Di(n, 3 * j + k);
            }
        }

        int j = 0;
        int k;
        for (auto vv: mesh.vertices(f)) {
            k = 0;
            for (auto h: mesh.halfedges(f)) {
                Vertex v = mesh.from_vertex(h);
                for (int i = 0; i < 3; i++) {
                    trip.emplace_back(v.idx(), 3 * s + i, Di(k, 3 * j + i));
                }
                k++;
            }
            j++;
            s++;
        }
    }

    D.resize(nv, 3 * nr_triangles);
    D.setFromTriplets(trip.begin(), trip.end());
}
//----------------------------------------------------------------------------------

void localCotanMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min,
                      Eigen::VectorXd &w, Eigen::MatrixXd &L) {
    const int n = (int) poly.rows();

    L.resize(n, n);
    L.setZero();

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i) {
        const int i1 = (i + 1) % n;

        l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
        l2[0] = (poly.row(i1) - min.transpose()).squaredNorm();
        l2[1] = (poly.row(i) - min.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.5 * sqrt(arg);
        if (area > 1e-9) {
            l[0] = 0.25 * (l2[1] + l2[2] - l2[0]) / area;
            l[1] = 0.25 * (l2[2] + l2[0] - l2[1]) / area;
            l[2] = 0.25 * (l2[0] + l2[1] - l2[2]) / area;

            L(i1, i1) += l[0];
            L(i, i) += l[1];
            L(i1, i) -= l[2];
            L(i, i1) -= l[2];
            L(i, i) += l[2];
            L(i1, i1) += l[2];

            ln(i1) -= l[0];
            ln(i) -= l[1];
            ln(n) += l[0] + l[1];
        }
    }

    // Sandwiching
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)

            L(i, j) += w(i) * ln(j) + w(j) * ln(i) + w(i) * w(j) * ln(n);
}

//----------------------------------------------------------------------------------

void localMassMatrix(const Eigen::MatrixXd &poly, const Eigen::Vector3d &min,
                     Eigen::VectorXd &w, Eigen::MatrixXd &M) {
    const int n = (int) poly.rows();

    M.resize(n, n);

    M.setZero();

    Eigen::VectorXd ln(n + 1);
    ln.setZero();

    double l[3], l2[3];

    for (int i = 0; i < n; ++i) {

        const int i1 = (i + 1) % n;

        l2[2] = (poly.row(i) - poly.row(i1)).squaredNorm();
        l2[0] = (poly.row(i1) - min.transpose()).squaredNorm();
        l2[1] = (poly.row(i) - min.transpose()).squaredNorm();

        l[0] = sqrt(l2[0]);
        l[1] = sqrt(l2[1]);
        l[2] = sqrt(l2[2]);

        const double arg = (l[0] + (l[1] + l[2])) * (l[2] - (l[0] - l[1])) *
                           (l[2] + (l[0] - l[1])) * (l[0] + (l[1] - l[2]));
        const double area = 0.25 * sqrt(arg);

        l[0] = 1.0 / 6.0 * area;
        l[1] = 1.0 / 12.0 * area;

        M(i1, i1) += 1.0 / 6.0 * area;
        M(i, i) += 1.0 / 6.0 * area;
        M(i1, i) += 1.0 / 12.0 * area;
        M(i, i1) += 1.0 / 12.0 * area;

        ln(i1) += l[1];
        ln(i) += l[1];
        ln(n) += l[0];
    }
    // Sandwiching
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            M(i, j) += w(i) * ln(j) + w(j) * ln(i) + w(i) * w(j) * ln(n);
}

//----------------------------------------------------------------------------------

void localGradientMatrix(const Eigen::MatrixXd &poly,
                         const Eigen::Vector3d &min,
                         Eigen::MatrixXd &G) {

    const int n = (int) poly.rows();

    G.resize(3 * n, n + 1);

    G.setZero();

    Eigen::Vector3d gradient_p, gradient_p0, gradient_p1, p, p0, p1;

    p = min;
    for (int i = 0; i < n; ++i) {

        const int i1 = (i + 1) % n;

        p0 = poly.row(i);
        p1 = poly.row(i1);

        gradient_p = gradient_hat_function(p, p0, p1);
        gradient_p0 = gradient_hat_function(p0, p1, p);
        gradient_p1 = gradient_hat_function(p1, p, p0);

        for (int j = 0; j < 3; j++) {
            G(3 * i + j, n) = gradient_p(j);
            G(3 * i + j, i) = gradient_p0(j);
            G(3 * i + j, i1) = gradient_p1(j);
        }
    }
}
//----------------------------------------------------------------------------------

Eigen::Vector3d gradient_hat_function(Eigen::Vector3d &i, Eigen::Vector3d &j, Eigen::Vector3d &k) {
    Eigen::Vector3d gradient, side, base, grad;
    double area = 0.5 * ((j - i).cross(k - i)).norm();
    side = i - j;
    base = k - j;
    const double eps = 1e-10;
    grad = side - (side.dot(base) / base.norm()) * base / base.norm();
    if (area < eps) {
        gradient = Eigen::Vector3d(0, 0, 0);
    } else {
        grad = base.norm() * grad / grad.norm();
        gradient = grad / (2.0 * area);
    }
    return gradient;
}
