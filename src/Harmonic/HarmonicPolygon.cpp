#include "HarmonicPolygon.h"
#include "PolygonSampler.h"
#include <fstream>
#include "GaussQuadrature.h"
#include <pmp/surface_mesh.h>
#include "../Surface/diffgeo.h"

Eigen::Vector3d HarmonicPolygon::getScaledNormal() const
{
    return area * a;
}

double HarmonicPolygon::getD() const
{
    return dPlane;
}

const Eigen::Vector3d& HarmonicPolygon::getNormal() const
{
    return a;
}

Eigen::Vector2d HarmonicPolygon::projectToPlane(const Eigen::Vector3d& p)
{
    return P * (p - mean) / scale;
}

Eigen::Vector3d HarmonicPolygon::unproject(const Eigen::Vector2d& p)
{
    return mean + scale * (P.transpose() * p);
}

std::vector<Eigen::Vector3d> HarmonicPolygon::unproject(
    const std::vector<Eigen::Vector2d>& p)
{
    std::vector<Eigen::Vector3d> ret(p.size());
    for (int i = 0; i < (int)p.size(); ++i)
        ret[i] = unproject(p[i]);
    return ret;
}

std::vector<Eigen::Vector2d> HarmonicPolygon::generateRandomSamples(
    const int nSamples) const
{
    return PolygonSampler(poly2d, area2d, nSamples).polySamples;
}

double HarmonicPolygon::quadrature(
    const std::function<double(double, double)>& g)
{
    double ret = .0;

    for (auto& t : quadratureTriangles)
        ret += gaussQuadratureTriangle(t, g, 4);

    return ret;
}

void HarmonicPolygon::initQuadrature()
{
    const int np = (int)poly2d.rows();

    //------- pmp version---------------------------
    pmp::SurfaceMesh polygon;
    Eigen::MatrixXd poly;
    poly.resize(np, 3);

    std::vector<pmp::Vertex> face_vertices;
    for (int i = 0; i < np; ++i)
    {
        poly.row(i) = Eigen::Vector3d(poly2d(i, 0), poly2d(i, 1), 0.0);
    }

    Eigen::VectorXd w;

    find_area_minimizer_weights(poly, w);
    Eigen::Vector3d point = poly.transpose() * w;
    quadratureTriangles.resize(np, Eigen::MatrixXd(3, 2));
    for (int i = 0; i < np; ++i)
    {
        quadratureTriangles[i].row(0) =
            Eigen::Vector2d(poly2d(i, 0), poly2d(i, 1));
        quadratureTriangles[i].row(1) = Eigen::Vector2d(point(0), point(1));
        quadratureTriangles[i].row(2) = Eigen::Vector2d(
            poly2d((i + 1) % (np), 0), poly2d((i + 1) % (np), 1));
    }

    //--------- libigl version ------------------------------
    //Flags meaning:
    // -a Imposes a maximum triangle area constraint
    // -q Quality mesh generation with no angles smaller than 20 degrees
    // -Q Quiet

    //    Eigen::MatrixXi E(np, 2);
    //
    //    for (int i = 0; i < np; ++i) {
    //        E(i, 0) = i;
    //        E(i, 1) = (i + 1) % np;
    //    }
    //
    //    Eigen::MatrixXd V2;
    //    Eigen::MatrixXi F2;
    //
    //    igl::triangle::triangulate(poly2d, E, Eigen::MatrixXi(0, 0), "-a0.05-q-Q", V2, F2);
    //
    //    quadratureTriangles.resize(F2.rows(), Eigen::MatrixXd(3, 2));
    //
    //    for (int i = 0; i < F2.rows(); ++i) {
    //        for (int j = 0; j < 3; ++j) {
    //            quadratureTriangles[i].row(j) = V2.row(F2(i, j));
    //        }
    //    }
}

double HarmonicPolygon::evaluate(const int id, const Eigen::Vector2d& p)
{
    double val = p(0) * coefficients(id, nnKernels) +
                 p(1) * coefficients(id, nnKernels + 1) +
                 coefficients(id, nnKernels + 2);

    for (int i = 0; i < nnKernels; ++i)
    {
        val +=
            coefficients(id, i) * log((kernels.row(i) - p.transpose()).norm());
    }

    return val;
}

Eigen::Vector2d HarmonicPolygon::evaluateGrad(const int id,
                                              const Eigen::Vector2d& p)
{
    Eigen::Vector2d ret(coefficients(id, nnKernels),
                        coefficients(id, nnKernels + 1));

    for (int i = 0; i < nnKernels; ++i)
    {
        ret += coefficients(id, i) /
               (kernels.row(i) - p.transpose()).squaredNorm() *
               (p.transpose() - kernels.row(i));
    }

    return ret;
}

void HarmonicPolygon::stiffnessMatrix(Eigen::MatrixXd& K)
{
    const int nv = (int)poly2d.rows();
    K.resize(nv, nv);

    for (int i = 0; i < nv; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double val = quadrature([&](const double x, const double y) {
                return evaluateGrad(i, Eigen::Vector2d(x, y))
                    .dot(evaluateGrad(j, Eigen::Vector2d(x, y)));
            });

            K(j, i) = K(i, j) = val;
        }
    }
}

void HarmonicPolygon::massMatrix(Eigen::MatrixXd& M)
{
    const int nv = (int)poly2d.rows();
    M.resize(nv, nv);

    for (int i = 0; i < nv; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            double val =
                scale * scale * quadrature([&](const double x, const double y) {
                    return evaluate(i, Eigen::Vector2d(x, y)) *
                           evaluate(j, Eigen::Vector2d(x, y));
                });

            M(j, i) = M(i, j) = val;
        }
    }
}

std::vector<double> HarmonicPolygon::evaluateAtPoints(
    const int id, const std::vector<Eigen::Vector2d>& points)
{
    if (id < 0 || id > coefficients.rows())
    {
        return std::vector<double>(points.size(), .0);
    }

    std::vector<double> ret;

    ret.reserve(points.size());

    for (auto& p : points)
    {
        double val = p(0) * coefficients(id, nnKernels) +
                     p(1) * coefficients(id, nnKernels + 1) +
                     coefficients(id, nnKernels + 2);

        for (int i = 0; i < nnKernels; ++i)
        {
            val += coefficients(id, i) *
                   log((kernels.row(i) - p.transpose()).norm());
        }

        ret.push_back(val);
    }

    return ret;
}

void HarmonicPolygon::dump2d(const std::string& fname)
{
    std::ofstream file(fname);
    file << n << " " << nKernels << " " << nProbes << std::endl;
    file << poly2d << std::endl;
    file << kernels << std::endl;
    file << coefficients << std::endl;
    file << probes;
    file.close();
}

HarmonicPolygon::HarmonicPolygon(const Eigen::MatrixXd& pts_)
{
    n = (int)pts_.rows();
    nnProbes = nProbes * n;
    nnKernels = nKernels * n;

    ///////////////////////////////////////////////
    // shift and rescale points

    mean = pts_.colwise().mean();
    Eigen::MatrixXd pts = pts_.rowwise() - mean.transpose();
    scale = pts.rowwise().norm().maxCoeff();
    pts /= scale;

    ///////////////////////////////////////////////
    // rotate to plane, normal vector given by vector area

    a.setZero();

    for (int i = 0; i < n; ++i)
    {
        a += Eigen::Vector3d(pts.row(i))
                 .cross(Eigen::Vector3d(pts.row((i + 1) % n)));
    }

    area = a.norm();

    a /= area;
    dPlane = a.dot(Eigen::Vector3d(pts_.row(0)));

    area *= 0.5;
    area2d = area;

    area *=
        scale * scale; // we want the original vector area: account for scaling

    ///////////////////////////////////////////////
    // span space orthogonal to 'a'
    P.resize(2, 3);
    P.setRandom();

    P.row(0) -= a * a.dot(P.row(0));
    P.row(0).normalize();

    P.row(1) -= a * a.dot(P.row(1));
    P.row(1) -= P.row(0) * P.row(0).dot(P.row(1));
    P.row(1).normalize();

    ///////////////////////////////////////////////
    // project pts to 2d
    poly2d = (P * pts.transpose()).transpose();

    ///////////////////////////////////////////////
    // setup kernels and probe matrix
    std::vector<Eigen::Triplet<double>> tripP;
    kernels.resize(n * nKernels, 2);

    for (int i = 0; i < n; ++i)
    {
        int i2 = i + 1 == n ? 0 : i + 1;

        Eigen::Vector2d normal(poly2d(i, 1) - poly2d(i2, 1),
                               poly2d(i2, 0) - poly2d(i, 0));
        normal.normalize();

        for (int j = 0; j < nKernels; ++j)
        {
            const double w = j / (double)nKernels;
            kernels.row(i * nKernels + j) = (1. - w) * poly2d.row(i) +
                                            w * poly2d.row(i2) +
                                            (eps)*normal.transpose();
        }

        for (int j = 0; j < nProbes; ++j)
        {
            const double w = j / (double)nProbes;
            tripP.emplace_back(i * nProbes + j, i, 1. - w);
            tripP.emplace_back(i * nProbes + j, i2, w);
        }
    }

    Pr.resize(n * nProbes, n);
    Pr.setFromTriplets(tripP.begin(), tripP.end());

    probes = Pr * poly2d;

    ///////////////////////////////////////////////
    // compute basis coefficients

    Eigen::MatrixXd A(nnProbes, nnKernels + 3);

    for (int i = 0; i < nnProbes; ++i)
    {
        for (int j = 0; j < nnKernels; ++j)
        {
            A(i, j) = log((kernels.row(j) - probes.row(i)).norm());
        }

        A(i, nnKernels + 0) = probes(i, 0);
        A(i, nnKernels + 1) = probes(i, 1);
        A(i, nnKernels + 2) = 1.;
    }

    // least squares solve: one for each node
    Eigen::MatrixXd ATA = A.transpose() * A;
    Eigen::LDLT<Eigen::MatrixXd> chol(ATA);

    coefficients =
        chol.solve(A.transpose() * (Pr * Eigen::MatrixXd::Identity(n, n)))
            .transpose();

    ///////////////////////////////////////////////
    // triangulate polygon for quadrature

    initQuadrature();
}

HarmonicPolygon::HarmonicPolygon(const Eigen::MatrixXd& pts_, int nKernel_,
                                 int nProbes_)
{
    nKernels = nKernel_;
    nProbes = nProbes_;
    n = (int)pts_.rows();
    nnProbes = nProbes * n;
    nnKernels = nKernels * n;

    ///////////////////////////////////////////////
    // shift and rescale points

    mean = pts_.colwise().mean();
    Eigen::MatrixXd pts = pts_.rowwise() - mean.transpose();
    scale = pts.rowwise().norm().maxCoeff();
    pts /= scale;

    ///////////////////////////////////////////////
    // rotate to plane, normal vector given by vector area

    a.setZero();

    for (int i = 0; i < n; ++i)
    {
        a += Eigen::Vector3d(pts.row(i))
                 .cross(Eigen::Vector3d(pts.row((i + 1) % n)));
    }

    area = a.norm();

    a /= area;
    dPlane = a.dot(Eigen::Vector3d(pts_.row(0)));

    area *= 0.5;
    area2d = area;

    area *=
        scale * scale; // we want the original vector area: account for scaling

    ///////////////////////////////////////////////
    // span space orthogonal to 'a'
    P.resize(2, 3);
    P.setRandom();

    P.row(0) -= a * a.dot(P.row(0));
    P.row(0).normalize();

    P.row(1) -= a * a.dot(P.row(1));
    P.row(1) -= P.row(0) * P.row(0).dot(P.row(1));
    P.row(1).normalize();

    ///////////////////////////////////////////////
    // project pts to 2d
    poly2d = (P * pts.transpose()).transpose();

    ///////////////////////////////////////////////
    // setup kernels and probe matrix
    std::vector<Eigen::Triplet<double>> tripP;
    kernels.resize(n * nKernels, 2);

    for (int i = 0; i < n; ++i)
    {
        int i2 = i + 1 == n ? 0 : i + 1;

        Eigen::Vector2d normal(poly2d(i, 1) - poly2d(i2, 1),
                               poly2d(i2, 0) - poly2d(i, 0));
        normal.normalize();

        for (int j = 0; j < nKernels; ++j)
        {
            const double w = j / (double)nKernels;
            kernels.row(i * nKernels + j) = (1. - w) * poly2d.row(i) +
                                            w * poly2d.row(i2) +
                                            (eps)*normal.transpose();
        }

        for (int j = 0; j < nProbes; ++j)
        {
            const double w = j / (double)nProbes;
            tripP.emplace_back(i * nProbes + j, i, 1. - w);
            tripP.emplace_back(i * nProbes + j, i2, w);
        }
    }

    Pr.resize(n * nProbes, n);
    Pr.setFromTriplets(tripP.begin(), tripP.end());

    probes = Pr * poly2d;

    ///////////////////////////////////////////////
    // compute basis coefficients

    Eigen::MatrixXd A(nnProbes, nnKernels + 3);

    for (int i = 0; i < nnProbes; ++i)
    {
        for (int j = 0; j < nnKernels; ++j)
        {
            A(i, j) = log((kernels.row(j) - probes.row(i)).norm());
        }

        A(i, nnKernels + 0) = probes(i, 0);
        A(i, nnKernels + 1) = probes(i, 1);
        A(i, nnKernels + 2) = 1.;
    }

    // least squares solve: one for each node
    Eigen::MatrixXd ATA = A.transpose() * A;
    Eigen::LDLT<Eigen::MatrixXd> chol(ATA);

    coefficients =
        chol.solve(A.transpose() * (Pr * Eigen::MatrixXd::Identity(n, n)))
            .transpose();

    ///////////////////////////////////////////////
    // triangulate polygon for quadrature

    initQuadrature();
}