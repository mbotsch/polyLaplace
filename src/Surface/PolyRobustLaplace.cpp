#include "PolyRobustLaplace.h"
#include <pmp/algorithms/normals.h>
#include <random>
#include "diffgeo.h"

#include <Eigen/IterativeLinearSolvers>


void find_trace_minimizer_weights(const Eigen::MatrixXd& polyVerts, Eigen::VectorXd &weights, bool use_fallback)
{

    Eigen::Vector3d fn;
    Eigen::MatrixXd projectedVerts;
    Eigen::MatrixXd projMatrix;
    computeMaxAreaPolygonProjection(polyVerts, fn, projectedVerts, projMatrix);

    int val = (int)polyVerts.rows();

    Eigen::VectorXd virtualVertex;
    getLeastSquaresPoint(polyVerts, virtualVertex);

    virtualVertex = projMatrix.transpose() * (virtualVertex - fn * virtualVertex.dot(fn)); // use projected as init
    double convergence_eps = 1e-10;

    double oldNd = 0;
    for (int i = 0; i < 1000; ++i) {
        auto [e, g, H_proj] = getEnergyAndGradientAndHessian(virtualVertex, projectedVerts);
        project_positive_definite(H_proj);
        Eigen::Vector2d d = - H_proj.inverse() * g;

        double nd = -0.5 * d.dot(g);

        if (nd < convergence_eps || oldNd == nd) {
            break;
        }
        virtualVertex = line_search(virtualVertex, d, e, g, projectedVerts);
        if (!(i % 10))
            oldNd = nd;

    }

    Eigen::Vector2d vV(virtualVertex);
    weights.resize(val);
    computeHarmonicWeights(projectedVerts, weights, vV);

    if(use_fallback)
    {
        virtualVertex = polyVerts.transpose() * weights;

        Eigen::VectorXd area_weights;
        find_area_minimizer_weights(polyVerts, area_weights);
        Eigen::Vector3d area_vVerts = polyVerts.transpose() * area_weights;

        double old_trace = getStiffnessTrace(area_vVerts, polyVerts, area_weights);

        double new_trace = getStiffnessTrace(virtualVertex, polyVerts, weights);

        if (new_trace > 1e-10 + old_trace) {
            weights = area_weights;
        }
    }
}



/**
 * Check if matrix is diagonally dominant and has positive diagonal entries.
 * This is a sufficient condition for positive-definiteness
 * and can be used as an early out to avoid eigen decomposition.
 */
bool positive_diagonally_dominant(
    Eigen::Matrix2d& H,
    double _eps)
{
    for (Eigen::Index i = 0; i < H.rows(); ++i)
    {
        double off_diag_abs_sum = 0.0;
        for(Eigen::Index j = 0; j < H.cols(); ++j)
        {
            if (i != j)
                off_diag_abs_sum += std::abs(H(i, j));
        }

        if (H(i, i) < off_diag_abs_sum + _eps)
            return false;
    }

    return true;
}

/**
 * Project symmetric matrix to positive-definite matrix
 * via eigen decomposition.
 */
void project_positive_definite(
    Eigen::Matrix2d& H,
    double _eigenvalue_eps)
{

    // Early out if sufficient condition is fulfilled
    if (positive_diagonally_dominant(H, _eigenvalue_eps))
        return;

    // Compute eigen-decomposition (of symmetric matrix)
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig;
    eig.computeDirect(H);
    Eigen::Matrix2d D = eig.eigenvalues().asDiagonal();

    // Clamp all eigenvalues to eps
    bool all_positive = true;
    for (Eigen::Index i = 0; i < H.rows(); ++i)
    {
        if (D(i, i) < _eigenvalue_eps)
        {
            D(i, i) = _eigenvalue_eps;
            all_positive = false;
        }
    }

    // Do nothing if all eigenvalues were already at least eps
    if (all_positive)
        return;

    // Re-assemble matrix using clamped eigenvalues
    H = eig.eigenvectors() * D * eig.eigenvectors().transpose();
}

double getEnergy(const Eigen::Vector2d& vVert,
                                  const Eigen::MatrixXd& poly)
{
    int n = (int)poly.rows();
    double sum = 0;

    double lenIV, lenIR, lenVR;

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector2d lV = poly.row(i);
        Eigen::Vector2d rV = poly.row((i + 1) % n);

        Eigen::Vector2d e1 = lV - vVert;
        Eigen::Vector2d e2 = rV - vVert;
        double area = (e1[0] * e2[1] - e2[0] * e1[1])/2.0;

        if (area < 0)
        { // flipped triangle
            return penalty_;
        }
        Eigen::Vector2d e = rV - lV;

        lenIV = e1.norm();
        lenVR = e2.norm();
        lenIR = e.norm();
        if ((lenIV+lenVR) / lenIR - 1 < 1e-6)
        {
            return penalty_;
        }

        sum += (lenIV*lenIV + lenVR*lenVR + lenIR*lenIR) / (4 * area);
    }
    return sum;
}

std::tuple<double, Eigen::Vector2d, Eigen::Matrix2d> getEnergyAndGradientAndHessian(const Eigen::Vector2d& vVert,
                                                                                                     const Eigen::MatrixXd& poly)
{
    int n = (int)poly.rows();
    Eigen::Vector2d gradientSum = Eigen::Vector2d::Zero();
    Eigen::Matrix2d hessianSum = Eigen::Matrix2d::Zero();
    double energySum = 0;

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector2d lV = poly.row(i);
        Eigen::Vector2d rV = poly.row((i + 1) % n);
        Eigen::Vector2d e1 = lV - vVert;
        Eigen::Vector2d e2 = rV - vVert;
        Eigen::Vector2d e3 = lV - rV;
        double area = (e1[0] * e2[1] - e2[0] * e1[1])/2.0;
        Eigen::Vector2d cr(-e3[1]*2*area, e3[0]*2*area);

        Eigen::Vector2d len_grad = -(e1+e2);
        double l1 = e1.squaredNorm();
        double l2 = e2.squaredNorm();
        double l3 = e3.squaredNorm();
        double sumSqN = l1 + l2 + l3;

        gradientSum += cr / (16 * pow(area, 3)) * sumSqN -
                       1 / (2 * area) * (e1 + e2);
        hessianSum += sumSqN / (16 * pow(area, 3)) * ((e3 * e3.transpose()) - e3.dot(e3) * Eigen::Matrix2d::Identity())
                      + 1/(8*pow(area, 3)) * (len_grad * cr.transpose() + cr * len_grad.transpose())
                      + 3 * sumSqN / (64 * pow(area, 5)) * (cr * cr.transpose())
                      + 1/area * Eigen::Matrix2d::Identity();

        if (area < 0 || (sqrt(l1) + sqrt(l2)) / sqrt(l3) - 1 < 1e-6)
        { // flipped triangle
            energySum += penalty_;
        }

        energySum += (l1 + l2 + l3) / (4 * area);
    }
    return {energySum, gradientSum, hessianSum};
}

bool armijo_condition(
    const double _f_curr,
    const double _f_new,
    const double _s,
    const Eigen::Vector2d& _d,
    const Eigen::Vector2d& _g,
    const double _armijo_const)
{
    return _f_new <= _f_curr + _armijo_const * _s * _d.dot(_g);
}

Eigen::Vector2d line_search(
    const Eigen::Vector2d& _x0, const Eigen::Vector2d& _d, const double _f,
    const Eigen::Vector2d& _g, const Eigen::MatrixXd& poly,
    const double _s_max, // Initial step size
    const double _shrink, const int _max_iters, const double _armijo_const)
{

    // Also try a step size of 1.0 (if valid)
    const bool try_one = _s_max > 1.0;

    Eigen::Vector2d x_new = _x0;
    double s = _s_max;
    for (int i = 0; i < _max_iters; ++i)
    {
        x_new = _x0 + s * _d;
        const double f_new = getEnergy(x_new, poly);

        if (armijo_condition(_f, f_new, s, _d, _g, _armijo_const))
            return x_new;

        if (try_one && s > 1.0 && s * _shrink < 1.0)
            s = 1.0;
        else
            s *= _shrink;
    }

    return _x0;
}

double cotan(Eigen::Vector2d& vL, Eigen::Vector2d& vC, Eigen::Vector2d& vR)
{
    Eigen::Vector2d eL = vL - vC;
    Eigen::Vector2d eR = vR - vC;
    double area = abs((-vL[1] * vR[0] + vC[1] * (-vL[0] + vR[0]) +
                       vC[0] * (vL[1] - vR[1]) + vL[0] * vR[1]));
    return eL.dot(eR) / area;
}

void computeHarmonicWeights(const Eigen::MatrixXd &poly, Eigen::VectorXd& weights, Eigen::Vector2d vV)
{
    int val = poly.rows();
    weights.resize(val);
    double wSum = 0;
    double w;

    Eigen::Vector2d vL, vI, vR;

    for (int i = 0; i < val; i++)
    {
        vL = poly.row((i + (val - 1)) % val);
        vI = poly.row(i);
        vR = poly.row((i + 1) % val);
        w = (cotan(vI, vL, vV) + cotan(vI, vR, vV)) / 2;
        weights[i] = w;
        wSum += w;
    }

    weights /= wSum;
}

void getLeastSquaresPoint(const Eigen::MatrixXd& poly, Eigen::VectorXd& x)
{
    int n = (int)poly.rows();

    Eigen::Vector3d vec = Eigen::Vector3d::Zero();
    Eigen::Matrix3d mat = Eigen::Matrix3d::Zero();
    Eigen::Vector3d p1, p2, e;
    for (int i = 0; i < n; ++i)
    {
        p1 = poly.row(i);
        p2 = poly.row((i + 1) % n);
        e = p1-p2;

        mat(0,0) += e[1]*e[1] + e[2]*e[2];
        mat(1,1) += e[0]*e[0] + e[2]*e[2];
        mat(2,2) += e[0]*e[0] + e[1]*e[1];

        mat(1,2) -= e[1]*e[2];
        mat(2,0) -= e[0]*e[2];
        mat(2,1) -= e[1]*e[2];
        mat(0,1) -= e[0]*e[1];
        mat(0,2) -= e[0]*e[2];
        mat(1,0) -= e[0]*e[1];

        vec += p1.cross(p2).cross(p1-p2);
    }
    x = mat.inverse() * vec;
}

void computeMaxAreaPolygonProjection(const Eigen::MatrixXd& poly, Eigen::Vector3d &n, Eigen::MatrixXd& pPoly, Eigen::MatrixXd& projMatrix)
{
    Eigen::MatrixXd E, B;
    E.resizeLike(poly);
    B.resizeLike(poly);

    int val = poly.rows();
    for (int i = 0; i < val; i++) {
        E.row(i) = poly.row((i + 1) % val) - poly.row(i);
        B.row(i) = 0.5 * (poly.row((i + 1) % val) + poly.row(i));
    }

    Eigen::Matrix3d A = E.transpose() * B;
    n[0] = -A(1, 2);
    n[1] = A(0, 2);
    n[2] = -A(0, 1);
    n.normalize();

    Eigen::MatrixXd tempPoly;
    tempPoly.resize(val, 3);
    for (int i = 0; i < val; i++) {
        tempPoly.row(i) = (Eigen::Vector3d) poly.row(i) - n * poly.row(i).dot(n);
    }

    Eigen::Vector3d e1 = (tempPoly.row(0) - tempPoly.row(1)).normalized();
    Eigen::Vector3d e2 = -e1.cross(n).normalized();
    Eigen::MatrixXd proj(3,2);
    proj.col(0) = e1;
    proj.col(1) = e2;

    pPoly = tempPoly * proj;
    projMatrix = proj;
}

double getStiffnessTrace(const Eigen::Vector3d& vVert,
                                          const Eigen::MatrixXd& poly,
                                          const Eigen::VectorXd& weights)
{
    int n = (int)poly.rows();
    Eigen::VectorXd s = Eigen::VectorXd::Zero(n);
    double sum = 0;
    Eigen::Vector3d lV, rV, e1, e2, e3;
    double area_coef;

    for (int i = 0; i < n; ++i)
    {
        lV = poly.row(i);
        rV = poly.row((i + 1) % n);
        e1 = lV - vVert;
        e2 = rV - vVert;
        e3 = rV - lV;
        area_coef = 2 * e1.cross(e2).norm();
        s(i) -= e2.dot(e3)/area_coef;
        s((i+1)%n) += e1.dot(e3)/area_coef;
        sum += 2*e1.dot(e2)/area_coef;
    }

    sum += 2*s.dot(weights) - s.sum() * (1 + weights.squaredNorm());
    return sum;
}