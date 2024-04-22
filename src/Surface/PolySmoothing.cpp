//=============================================================================
// Copyright 2024 Astrid Bunge, Sven Wagner, Dennis Bukenberger, Mario Botsch, Marc Alexa
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "../common_util.h"

#include "SpectralProcessing.h"
#include "PolySmoothing.h"
#include <pmp/algorithms/normals.h>
#include <pmp/algorithms/utilities.h>
#include <random>
#include "diffgeo.h"

#include <iostream>
#include <TinyAD/Support/pmp.hh>

#include <TinyAD/Scalar.hh>
#include <TinyAD/ScalarFunction.hh>
#include <TinyAD/Utils/NewtonDirection.hh>
#include <TinyAD/Utils/NewtonDecrement.hh>
#include <TinyAD/Utils/LineSearch.hh>
#include <TinyAD/Utils/Timer.hh>

//=============================================================================

using namespace pmp;

//=============================================================================

template <typename T>
T cotan(Eigen::Vector3<T> vL, Eigen::Vector3<T> vC, Eigen::Vector3<T> vR)
{
    Eigen::Vector3<T> eL = vL - vC;
    Eigen::Vector3<T> eR = vR - vC;
    return eL.dot(eR) / eL.cross(eR).norm();
}

PolySmoothing::PolySmoothing(SurfaceMesh& mesh, SmoothingConfigs oConf)
    : mesh_(mesh), oConf_(oConf), penalty_(100000)
{
    // Setup Bounding Boxes
    BoundingBox bb = bounds(mesh_);
    is2D_ = bb.min()[2] == bb.max()[2];

    setupQuadrics();
    computeVirtualVertices();
    if (oConf_.fixBoundary)
    {
        auto originalPosition =
            mesh_.vertex_property<Eigen::Vector3d>("v:originalPos");
        for (auto v : mesh_.vertices())
        {
            originalPosition[v] = mesh_.position(v);
        }
    }
}

void PolySmoothing::optimize(int quadricsTau)
{
    auto faceVirtuals = mesh_.get_face_property<Eigen::Vector3d>("f:Virtuals");
    auto vertexQuadrics =
        mesh_.get_vertex_property<Eigen::Matrix4d>("v:Quadrics");
    auto originalPosition =
        mesh_.vertex_property<Eigen::Vector3d>("v:originalPos");

    auto func = TinyAD::scalar_function<3>(mesh_.vertices());

    func.add_elements<2>(
        mesh_.halfedges(), [&](auto& element) -> TINYAD_SCALAR_TYPE(element) {
            using T = TINYAD_SCALAR_TYPE(element);
            Halfedge hIdx = element.handle;
            if (!mesh_.is_boundary(hIdx)){
                Vertex v1 = mesh_.from_vertex(hIdx);
                Vertex v2 = mesh_.to_vertex(hIdx);
                Face f = mesh_.face(hIdx);

                Eigen::Vector3<T> a = faceVirtuals[f]; // virtual face vertex
                Eigen::Vector3<T> b = element.variables(v1);
                Eigen::Vector3<T> c = element.variables(v2);

                Eigen::Vector3<T> n = (b - a).cross(c - a);
                Eigen::Vector3d fn = face_normal(mesh_, f);

                T area2 = n.squaredNorm();

                if (area2 < 1e-12)
                    return penalty_;

                if (fn.dot(n) < 0)
                    return penalty_;

                return cotan(c, a, b) + cotan(a, b, c) + cotan(b, c, a);
            }
            return 0.0;
        });

    if (quadricsTau >= 0)
    {
        func.add_elements<1>(
            mesh_.vertices(),
            [&](auto& element) -> TINYAD_SCALAR_TYPE(element) {
                using T = TINYAD_SCALAR_TYPE(element);
                Vertex v = element.handle;

                Eigen::Vector3<T> p = element.variables(v);
                Eigen::Vector4<T> ph(p.data());
                ph[3] = 1;

                if (quadricsTau < 0)
                {
                    return pow(1 - p.norm(), 2) * penalty_;
                }

                T err = 0;
                err += ph.transpose() * vertexQuadrics[v] * ph;

                if (is2D_ && mesh_.is_boundary(v) && oConf_.fixBoundary)
                {
                    Eigen::Vector3<T> q = originalPosition[v];
                    err += (p - q).dot(p - q) * penalty_;
                }
                return err * pow(10, quadricsTau);
            });
    }

    Eigen::VectorXd x =
        func.x_from_data([&](Vertex v) { return mesh_.position(v); });

    Eigen::Vector3d values(0, 0, 0);

    double convergence_eps = 1e-4;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_solver;

    double oldNd = 0;
    for (int i = 0; i < oConf_.numIters; ++i)
    {
        auto [e, g, H_proj] = func.eval_with_hessian_proj(x);

        Eigen::VectorXd d = llt_solver.compute(H_proj).solve(-g);

        double nd = TinyAD::newton_decrement(d, g);

        if (oConf_.withCnum)
            condition_number(mesh_, PolySimpleLaplace, TraceMinimizer, values,
                             oConf_.generalizedCnum);

        TINYAD_DEBUG_OUT("It: " << i << ", E: " << e << ", nd: " << nd
                                << ", evMx: " << values[0] << ", evMn: "
                                << values[1] << ", cNum: " << values[2])

        if (nd < convergence_eps || oldNd == nd)
            break;

        x = TinyAD::line_search(x, d, e, g, func);

        func.x_to_data(x, [&](Vertex v, const Eigen::Vector3d& p) {
            mesh_.position(v) = p;
        });

        computeVirtualVertices(false);

        if (i && oConf_.updateQuadrics)
            setupQuadrics();

        // check if iteration is stuck
        if (!(i % 10))
            oldNd = nd;
    }

    TINYAD_DEBUG_OUT("final energy: " << func.eval(x))

    mesh_.remove_vertex_property(originalPosition);
    mesh_.remove_vertex_property(vertexQuadrics);
    mesh_.remove_face_property(faceVirtuals);
}

void PolySmoothing::computeVirtualVertices(bool use_fallback)
{
    auto faceVirtuals = mesh_.face_property<Eigen::Vector3d>("f:Virtuals");
    Eigen::MatrixXd poly;
    Eigen::Vector3d vV;
    Eigen::VectorXd weights;

    for (auto f : mesh_.faces())
    {
        get_polygon_from_face(mesh_, f, poly);
        find_trace_minimizer_weights(poly, weights, use_fallback);
        vV = poly.transpose() * weights;
        faceVirtuals[f] = vV;
    }
}

void PolySmoothing::setupQuadrics()
{
    auto vertexQuadrics = mesh_.vertex_property<Eigen::Matrix4d>("v:Quadrics");

    for (Vertex v : mesh_.vertices())
    {
        vertexQuadrics[v].setZero();
    }

    Eigen::Vector3d n, c;

    for (Face f : mesh_.faces())
    {
        n = face_normal(mesh_, f);
        Eigen::Vector4d qv(n.data());
        for (Vertex v : mesh_.vertices(f))
        {
            qv[3] = -n.dot((Eigen::Vector3d)mesh_.position(v));
            vertexQuadrics[v] += qv * qv.transpose();
        }

        for (Halfedge he : mesh_.halfedges(f))
        {
            Edge e = mesh_.edge(he);
            if (mesh_.is_boundary(e))
            {
                Vertex v0 = mesh_.vertex(e, 0);
                Vertex v1 = mesh_.vertex(e, 1);
                Eigen::Vector3d p0 = mesh_.position(v0);
                Eigen::Vector3d p1 = mesh_.position(v1);
                Eigen::Vector3d en = n.cross(p1 - p0).normalized();

                Eigen::Vector4d qv_(en.data());
                qv_[3] = -en.dot(p0);
                vertexQuadrics[v0] += qv_ * qv_.transpose();
                vertexQuadrics[v1] += qv_ * qv_.transpose();
            }
        }
    }
}