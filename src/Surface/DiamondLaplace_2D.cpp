//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "DiamondLaplace_2D.h"
#include "diffgeo.h"
#include <pmp/algorithms/DifferentialGeometry.h>

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

using namespace std;

enum InsertedPoint
{
    Centroid = 0,
    AreaMinimizer = 2,
};

//----------------------------------------------------------------------------------

void compute_primal_points(SurfaceMesh &mesh, int minpoint)
{

    FaceProperty<pmp::Point> area_points = mesh.get_face_property<pmp::Point>("f:point");
    Eigen::VectorXd w;
    Eigen::MatrixXd poly;
    Eigen::Vector3d point;
    for (Face f : mesh.faces())
    {
        const int n = mesh.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh.vertices(f))
        {
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        if (minpoint == Centroid)
        {
            int val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min = poly.transpose() * w;
        Point p = Point(min[0], min[1], min[2]);
        area_points[f] = p;
    }
}

//=============================================================================

void setup_diamond_mass_matrix(SurfaceMesh &mesh,
                               Eigen::SparseMatrix<double> &M)
{

    auto primal_points = mesh.get_face_property<Point>("f:point");
    std::vector<Triplet> triplets_area;
    const int nv = mesh.n_vertices();
    const int nf = mesh.n_faces();
    Point F1, F2, T1, T2;
    M.resize(nv + nf, nv + nf);
    Vertex v, vv;
    Face f, ff;
    int f1_idx, f2_idx, t1_idx, t2_idx;
    for (Edge e : mesh.edges())
    {

        if (!mesh.is_boundary(e))
        {
            Halfedge he = mesh.halfedge(e, 0);

            v = mesh.from_vertex(he);
            vv = mesh.to_vertex(he);

            f = mesh.face(he);
            ff = mesh.face(mesh.opposite_halfedge(he));

            F1 = mesh.position(v);
            F2 = mesh.position(vv);

            f1_idx = v.idx();
            f2_idx = vv.idx();

            T2 = primal_points[f];
            T1 = primal_points[ff];

            t2_idx = nv + f.idx();
            t1_idx = nv + ff.idx();
        }
        else
        {
            Halfedge he = mesh.halfedge(e, 0);
            v = mesh.from_vertex(he);
            vv = mesh.to_vertex(he);
            F1 = mesh.position(v);
            F2 = mesh.position(vv);

            f1_idx = v.idx();
            f2_idx = vv.idx();

            if (!mesh.is_boundary(he))
            {
                f = mesh.face(he);
                T2 = primal_points[f];
                T1 = 0.5 * (F1 + F2);
                t2_idx = nv + f.idx();
            }
            else
            {
                ff = mesh.face(mesh.opposite_halfedge(he));
                T2 = primal_points[ff];
                T1 = 0.5 * (F1 + F2);
                t2_idx = nv + ff.idx();
            }
        }

        Point e2 = normalize(F2 - F1);
        Point n, e1;
        vec2 f1, f2;

        n = normalize(cross(T2 - F2, F2 - F1));
        e1 = cross(n, e2);

        f1 = vec2(0.0, 0.0);
        f2 = vec2(0.0, norm(F1 - F2));

        double d_area;

        if (!mesh.is_boundary(e))
        {
            d_area =
                (triangle_area(T1, F2, F1) + triangle_area(F1, F2, T2)) / 4.0;

            triplets_area.emplace_back(f1_idx, f1_idx, d_area);
            triplets_area.emplace_back(f2_idx, f2_idx, d_area);
            triplets_area.emplace_back(t2_idx, t2_idx, d_area);
            triplets_area.emplace_back(t1_idx, t1_idx, d_area);
        }
        else
        {

            d_area =
                (triangle_area(T1, F2, F1) + triangle_area(F1, F2, T2)) / 3.0;
            triplets_area.emplace_back(f1_idx, f1_idx, d_area);
            triplets_area.emplace_back(f2_idx, f2_idx, d_area);
            triplets_area.emplace_back(t2_idx, t2_idx, d_area);
        }
    }

    M.setFromTriplets(triplets_area.begin(), triplets_area.end());
}

//====================================================================================================

void setup_diamond_gradient_divergence_intrinsic(SurfaceMesh &mesh,
                                                 Eigen::SparseMatrix<double> &G,
                                                 Eigen::SparseMatrix<double> &D)
{
    auto primal_points = mesh.get_face_property<Point>("f:point");
    Eigen::SparseMatrix<double> M;
    std::vector<Triplet> triplets, triplets_area;
    const int nv = mesh.n_vertices();
    const int ne = mesh.n_edges();
    const int nf = mesh.n_faces();
    Point F1, F2, T1, T2;

    G.resize(2 * ne, nv + nf);
    D.resize(nv + nf, 2 * ne);
    M.resize(2 * ne, 2 * ne);

    Vertex v, vv;
    Face f, ff;

    for (Edge e : mesh.edges())
    {

        if (!mesh.is_boundary(e))
        {
            Halfedge he = mesh.halfedge(e, 0);

            v = mesh.from_vertex(he);
            vv = mesh.to_vertex(he);

            f = mesh.face(he);
            ff = mesh.face(mesh.opposite_halfedge(he));

            F1 = mesh.position(v);
            F2 = mesh.position(vv);

            T2 = primal_points[f];
            T1 = primal_points[ff];
        }
        else
        {
            Halfedge he = mesh.halfedge(e, 0);
            v = mesh.from_vertex(he);
            vv = mesh.to_vertex(he);
            F1 = mesh.position(v);
            F2 = mesh.position(vv);

            if (!mesh.is_boundary(he))
            {
                f = mesh.face(he);
                T2 = primal_points[f];
                T1 = 0.5 * (F1 + F2);
            }
            else
            {
                ff = mesh.face(mesh.opposite_halfedge(he));
                T2 = primal_points[ff];
                T1 = 0.5 * (F1 + F2);
            }
        }

        //=============================2D point projection==================================

        Point e2 = normalize(F2 - F1);
        Point n, e1;
        vec2 f1, f2, t1, t2;

        n = normalize(cross(T2 - F2, F2 - F1));
        e1 = cross(n, e2);

        f1 = vec2(0.0, 0.0);
        f2 = vec2(0.0, norm(F1 - F2));

        t2 = dot((T2 - F1), e1) * vec2(1.0, 0) +
             dot((T2 - F1), e2) * vec2(0.0, 1.0);

        if (!mesh.is_boundary(e))
        {
            t1 = dot((T1 - F1), e1) * vec2(1.0, 0.0) +
                 dot((T1 - F1), e2) * vec2(0.0, 1.0);
        }
        else
        {
            t1 = 0.5 * (f1 + f2);
        }

        //=============================Gradient Construction ==================================
        Eigen::MatrixXd Diff, X, Diff_X, A, A_inv, G_, G2, E;

        X.resize(4, 2);
        G2.resize(2, 4);

        for (int i = 0; i < 2; i++)
        {
            X(0, i) = f1[i];
            X(1, i) = t2[i];
            X(2, i) = f2[i];
            X(3, i) = t1[i];
        }

        double denom = (t2[0] - t1[0]) * (f1[1] - f2[1]) -
                       (f1[0] - f2[0]) * (t2[1] - t1[1]);

        G2.row(0) << (t1[1] - t2[1]), (f1[1] - f2[1]), (t2[1] - t1[1]),
            (f2[1] - f1[1]);
        G2.row(1) << (t2[0] - t1[0]), (f2[0] - f1[0]), (t1[0] - t2[0]),
            (f1[0] - f2[0]);

        G2 /= denom;

        double area1 = 0.5 * (-f1[1] * f2[0] + t2[1] * (-f1[0] + f2[0]) +
                              t2[0] * (f1[1] - f2[1]) + f1[0] * f2[1]);
        double area2 = 0.5 * (-f1[1] * t1[0] + f2[1] * (-f1[0] + t1[0]) +
                              f2[0] * (f1[1] - t1[1]) + f1[0] * t1[1]);

        if (abs(denom) < 1e-8)
        {
            G2.setZero();
        }

        double diamond_area = area1 + area2;
        for (int i = 0; i < 2; i++)
        {
            int row = 2 * e.idx() + i;
            triplets_area.emplace_back(row, row, diamond_area);

            if (!mesh.is_boundary(e))
            {
                triplets.emplace_back(row, nv + ff.idx(), G2(i, 3));
                triplets.emplace_back(row, nv + f.idx(), G2(i, 1));
            }
            else
            {
                Halfedge he = mesh.halfedge(e, 0);
                if (!mesh.is_boundary(he))
                {
                    f = mesh.face(he);
                    triplets.emplace_back(row, nv + f.idx(), G2(i, 1));
                    triplets.emplace_back(row, vv.idx(), 0.5 * G2(i, 3));
                    triplets.emplace_back(row, v.idx(), 0.5 * G2(i, 3));
                }
                else
                {
                    ff = mesh.face(mesh.opposite_halfedge(he));
                    triplets.emplace_back(row, nv + ff.idx(), G2(i, 1));
                    triplets.emplace_back(row, vv.idx(), 0.5 * G2(i, 3));
                    triplets.emplace_back(row, v.idx(), 0.5 * G2(i, 3));
                }
            }
            triplets.emplace_back(row, vv.idx(), G2(i, 2));
            triplets.emplace_back(row, v.idx(), G2(i, 0));
        }
    }
    G.setFromTriplets(triplets.begin(), triplets.end());
    M.setFromTriplets(triplets_area.begin(), triplets_area.end());
    D = -G.transpose() * M;
}

//====================================================================================================
