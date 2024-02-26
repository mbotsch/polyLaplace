//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "[AW11]Laplace.h"
#include "GeodesicsInHeat.h"
#include <Eigen/Sparse>
#include <iostream>
#include <pmp/algorithms/normals.h>
#include <cfloat>
#include "DiamondLaplace_2D.h"
#include "[dGBD20]Laplace.h"
#include "diffgeo.h"
#include "LaplaceConstruction.h"
#include "SpectralProcessing.h"

//=============================================================================

GeodesicsInHeat::GeodesicsInHeat(pmp::SurfaceMesh& mesh, int laplace,
                                 int min_point, bool geodist, bool euklid,
                                 DiffusionStep diffusion)
    : mesh_(mesh),
      laplace_(laplace),
      min_point_(min_point),
      geodist_sphere_(geodist),
      geodist_cube_(euklid),
      diffusionStep_(diffusion)
{
    if (!mesh.has_face_property("f:point") ||
        !mesh.has_face_property("f:weights"))
    {
        mesh_.add_face_property<pmp::Point>("f:point");
        mesh_.add_face_property<Eigen::VectorXd>("f:weights");
    }
    setup_face_point_properties(mesh_, min_point);
}

//-----------------------------------------------------------------------------

GeodesicsInHeat::~GeodesicsInHeat()
{
    auto area_points = mesh_.face_property<Point>("f:point");
    auto area_weights = mesh_.face_property<Eigen::VectorXd>("f:weights");

    mesh_.remove_face_property(area_points);
    mesh_.remove_face_property(area_weights);
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::edgeLength(const pmp::SurfaceMesh& mesh, pmp::Edge e)
{
    auto v0 = mesh.vertex(e, 0);
    auto v1 = mesh.vertex(e, 1);
    return pmp::distance(mesh.position(v0), mesh.position(v1));
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::averageEdgeLength(const pmp::SurfaceMesh& mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {
        avgLen += edgeLength(mesh, e);
    }

    return avgLen / (double)mesh.n_edges();
}
//-----------------------------------------------------------------------------

double GeodesicsInHeat::maxEdgeLength(const pmp::SurfaceMesh& mesh)
{
    double maxLen = 0.;
    double currLen;
    for (auto e : mesh.edges())
    {
        currLen = edgeLength(mesh, e);
        if (currLen > maxLen)
        {
            maxLen = currLen;
        }
    }

    return maxLen;
}
//-----------------------------------------------------------------------------

double GeodesicsInHeat::maxDiagonalLength(const pmp::SurfaceMesh& mesh)
{
    double maxDiag = 0.;
    double currLen;
    for (auto f : mesh.faces())
    {
        for (auto v : mesh.vertices(f))
        {
            for (auto vv : mesh.vertices(f))
            {
                currLen = distance(mesh.position(v), mesh.position(vv));
                if (currLen > maxDiag)
                {
                    maxDiag = currLen;
                }
            }
        }
    }
    return maxDiag;
}
//-----------------------------------------------------------------------------

void GeodesicsInHeat::compute_geodesics(double& condition_number, bool lumped)
{
    pos.resize((int)mesh_.n_vertices(), 3);

    for (int i = 0; i < (int)mesh_.n_vertices(); ++i)
        for (int j = 0; j < 3; ++j)
            pos(i, j) = mesh_.positions()[i][j];

    Eigen::SparseMatrix<double> P;
    setup_prolongation_matrix(mesh_, P);

    if (laplace_ == Diamond)
    {
        Eigen::SparseMatrix<double> G, D;
        compute_primal_points(mesh_, min_point_);
        setup_diamond_gradient_divergence_intrinsic(mesh_, G, D);
        gradOperator = G * P;
        divOperator = P.transpose() * D;
    }
    else
    {
        setup_gradient(mesh_, gradOperator, laplace_, min_point_);
        setup_divergence(mesh_, divOperator, laplace_, min_point_);
    }
    Eigen::SparseMatrix<double> S, A, M;

    setup_stiffness_matrices(mesh_, S, laplace_, min_point_);
    setup_mass_matrices(mesh_, M, laplace_, min_point_, lumped);

    double h;
    if (diffusionStep_ == MeanEdge)
    {
        h = pow(averageEdgeLength(mesh_), 2);
        std::cout << "Mean edge length diffusion \n";
    }
    else if (diffusionStep_ == MaxEdge)
    {
        h = pow(maxEdgeLength(mesh_), 2);
        std::cout << "Max edge length diffusion \n";
    }
    else
    {
        h = pow(maxDiagonalLength(mesh_), 2);
        std::cout << "Max diagonal length diffusion \n";
    }

    A = M - h * S;

    condition_number = get_condition_number(A, false);

    cholL.analyzePattern(S);
    cholL.factorize(S);

    cholA.analyzePattern(A);
    cholA.factorize(A);
}
//-----------------------------------------------------------------------------

double GeodesicsInHeat::getDistance(const int vertex, Eigen::VectorXd& dist,
                                    Eigen::VectorXd& orthodist, bool verbose)
{
    // diffuse heat
    const int N = (int)mesh_.n_vertices();

    auto distances = mesh_.add_vertex_property<Scalar>("v:dist");

    Eigen::SparseVector<double> b(N);
    b.coeffRef(vertex) = 1.;

    // compute gradients
    Eigen::VectorXd heat = cholA.solve(b);
    Eigen::VectorXd grad = gradOperator * heat;

    // normalize gradients
    if (laplace_ == AlexaWardetzkyLaplace)
    {
        normalize_poly_gradients(mesh_, grad, heat);
    }
    else if (laplace_ == deGoesLaplace)
    {
        for (int i = 0; i < grad.rows() / 3; i++)
        {
            dvec3 g = dvec3(grad(3 * i), grad(3 * i + 1), grad(3 * i + 2));
            double n = norm(g);
            if (n > DBL_MIN)
                g /= n;

            grad(3 * i) = g[0];
            grad(3 * i + 1) = g[1];
            grad(3 * i + 2) = g[2];
        }
    }
    else if (laplace_ == Diamond)
    {
        for (int i = 0; i < grad.rows(); i += 2)
        {
            dvec2& g = *reinterpret_cast<dvec2*>(&grad[i]);
            double n = norm(g);
            if (n > DBL_MIN)
                g /= n;
        }
    }
    else
    {
        for (int i = 0; i < grad.rows(); i += 3)
        {
            dvec3& g = *reinterpret_cast<dvec3*>(&grad[i]);
            double n = norm(g);
            if (n > DBL_MIN)
                g /= n;
        }
    }
    dist = cholL.solve(divOperator * (-grad));

    orthodist.resize(dist.size());

    double mi = dist.minCoeff();
    for (int i = 0; i < dist.rows(); ++i)
        dist[i] -= mi;

    int k = 0;
    Vertex v0 = Vertex(vertex);
    double rms = 0.0;
    double radius = norm(mesh_.position(v0));
    for (auto v : mesh_.vertices())
    {
        distances[v] = dist[k];

        if (geodist_sphere_)
        {
            orthodist(k) = great_circle_distance(v0, v, radius);
            rms += (dist(k) - orthodist(k)) * (dist(k) - orthodist(k));
        }

        if (geodist_cube_)
        {
            orthodist(k) = norm(mesh_.position(v0) - mesh_.position(v));
            rms += (dist(k) - orthodist(k)) * (dist(k) - orthodist(k));
        }

        k++;
    }
    std::string meshtype;
    if (geodist_sphere_)
    {
        rms /= (double)mesh_.n_vertices();
        rms = sqrt(rms);
        rms /= radius;
        meshtype = " sphere ";
    }
    else if (geodist_cube_)
    {
        rms /= (double)mesh_.n_vertices();
        rms = sqrt(rms);
        meshtype = " plane ";
    }
    if ((geodist_cube_ || geodist_sphere_) && verbose)
    {
        if (laplace_ == AlexaWardetzkyLaplace)
        {
            std::cout << "Distance deviation" + meshtype +
                             "(AlexaWardetzky Laplace, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == deGoesLaplace)
        {
            std::cout << "Distance deviation" + meshtype + "(deGoes Laplace, l="
                      << deGoes_laplace_lambda_ << "): " << rms << std::endl;
        }
        else
        {
            if (laplace_ == PolySimpleLaplace)
            {
                std::cout << "Distance deviation" + meshtype + "(Polysimple ";
            }
            else if (laplace_ == Diamond)
            {
                std::cout << "Distance deviation" + meshtype + "(Diamond ";
            }

            if (min_point_ == AreaMinimizer)
            {
                std::cout << "squared area minimizer): " << rms << std::endl;
            }
            else if (min_point_ == Centroid_)
            {
                std::cout << "centroid): " << rms << std::endl;
            }
            else
            {
                std::cout << "trace minimizer): " << rms << std::endl;
            }
        }
    }
    distance_to_texture_coordinates();
    mesh_.remove_vertex_property<Scalar>(distances);
    return rms;
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::distance_to_texture_coordinates() const
{
    //     remove per-halfedge texture coordinates
    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh_.remove_halfedge_property(htex);

    auto distances = mesh_.get_vertex_property<Scalar>("v:dist");
    assert(distances);

    // find maximum distance
    Scalar maxdist(0);
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            maxdist = std::max(maxdist, distances[v]);
        }
    }

    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            tex[v] = TexCoord(distances[v] / maxdist, 0.0);
        }
        else
        {
            tex[v] = TexCoord(1.0, 0.0);
        }
    }
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::great_circle_distance(Vertex v, Vertex vv, double r)
{
    double dis;
    if (v == vv)
    {
        return 0.0;
    }
    Normal n = pmp::vertex_normal(mesh_, v);
    Normal nn = pmp::vertex_normal(mesh_, vv);
    double delta_sigma = acos(dot(n, nn));
    if (std::isnan(delta_sigma))
    {
        dis = haversine_distance(v, vv, r);
        if (std::isnan(delta_sigma))
        {
            dis = vincenty_distance(v, vv, r);
        }
        return dis;
    }

    return r * delta_sigma;
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::haversine_distance(Vertex v, Vertex vv, double r)
{
    Point p = mesh_.position(v);
    Point pp = mesh_.position(vv);

    double lamda1 = atan2(p[1], p[0]) + M_PI;
    double phi1 = M_PI / 2.0 - acos(p[2] / r);

    double lamda2 = atan2(pp[1], pp[0]) + M_PI;
    double phi2 = M_PI / 2.0 - acos(pp[2] / r);

    double d_lamda = fabs(lamda1 - lamda2);
    double d_phi = fabs(phi1 - phi2);

    double a = pow(sin(d_phi / 2), 2) +
               cos(phi1) * cos(phi2) * pow(sin(d_lamda / 2), 2);

    double d_sigma = 2 * asin(sqrt(a));

    return r * d_sigma;
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::vincenty_distance(Vertex v, Vertex vv, double r)
{
    //  special case of the Vincenty formula for an ellipsoid with equal major and minor axes
    Point p = mesh_.position(v);
    Point pp = mesh_.position(vv);

    double lamda1 = atan2(p[1], p[0]) + M_PI;
    double phi1 = M_PI / 2.0 - acos(p[2] / r);

    double lamda2 = atan2(pp[1], pp[0]) + M_PI;
    double phi2 = M_PI / 2.0 - acos(pp[2] / r);

    double d_lamda = fabs(lamda1 - lamda2);

    // Numerator
    double a = pow(cos(phi2) * sin(d_lamda), 2);

    double b = cos(phi1) * sin(phi2);
    double c = sin(phi1) * cos(phi2) * cos(d_lamda);
    double d = pow(b - c, 2);

    double e = sqrt(a + d);

    // Denominator
    double f = sin(phi1) * sin(phi2);
    double g = cos(phi1) * cos(phi2) * cos(d_lamda);

    double h = f + g;

    double d_sigma = atan2(e, h);

    return r * d_sigma;
}

//=============================================================================
