
//=============================================================================

#include <pmp/SurfaceMesh.h>
#include "PolyLaplace.h"
#include "GeodesicsInHeat.h"
#include <Eigen/Sparse>
#include <iostream>
#include <pmp/algorithms/SurfaceNormals.h>
#include <cfloat>
#include "Diamond_2D.h"
#include "DisneyLaplace.h"
#include "SECLaplace.h"
#include "diffgeo.h"
#include "LaplaceConstruction.h"

//=============================================================================

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    Disney = 5,
    SEC = 6
};

enum InsertedPoint
{
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2
};

GeodesicsInHeat::GeodesicsInHeat(pmp::SurfaceMesh &mesh, int laplace,
                                 int min_point, bool geodist, bool euklid,
                                 bool mean_edge)
    : mesh_(mesh),
      laplace_(laplace),
      min_point_(min_point),
      geodist_sphere_(geodist),
      geodist_cube_(euklid),
      mean_edge_(mean_edge)
{
    SurfaceNormals::compute_face_normals(mesh_);
    mesh_.add_face_property<Point>("f:point");
    mesh_.add_face_property<Eigen::VectorXd>("f:weights");
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

double GeodesicsInHeat::averageEdgeLength(const pmp::SurfaceMesh &mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {
        avgLen += mesh.edge_length(e);
    }

    return avgLen / mesh.n_edges();
}
//-----------------------------------------------------------------------------

double GeodesicsInHeat::maxEdgeLength(const pmp::SurfaceMesh &mesh)
{
    double maxLen = 0.;
    double currLen;
    for (auto e : mesh.edges())
    {
        currLen = mesh.edge_length(e);
        if (currLen > maxLen)
        {
            maxLen = currLen;
        }
    }

    return maxLen;
}
//-----------------------------------------------------------------------------

void GeodesicsInHeat::buildGradientOperator()
{

    // Polylaplace Gradient Operator
    if (laplace_ == AlexaLaplace)
    {
        setup_poly_gradient_operator(mesh_, gradOperator);
        gradOperator *= 0.5;
    }
    else if (laplace_ == Disney)
    {
        setup_disney_gradient_operator(mesh_, gradOperator);
    }
    else if (laplace_ == SEC)
    {
        setup_sec_gradient_operator(mesh_, gradOperator);
    }
    else
    {
        gradOperator.resize(3 * mesh_.n_faces(), mesh_.n_vertices());
        setup_sandwich_gradient_matrix(mesh_, gradOperator, min_point_);
    }
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::buildDivOperator()
{

    if (laplace_ == AlexaLaplace)
    {
        setup_poly_divergence_operator(mesh_, divOperator);
    }
    else if (laplace_ == Disney)
    {
        setup_disney_divergence_operator(mesh_, divOperator);
    }
    else if (laplace_ == SEC)
    {
        setup_sec_divergence_operator(mesh_, divOperator);
    }
    else
    {
        divOperator.resize(mesh_.n_vertices(), 3 * mesh_.n_faces());
        setup_sandwich_divergence_matrix(mesh_, divOperator, min_point_);
    }
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::compute_geodesics(bool lumped)
{

    pos.resize(mesh_.n_vertices(), 3);

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
        buildGradientOperator();
        buildDivOperator();
    }

    Eigen::SparseMatrix<double> L, A, M, M_bar;

    if (laplace_ == CotanLaplace)
    {
        //boundary handling in triangle laplace matrix
        L = divOperator * gradOperator;
    }
    else
    {
        setup_stiffness_matrices(mesh_, L, laplace_, min_point_);
    }
    setup_mass_matrices(mesh_, M, laplace_, min_point_, lumped);

    double h;
    if (mean_edge_)
    {
        h = pow(averageEdgeLength(mesh_), 2);
    }
    else
    {
        h = pow(maxEdgeLength(mesh_), 2);
    }
    A = M - h * L;

    cholL.analyzePattern(L);
    cholL.factorize(L);

    cholA.analyzePattern(A);
    cholA.factorize(A);
}
//-----------------------------------------------------------------------------

void GeodesicsInHeat::getDistance(const int vertex, Eigen::VectorXd &dist,
                                  Eigen::VectorXd &orthodist)
{

    // diffuse heat
    const int N = mesh_.n_vertices();

    auto distances = mesh_.add_vertex_property<Scalar>("v:dist");

    Eigen::SparseVector<double> b(N);
    b.coeffRef(vertex) = 1.;

    // compute gradients
    Eigen::VectorXd heat = cholA.solve(b);
    Eigen::VectorXd grad = gradOperator * heat;

    // normalize gradients
    if (laplace_ == AlexaLaplace)
    {
        normalize_poly_gradients(mesh_, grad, heat);
    }
    else if (laplace_ == SEC)
    {
        normalize_sec_gradients(mesh_, grad, heat);
    }
    else if (laplace_ == Disney)
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
            dvec2 &g = *reinterpret_cast<dvec2 *>(&grad[i]);
            double n = norm(g);
            if (n > DBL_MIN)
                g /= n;
        }
    }
    else
    {
        for (int i = 0; i < grad.rows(); i += 3)
        {
            dvec3 &g = *reinterpret_cast<dvec3 *>(&grad[i]);
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

    if (geodist_sphere_)
    {
        rms /= mesh_.n_vertices();
        rms = sqrt(rms);
        rms /= radius;
        if (laplace_ == AlexaLaplace)
        {
            std::cout << "Distance deviation sphere (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == SEC)
        {
            std::cout << "Distance deviation sphere (SEC, l="
                      << poly_laplace_lambda_ << ", lvl=" << sec_laplace_lvl
                      << ") : " << rms << std::endl;
        }
        else if (laplace_ == Disney)
        {
            std::cout << "Distance deviation sphere (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == CotanLaplace)
        {
            std::cout << "Distance deviation sphere (Cotan): " << rms
                      << std::endl;
        }
        else
        {
            if (laplace_ == SandwichLaplace)
            {
                std::cout << "Distance deviation sphere (Sandwich laplace ";
            }
            else if (laplace_ == Diamond)
            {
                std::cout << "Distance deviation sphere (Diamond laplace ";
            }
            if (min_point_ == AbsAreaMinimizer)
            {
                std::cout << "abs. area minimizer): " << rms << std::endl;
            }
            else if (min_point_ == AreaMinimizer)
            {
                std::cout << "area minimizer): " << rms << std::endl;
            }
            else
            {
                std::cout << "centroid): " << rms << std::endl;
            }
        }
    }
    if (geodist_cube_)
    {
        rms /= mesh_.n_vertices();
        rms = sqrt(rms);
        if (laplace_ == AlexaLaplace)
        {

            std::cout << "Distance deviation plane (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == SEC)
        {
            std::cout << "Distance deviation plane (SEC, l="
                      << poly_laplace_lambda_ << ", lvl=" << sec_laplace_lvl
                      << ") : " << rms << std::endl;
        }
        else if (laplace_ == Disney)
        {
            std::cout << "Distance deviation plane (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == CotanLaplace)
        {
            std::cout << "Distance deviation sphere (Cotan): " << rms
                      << std::endl;
        }
        else
        {
            if (laplace_ == SandwichLaplace)
            {
                std::cout << "Distance deviation plane (Sandwich laplace ";
            }
            else if (laplace_ == Diamond)
            {
                std::cout << "Distance deviation plane (Diamond laplace ";
            }
            if (min_point_ == AbsAreaMinimizer)
            {
                std::cout << "abs. area minimizer): " << rms << std::endl;
            }
            else if (min_point_ == AreaMinimizer)
            {
                std::cout << "area minimizer): " << rms << std::endl;
            }
            else
            {
                std::cout << "centroid): " << rms << std::endl;
            }
        }
    }

    distance_to_texture_coordinates();

    mesh_.remove_vertex_property<Scalar>(distances);
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::write_Distance_error(const int vertex,
                                           Eigen::VectorXd &dist,
                                           Eigen::VectorXd &orthodist,
                                           bool lumped, std::ofstream &file)
{

    // diffuse heat
    const int N = mesh_.n_vertices();

    auto distances = mesh_.add_vertex_property<Scalar>("v:dist");

    Eigen::SparseVector<double> b(N);
    b.coeffRef(vertex) = 1.;
    Eigen::SparseMatrix<double> M;
    setup_mass_matrices(mesh_, M, laplace_, min_point_, lumped);
    // compute gradients
    Eigen::VectorXd heat = cholA.solve(b);
    Eigen::VectorXd grad = gradOperator * heat;
    // normalize gradients
    if (laplace_ == AlexaLaplace)
    {
        normalize_poly_gradients(mesh_, grad, heat);
    }
    else if (laplace_ == SEC)
    {
        normalize_sec_gradients(mesh_, grad, heat);
    }
    else if (laplace_ == Disney)
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
            dvec2 &g = *reinterpret_cast<dvec2 *>(&grad[i]);
            double n = norm(g);
            if (n > DBL_MIN)
                g /= n;
        }
    }
    else
    {
        for (int i = 0; i < grad.rows(); i += 3)
        {
            dvec3 &g = *reinterpret_cast<dvec3 *>(&grad[i]);
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
        }
        if (geodist_cube_)
        {
            orthodist(k) = norm(mesh_.position(v0) - mesh_.position(v));
        }
        rms += (dist(k) - orthodist(k)) * (dist(k) - orthodist(k));

        k++;
    }

    if (geodist_sphere_)
    {
        rms /= mesh_.n_vertices();
        rms = sqrt(rms);
        rms /= radius;
        if (laplace_ == AlexaLaplace)
        {
            std::cout << "Distance deviation sphere (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == SEC)
        {
            std::cout << "Distance deviation sphere (SEC, l="
                      << poly_laplace_lambda_ << ", lvl=" << sec_laplace_lvl
                      << ") : " << rms << std::endl;
        }
        else if (laplace_ == Disney)
        {
            std::cout << "Distance deviation sphere (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == CotanLaplace)
        {
            std::cout << "Distance deviation sphere (Cotan): " << rms
                      << std::endl;
        }
        else
        {
            if (laplace_ == SandwichLaplace)
            {
                std::cout << "Distance deviation sphere (Sandwich laplace ";
            }
            else if (laplace_ == Diamond)
            {
                std::cout << "Distance deviation sphere (Diamond laplace ";
            }
            if (min_point_ == AbsAreaMinimizer)
            {
                std::cout << "abs. area minimizer): " << rms << std::endl;
            }
            else if (min_point_ == AreaMinimizer)
            {
                std::cout << "area minimizer): " << rms << std::endl;
            }
            else
            {
                std::cout << "centroid): " << rms << std::endl;
            }
        }
    }
    if (geodist_cube_)
    {
        rms /= mesh_.n_vertices();
        rms = sqrt(rms);
        if (laplace_ == AlexaLaplace)
        {

            std::cout << "Distance deviation plane (Alexa, l="
                      << poly_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == SEC)
        {
            std::cout << "Distance deviation plane (SEC, l="
                      << poly_laplace_lambda_ << ", lvl=" << sec_laplace_lvl
                      << ") : " << rms << std::endl;
        }
        else if (laplace_ == Disney)
        {
            std::cout << "Distance deviation plane (Disney, l="
                      << disney_laplace_lambda_ << "): " << rms << std::endl;
        }
        else if (laplace_ == CotanLaplace)
        {
            std::cout << "Distance deviation sphere (Cotan): " << rms
                      << std::endl;
        }
        else
        {
            if (laplace_ == SandwichLaplace)
            {
                std::cout << "Distance deviation plane (Sandwich laplace ";
            }
            else if (laplace_ == Diamond)
            {
                std::cout << "Distance deviation plane (Diamond laplace ";
            }
            if (min_point_ == AbsAreaMinimizer)
            {
                std::cout << "abs. area minimizer): " << rms << std::endl;
            }
            else if (min_point_ == AreaMinimizer)
            {
                std::cout << "area minimizer): " << rms << std::endl;
            }
            else
            {
                std::cout << "centroid): " << rms << std::endl;
            }
        }
    }
    file << rms << ",";
    mesh_.remove_vertex_property<Scalar>(distances);
}

//-----------------------------------------------------------------------------

void GeodesicsInHeat::distance_to_texture_coordinates() const
{

    auto distances = mesh_.get_vertex_property<Scalar>("v:dist");
    assert(distances);
    // scale with boundingbox size for comparable geodesic rings
    Scalar bb_size(0);
    bb_size = mesh_.bounds().size();
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    for (auto v : mesh_.vertices())
    {
        if (distances[v] <= FLT_MAX)
        {
            tex[v] = TexCoord(distances[v] / bb_size, 0.0);
        }
        else
        {
            tex[v] = TexCoord(1.0, 0.0);
        }
    }

    // remove per-halfedge texture coordinates
    auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
    if (htex)
        mesh_.remove_halfedge_property(htex);
}

//-----------------------------------------------------------------------------

double GeodesicsInHeat::great_circle_distance(Vertex v, Vertex vv, double r)
{
    double dis;
    if (v == vv)
    {
        return 0.0;
    }
    Normal n = pmp::SurfaceNormals::compute_vertex_normal(mesh_, v);
    Normal nn = pmp::SurfaceNormals::compute_vertex_normal(mesh_, vv);
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
