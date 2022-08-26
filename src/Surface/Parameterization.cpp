//=============================================================================

#include "Parameterization.h"
#include "diffgeo.h"
#include "PolyLaplace.h"
#include "LaplaceConstruction.h"
#include <igl/slice.h>

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

bool Parameterization::harmonic(unsigned int laplace, unsigned int min_point)
{
    // map boundary to circle
    if (!setup_boundary_constraints())
    {
        std::cerr << "Could not perform setup of boundary constraints.\n";
        return false;
    }

    // get property for 2D texture coordinates
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");
    // Compute the chosen implicit points and its convex combination weights for each face

    int nb = 0;
    for (auto v : mesh_.vertices())
    {
        if (mesh_.is_boundary(v))
        {
            nb++;
        }
    }
    const unsigned int nv = mesh_.n_vertices();

    unsigned k = 0;
    unsigned ins = 0;
    unsigned out = 0;

    Eigen::MatrixXd B(nv, 2);
    std::vector<Vertex> vertices;
    Vertex v;
    vertices.reserve(mesh_.n_vertices());
    Eigen::VectorXi in(nv - nb), b(nb), x(2);
    x << 0, 1;
    for (auto v : mesh_.vertices())
    {
        vertices.push_back(v);
        if (!mesh_.is_boundary(v))
        {
            in(ins) = k;
            B(k, 0) = 0.0;
            B(k, 1) = 0.0;
            ins++;
        }
        else
        {
            b(out) = k;
            B(k, 0) = tex[v][0];
            B(k, 1) = tex[v][1];
            out++;
        }
        k++;
    }

    SparseMatrix L;
    Eigen::SparseMatrix<double> L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out, X;

    setup_stiffness_matrices(mesh_,L, laplace,min_point);

    igl::slice(L, in, in, L_in_in);
    igl::slice(L, in, b, L_in_b);
    igl::slice(B, in, x, b_in);
    igl::slice(B, b, x, b_out);

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L_in_in);
    X = solver.solve(b_in - L_in_b * b_out);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "harmonic(): Could not solve linear system\n";
    }
    else
    {
        // copy solution
        k = 0;
        for (unsigned int i = 0; i < nv; ++i)
        {
            v = vertices[i];
            if (!mesh_.is_boundary(v))
            {
                tex[v][0] = X(k, 0);
                tex[v][1] = X(k, 1);
                k++;
            }
        }
    }
    return true;
}

//-----------------------------------------------------------------------------

bool Parameterization::setup_boundary_constraints()
{
    // get properties
    auto points = mesh_.vertex_property<Point>("v:point");
    auto tex = mesh_.vertex_property<TexCoord>("v:tex");

    SurfaceMesh::VertexIterator vit, vend = mesh_.vertices_end();
    Vertex vh;
    Halfedge hh;
    std::vector<Vertex> loop;

    // Initialize all texture coordinates to the origin.
    for (auto v : mesh_.vertices())
        tex[v] = TexCoord(0.5, 0.5);

    // find 1st boundary vertex
    for (vit = mesh_.vertices_begin(); vit != vend; ++vit)
        if (mesh_.is_boundary(*vit))
            break;

    // no boundary found ?
    if (vit == vend)
    {
        std::cerr << "Mesh has no boundary." << std::endl;
        return false;
    }

    // collect boundary loop
    vh = *vit;
    hh = mesh_.halfedge(vh);
    do
    {
        loop.push_back(mesh_.to_vertex(hh));
        hh = mesh_.next_halfedge(hh);
    } while (hh != mesh_.halfedge(vh));

    // map boundary loop to unit circle in texture domain
    unsigned int i, n = loop.size();
    Scalar angle, l, length;
    TexCoord t;

    // compute length of boundary loop
    for (i = 0, length = 0.0; i < n; ++i)
        length += distance(points[loop[i]], points[loop[(i + 1) % n]]);

    // map length intervalls to unit circle intervals
    for (i = 0, l = 0.0; i < n;)
    {
        // go from 2pi to 0 to preserve orientation
        angle = 2.0 * M_PI * (1.0 - l / length);

        t[0] = 0.5 + 0.5 * cosf(angle);
        t[1] = 0.5 + 0.5 * sinf(angle);

        tex[loop[i]] = t;

        ++i;
        if (i < n)
        {
            l += distance(points[loop[i]], points[loop[(i + 1) % n]]);
        }
    }

    return true;
}

//=============================================================================
