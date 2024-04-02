//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#include "Parameterization.h"
#include "diffgeo.h"
#include "LaplaceConstruction.h"

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

    const unsigned int nv = mesh_.n_vertices();

    Eigen::MatrixXd B(nv, 2);
    SparseMatrix L;
    Eigen::MatrixXd X;
    for (auto v : mesh_.vertices())
    {
        if (!mesh_.is_boundary(v))
        {
            B(v.idx(), 0) = 0.0;
            B(v.idx(), 1) = 0.0;
        }
        else
        {
            B(v.idx(), 0) = tex[v][0];
            B(v.idx(), 1) = tex[v][1];
        }
    }
    setup_stiffness_matrices(mesh_, L, laplace, min_point);

    // Adjust the right-hand-side to account for the locked nodes
    for (unsigned int i = 0; i < L.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(L, i); iter;
             ++iter)
        {
            Vertex row = pmp::Vertex(iter.row());
            Vertex col = pmp::Vertex(iter.col());
            if (!mesh_.is_boundary(row) && mesh_.is_boundary(col))
            {
                B(iter.row(), 0) -= B(iter.col(), 0) * iter.value();
                B(iter.row(), 1) -= B(iter.col(), 1) * iter.value();
            }
        }

    // Adjust the system matrix to account for the locked nodes
    for (unsigned int i = 0; i < L.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(L, i); iter;
             ++iter)
        {
            Vertex row = pmp::Vertex(iter.row());
            Vertex col = pmp::Vertex(iter.col());
            if (mesh_.is_boundary(row))
                iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
            else if (mesh_.is_boundary(col))
                iter.valueRef() = 0;
        }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
    solver.compute(L);
    X = solver.solve(B);
    if (solver.info() != Eigen::Success)
    {
        std::cerr << "harmonic(): Could not solve linear system\n";
    }
    else
    {
        // copy solution
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v))
            {
                tex[v][0] = X(v.idx(), 0);
                tex[v][1] = X(v.idx(), 1);
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
        angle = 2.0 * std::numbers::pi * (1.0 - l / length);

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
