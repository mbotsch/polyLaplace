//=============================================================================
#include "Smoothing.h"
//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

void Smoothing::implicit_smoothing(Scalar timestep, unsigned int laplace,
                                   unsigned int min_point, bool rescale)
{
    if (!mesh_.n_vertices())
        return;

    // properties
    auto points = mesh_.vertex_property<Point>("v:point");

    // update Laplace matrix (if required)
    update_Laplaces(laplace, min_point);

    // store center and area
    Point center_before(0, 0, 0);
    Scalar area_before(0);
    if (rescale)
    {
        center_before = poly_centroid(mesh_);
        area_before = poly_surface_area(mesh_);
        for (auto v : mesh_.vertices())
        {
            mesh_.position(v) -= center_before;
        }
        for (auto v : mesh_.vertices())
        {
            mesh_.position(v) *= sqrt(1.0 / area_before);
        }
    }

    const unsigned int nv = mesh_.n_vertices();

    unsigned k = 0;
    SparseMatrix L, I(nv, nv), M;
    Eigen::MatrixXd B(nv, 3);

    for (auto v : mesh_.vertices())
    {
        B(k, 0) = points[v][0];
        B(k, 1) = points[v][1];
        B(k, 2) = points[v][2];
        k++;
    }

    if (laplace == 2)
    {
        if (!mesh_.has_face_property("f:point") ||
            !mesh_.has_face_property("f:weights"))
        {
            mesh_.add_face_property<pmp::Point>("f:point");
            mesh_.add_face_property<Eigen::VectorXd>("f:weights");
        }
        setup_face_point_properties(mesh_, min_point);
    }
    setup_mass_matrices(mesh_, M, laplace, min_point);
    lump_matrix(M);
    L = M - timestep * S_;

    Eigen::MatrixXd M_B = M * B;
    Eigen::MatrixXd X;
    static Eigen::SimplicialLLT<SparseMatrix> solver;
    solver.compute(L);
    X = solver.solve(M_B);

    if (solver.info() != Eigen::Success)
    {
        std::cerr << "SurfaceSmoothing: Could not solve linear system\n";
    }
    else
    {
        k = 0;
        for (unsigned int i = 0; i < nv; ++i)
        {
            Vertex v(i);
            points[v][0] = X(k, 0);
            points[v][1] = X(k, 1);
            points[v][2] = X(k, 2);
            k++;
        }
    }

    if (rescale)
    {
        Scalar area_after = poly_surface_area(mesh_);
        Scalar scale = sqrt(area_before / area_after);
        for (auto v : mesh_.vertices())
        {
            mesh_.position(v) *= scale;
        }

        Point center_after = poly_centroid(mesh_);
        Point trans = center_before - center_after;
        for (auto v : mesh_.vertices())
        {
            mesh_.position(v) += trans;
        }
    }
}