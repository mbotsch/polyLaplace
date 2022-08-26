//=============================================================================

#include "Smoothing.h"
#include <igl/slice.h>
//=============================================================================

using namespace pmp;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

void Smoothing::implicit_smoothing_misha(Scalar timestep, unsigned int laplace, unsigned int min_point)
{

    if (!mesh_.n_vertices())
        return;

    // properties
    auto points = mesh_.vertex_property<Point>("v:point");
    auto area_points = mesh_.add_face_property<Point>("f:point");
    auto area_weights = mesh_.add_face_property<Eigen::VectorXd>("f:weights");

    setup_face_point_properties(mesh_, min_point);
    // update Laplace matrix (if required)
    update_Laplaces(laplace,min_point);

    // store center and area
    Point center_before = my_centroid(mesh_);
    Scalar area_before = my_surface_area(mesh_);

    for (auto v : mesh_.vertices()) {
        mesh_.position(v) -= center_before;
    }
    for (auto v : mesh_.vertices()) {
        mesh_.position(v) *= sqrt(1.0 / area_before);
    }
    // Compute the chosen implicit points and its convex combination weights for each face
    setup_face_point_properties(mesh_, min_point);


    const unsigned int nv = mesh_.n_vertices();

    unsigned k = 0;
    SparseMatrix L, I(nv, nv), M;
    Eigen::MatrixXd B(nv, 3);

    I.setIdentity();

    for (auto v : mesh_.vertices()) {
        B(k, 0) = points[v][0];
        B(k, 1) = points[v][1];
        B(k, 2) = points[v][2];
        k++;
    }

    setup_mass_matrices(mesh_,M,laplace, min_point,true);
    L = M - timestep * S_;

    Eigen::MatrixXd M_B = M * B;
    SparseMatrix L_in_in, L_in_b;
    Eigen::MatrixXd b_in, b_out, X;
    static Eigen::SimplicialLLT<SparseMatrix> solver;
    solver.compute(L);
    X = solver.solve(M_B);

    if (solver.info() != Eigen::Success) {
        std::cerr << "SurfaceSmoothing: Could not solve linear system\n";
    } else {
        //                // copy solution
        k = 0;
        for (unsigned int i = 0; i < nv; ++i) {
            Vertex v(i);
            if (!mesh_.is_boundary(v)) {
                points[v][0] = X(k, 0);
                points[v][1] = X(k, 1);
                points[v][2] = X(k, 2);
                k++;
            }
        }
    }

    //update point properties for rescaling
    setup_face_point_properties(mesh_, min_point);

    //     restore original surface area
    Scalar area_after = my_surface_area(mesh_);

    Scalar scale = sqrt(1 / area_after);
    Point center_after = my_centroid(mesh_);

    for (auto v : mesh_.vertices()) {
        mesh_.position(v) -= center_after;
    }


    for (auto v : mesh_.vertices()) {
        mesh_.position(v) *= scale;
    }

    mesh_.remove_face_property<Point>(area_points);
    mesh_.remove_face_property<Eigen::VectorXd>(area_weights);

}