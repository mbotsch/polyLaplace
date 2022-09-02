#include "LaplaceConstruction_3D.h"
#include "diffgeo_3D.h"
#include "../Surface/diffgeo.h"
#include <Eigen/Sparse>

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================
enum LaplaceMethods {
    Diamond = 0,
    PolySimpleLaplace = 2,
};

void setup_3D_stiffness_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &S, int Laplace, int face_point,
                               int cell_point) {
    if (Laplace == Diamond) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        SparseMatrix G, V, Div, P, Pc, Pf;
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_diamond_gradient(mesh, G, V);
        Div = -G.transpose() * V;
        S = P.transpose() * Div * G * P;

    } else if (Laplace == PolySimpleLaplace) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        SparseMatrix Sf, G, Grad, V, Div, P, Pc, Pf;
        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_polysimple_stiffness_matrix(mesh, S, Sf, face_point, cell_point);
    }
}

//-----------------------------------------------------------------------------
void
setup_3D_mass_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &M, int Laplace, int face_point,
                     int cell_point) {
    if (Laplace == Diamond) {
        auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
        auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");
        auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
        auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");

        compute_3D_virtual_points_and_weights(mesh, face_point, cell_point);
        SparseMatrix P, Pc, Pf, M2;
        setup_3D_cell_face_prolongation_matrix(mesh, P, Pc, Pf);
        setup_3D_diamond_mass_matrix(mesh, M2);
        M = P.transpose() * M2 * P;
        M /= 2.0;

    } else if (Laplace == PolySimpleLaplace) {
        setup_3D_polysimple_mass_matrix(mesh, M, face_point, cell_point);
        lump_matrix(M);
    }
    double count = 0;
    for (int k = 0; k < M.outerSize(); ++k) {
        for (SparseMatrix::InnerIterator it(M, k); it; ++it) {
            count += it.value();
        }
    }
    std::cout << "Volume mass matrix: " << count << std::endl;
}


//-----------------------------------------------------------------------------

void setup_3D_cell_face_prolongation_matrix(VolumeMesh &mesh, Eigen::SparseMatrix<double> &P,
                                            Eigen::SparseMatrix<double> &P_cells,
                                            Eigen::SparseMatrix<double> &P_faces) {

    std::vector<T> triplet_face;
    std::vector<T> triplet_cells;

    auto f_w_prop = mesh.request_face_property<Eigen::VectorXd>("face weights");
    auto c_w_prop = mesh.request_cell_property<Eigen::VectorXd>("cell weights");


    int n_v = (int) mesh.n_vertices();
    int n_f = (int) mesh.n_faces();
    int n_c = (int) mesh.n_cells();


    for (auto v_it = mesh.v_iter(); v_it.valid(); ++v_it) {
        triplet_cells.emplace_back(T(v_it->idx(), v_it->idx(), 1.0));
        triplet_face.emplace_back(T(v_it->idx(), v_it->idx(), 1.0));
    }

    for (auto f_it = mesh.f_iter(); f_it.valid(); ++f_it) {
        Eigen::VectorXd face_w = f_w_prop[*f_it];
        auto f_v_it_pair = mesh.face_vertices(*f_it);
        int i = 0;
        for (auto f_v_it = f_v_it_pair.first;
             f_v_it != f_v_it_pair.second; ++f_v_it) {

            //add weight entries of the virtual face vertices in first prolongation matrix
            triplet_face.emplace_back(T(n_v + f_it->idx(), f_v_it->idx(), face_w(i)));
            i++;
        }
        //second prolongation matrix considers virtual face vertices as existent in the mesh
        triplet_cells.emplace_back(T(n_v + f_it->idx(), n_v + f_it->idx(), 1.0));
    }

    for (auto c_it = mesh.c_iter(); c_it.valid(); ++c_it) {

        Eigen::VectorXd cell_w = c_w_prop[*c_it];
        //face vertices are pushed back first during computation, so their weights make up the first part of the vector
        auto c_f_it_pair = mesh.cell_faces(*c_it);
        auto c_v_it_pair = mesh.cell_vertices(*c_it);

        int i = 0;
        for (auto c_f_it = c_f_it_pair.first;
             c_f_it != c_f_it_pair.second; ++c_f_it) {

            triplet_cells.emplace_back(T(n_v + n_f + c_it->idx(), n_v + c_f_it->idx(), cell_w(i)));
            i++;
        }
        for (auto c_v_it = c_v_it_pair.first;
             c_v_it != c_v_it_pair.second; ++c_v_it) {

            triplet_cells.emplace_back(T(n_v + n_f + c_it->idx(), c_v_it->idx(), cell_w(i)));
            i++;
        }

    }

    P_faces.resize(n_v + n_f, n_v);
    P_cells.resize(n_v + n_f + n_c, n_v + n_f);
    P_faces.setFromTriplets(triplet_face.begin(), triplet_face.end());
    P_cells.setFromTriplets(triplet_cells.begin(), triplet_cells.end());
    P = P_cells * P_faces;
}
//-----------------------------------------------------------------------------


double
volume(const Eigen::Vector3d &a, const Eigen::Vector3d &b, const Eigen::Vector3d &c, const Eigen::Vector3d &d) {
    Eigen::Matrix3d A;
    A.col(0) = b - a;
    A.col(1) = c - a;
    A.col(2) = d - a;

    return A.determinant() / 6.;
}

