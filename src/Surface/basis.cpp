
#include "basis.h"

//----------------------------------------------------------------------------------


void linear_basis_gradient(Eigen::MatrixXd &basis) {
    basis.resize(3, 2);
    Eigen::Vector2d g_phi1(-1, -1);
    Eigen::Vector2d g_phi2(1, 0);
    Eigen::Vector2d g_phi3(0, 1);

    basis.row(0) = g_phi1;
    basis.row(1) = g_phi2;
    basis.row(2) = g_phi3;
}
//----------------------------------------------------------------------------------

void grad_phi_1(Eigen::Vector2d &val, double x, double y) {
    val[0] = 4.0 * x + 4.0 * y - 3.0;
    val[1] = 4.0 * x + 4.0 * y - 3.0;
}

void grad_phi_2(Eigen::Vector2d &val, double x, double y) {
    val[0] = 4.0 * x - 1.0;
    val[1] = 0.0;
}

void grad_phi_3(Eigen::Vector2d &val, double x, double y) {
    val[0] = 0.0;
    val[1] = 4.0 * y - 1.0;
}

void grad_phi_4(Eigen::Vector2d &val, double x, double y) {
    val[0] = -8.0 * x - 4.0 * y + 4.0;
    val[1] = -4.0 * x;
}

void grad_phi_5(Eigen::Vector2d &val, double x, double y) {
    val[0] = 4.0 * y;
    val[1] = 4.0 * x;
}

void grad_phi_6(Eigen::Vector2d &val, double x, double y) {
    val[0] = -4.0 * y;
    val[1] = -8.0 * y - 4.0 * x + 4.0;
}

//----------------------------------------------------------------------------------

int local_to_global_index(pmp::SurfaceMesh &mesh, pmp::Face &f, pmp::Halfedge &he, int i) {
    int idx, idxu3, idxu5;

    auto halfedge_idx = mesh.get_halfedge_property<int>("he:idx");

    pmp::Edge e = mesh.edge(he);


    pmp::Vertex v = mesh.vertex(e, 1);
    pmp::Vertex vv = mesh.vertex(e, 0);
    pmp::Halfedge next_he = mesh.next_halfedge(he);
    if ((mesh.from_vertex(he).idx() == v.idx()) &&
        (mesh.to_vertex(he).idx() == vv.idx())) {
        idxu3 = halfedge_idx[he];
        idxu5 = halfedge_idx[next_he];
    } else {
        idxu3 = halfedge_idx[next_he];
        idxu5 = halfedge_idx[he];
    }

    if (i == 0) {
        idx = mesh.n_vertices() + f.idx();
    } else if (i == 1) {
        idx = v.idx();
    } else if (i == 2) {
        idx = vv.idx();
    } else if (i == 3) {
        idx = mesh.n_vertices() + mesh.n_faces() + mesh.n_edges() + idxu3;
    } else if (i == 4) {
        idx = mesh.n_vertices() + mesh.n_faces() + e.idx();
    } else {
        idx = mesh.n_vertices() + mesh.n_faces() + mesh.n_edges() + idxu5;
    }
    return idx;
}

//----------------------------------------------------------------------------------

double gauss_quadrature_gradient(Eigen::MatrixXd &jacobian, int bi, int bj) {
    double val = 0;
    Eigen::Vector2d gi, gj;

    Eigen::Vector3d gauss_coordinates(1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0);
    //    Eigen::Vector3d gauss_coordinates(0.5, 0.0, 0.5);

    for (int k = 0; k < 3; k++) {
        //        for(int j = 0 ; j<3; j++){
        double x = gauss_coordinates(k);
        int yidx = (k + 2) % 3;
        double y = gauss_coordinates(yidx);
        compute_basis_gradient(gi, bi, x, y);

        compute_basis_gradient(gj, bj, x, y);
        val += gi.transpose() * jacobian * jacobian.transpose() * gj;
    }

    // 1/6 ist gaussian weighting term for all three points
    val /= 6.0;
    return val;
}
//----------------------------------------------------------------------------------

double gauss_quadrature_mass(int bi, int bj) {

    //--------------------Quintic Gauss Quadrature----------------------
    double val = 0.0;
    Eigen::VectorXd gauss_coordinates_x(7), gauss_coordinates_y(7);
    gauss_coordinates_x << 1.0 / 3.0, 0.0597158717, 0.4701420641, 0.4701420641,
            0.7974269853, 0.1012865073, 0.1012865073;
    gauss_coordinates_y << 1.0 / 3.0, 0.4701420641, 0.0597158717, 0.4701420641,
            0.1012865073, 0.7974269853, 0.1012865073;

    bool linear = true;
    for (int k = 0; k < 7; k++)
    {
        double x = gauss_coordinates_x(k);
        double y = gauss_coordinates_y(k);
        double fi = compute_triangle_basis(bi, x, y, !linear);
        double fj = compute_triangle_basis(bj, x, y, !linear);
        if (k == 0)
            val += fi * fj * 0.225*0.5;
        else if( k < 4)
        {
            val += fi * fj * 0.1323941527*0.5;
        }else{
            val += fi * fj * 0.1259391805*0.5;
        }
    }

    return val;
}
//----------------------------------------------------------------------------------

void compute_basis_gradient(Eigen::Vector2d &grad, int idx, double x, double y) {
    if (idx == 0) {
        grad_phi_1(grad, x, y);
    } else if (idx == 1) {
        grad_phi_2(grad, x, y);
    } else if (idx == 2) {
        grad_phi_3(grad, x, y);
    } else if (idx == 3) {
        grad_phi_4(grad, x, y);
    } else if (idx == 4) {
        grad_phi_5(grad, x, y);
    } else {
        grad_phi_6(grad, x, y);
    }
}
//----------------------------------------------------------------------------------

double compute_triangle_basis(int basis_idx, double x, double y, bool lin) {
    double val;
    if (lin) {
        if (basis_idx == 0) {
            val = phi_1_lin(x, y);
        } else if (basis_idx == 1) {
            val = phi_2_lin(x, y);
        } else if (basis_idx == 2) {
            val = phi_3_lin(x, y);
        }
    } else {
        if (basis_idx == 0) {
            val = phi_1_quad(x, y);
        } else if (basis_idx == 1) {
            val = phi_2_quad(x, y);
        } else if (basis_idx == 2) {
            val = phi_3_quad(x, y);
        } else if (basis_idx == 3) {
            val = phi_4_quad(x, y);
        } else if (basis_idx == 4) {
            val = phi_5_quad(x, y);
        } else {
            val = phi_6_quad(x, y);
        }
    }
    return val;
}

//----------------------------------------------------------------------------------
double compute_quad_basis(int basis_idx, double x, double y, bool lin) {
    // Followed this https://www.mathematik.uni-dortmund.de/~kuzmin/cfdintro/lecture7.pdf
    double val;
    if (lin) {
        if (basis_idx == 0) {
            //phi 1 linear
            val = 1.0 - x - y + x * y;
        } else if (basis_idx == 1) {
            val = x - x * y;
        } else if (basis_idx == 2) {
            val = x * y;
        } else if (basis_idx == 3) {
            val = y - x * y;
        }
    } else {
        val = 0;
    }
    return val;
}
//----------------------------------------------------------------------------------

double phi_1_lin(double x, double y) {
    return 1.0 - x - y;
}

double phi_2_lin(double x, double y) {
    return x;
}

double phi_3_lin(double x, double y) {
    return y;
}

//----------------------------------------------------------------------------------

double phi_1_quad(double x, double y) {
    return 2.0 * x * x + 4.0 * x * y - 3.0 * x + 2.0 * y * y - 3.0 * y + 1.0;
}

double phi_2_quad(double x, double y) {
    return 2.0 * x * x - x;
}

double phi_3_quad(double x, double y) {
    return 2.0 * y * y - y;
}

double phi_4_quad(double x, double y) {
    return -4.0 * x * x - 4.0 * x * y + 4.0 * x;
}

double phi_5_quad(double x, double y) {
    return 4.0 * x * y;
}

double phi_6_quad(double x, double y) {
    return -4.0 * x * y - 4.0 * y * y + 4.0 * y;
}

//----------------------------------------------------------------------------------
