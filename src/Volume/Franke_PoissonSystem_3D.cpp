
#include "Franke_PoissonSystem_3D.h"
#include "LaplaceConstruction_3D.h"
#include <igl/slice.h>
#include <fstream>

//=============================================================================

enum LaplaceMethods {
    Diamond = 0,
    PolySimpleLaplace = 2,
};

double franke3d(VolumeMesh::PointT vec) {
    double x, y, z;
    x = vec[0];
    y = vec[1];
    z = vec[2];

    double cx2 = (9. * x - 2.) * (9. * x - 2.);
    double cy2 = (9. * y - 2.) * (9. * y - 2.);
    double cz2 = (9. * z - 2.) * (9. * z - 2.);

    double cx1 = (9. * x + 1.) * (9. * x + 1.);
    double cx7 = (9. * x - 7.) * (9. * x - 7.);

    double cy3 = (9. * y - 3.) * (9. * y - 3.);
    double cx4 = (9. * x - 4.) * (9. * x - 4.);

    double cy7 = (9. * y - 7.) * (9. * y - 7.);
    double cz5 = (9. * z - 5.) * (9. * z - 5.);

    return (3. / 4.) *
           exp(-(1. / 4.) * cx2 - (1. / 4.) * cy2 - (1. / 4.) * cz2) +
           (3. / 4.) * exp(-(1. / 49.) * cx1 - (9. / 10.) * y - 1. / 10. -
                           (9. / 10.) * z - 1. / 10.) +
           (1. / 2.) *
           exp(-(1. / 4.) * cx7 - (1. / 4.) * cy3 - (1. / 4.) * cz5) -
           (1. / 5.) * exp(-cx4 - cy7 - cz5);
}

//-----------------------------------------------------------------------------

double laplace_franke3d(VolumeMesh::PointT vec) {
    double x, y, z;
    x = vec[0];
    y = vec[1];
    z = vec[2];

    return (243. * (-2299. + 1800. * x * (2. + 9. * x))) /
           (480200. *
            exp((9. * (12. + 10. * x * (2. + 9. * x) + 49. * y + 49. * z)) /
                490.)) -
           (486. *
            exp(-pow(4. - 9. * x, 2) - pow(7. - 9. * y, 2) -
                pow(5. - 9. * z, 2)) *
            (59. + 6. * x * (-8. + 9. * x) + 6. * y * (-14. + 9. * y) +
             6. * z * (-10 + 9 * z))) /
           5. +
           (81. *
            exp((-pow(7. - 9. * x, 2) - 9 * pow(1. - 3. * y, 2) -
                 pow(5. - 9. * z, 2)) /
                4.) *
            (77. + 9. * x * (-14. + 9. * x) + 27. * y * (-2. + 3. * y) +
             9 * z * (-10. + 9. * z))) /
           8. +
           (729. * (2. + 3. * x * (-4. + 9. * x) + 3. * y * (-4. + 9. * y) +
                    3. * z * (-4. + 9. * z))) /
           (16. *
            exp((3. * (4. + 3. * x * (-4. + 9. * x) +
                       3. * y * (-4. + 9. * y) + 3. * z * (-4. + 9. * z))) /
                4.0));
}

//-----------------------------------------------------------------------------

double solve_franke_poisson(VolumeMesh &mesh, int Laplace, int face_point,
                            int cell_point) {


    Eigen::SparseMatrix<double> M, S, S_f;
    Eigen::VectorXd b(mesh.n_vertices());

    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    setup_3D_stiffness_matrix(mesh, S, Laplace, face_point, cell_point);
    setup_3D_mass_matrix(mesh, M, Laplace, face_point, cell_point);

    for (auto v: mesh.vertices()) {
        b(v.idx()) = laplace_franke3d(mesh.vertex(v));
    }

    b = M * b;

    // Set the constraints at the locked vertices to the evluation of the Franke function
    for (auto v: mesh.vertices()) {
        if (mesh.is_boundary(v)) {
            // right-hand side: fix boundary values with franke function of the vertices
            b(v.idx()) = franke3d(mesh.vertex(v));
        }
    }

    // Adjust the right-hand-side to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter; ++iter) {
            OpenVolumeMesh::VertexHandle row = OpenVolumeMesh::VertexHandle((int) iter.row());
            OpenVolumeMesh::VertexHandle col = OpenVolumeMesh::VertexHandle((int) iter.col());
            if (!mesh.is_boundary(row) && mesh.is_boundary(col)) {
                b[iter.row()] -= b[iter.col()] * iter.value();
            }
        }

    // Adjust the system matrix to account for the locked nodes
    for (unsigned int i = 0; i < S.outerSize(); i++)
        for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i); iter; ++iter) {
            OpenVolumeMesh::VertexHandle row = OpenVolumeMesh::VertexHandle((int) iter.row());
            OpenVolumeMesh::VertexHandle col = OpenVolumeMesh::VertexHandle((int) iter.col());
            if (mesh.is_boundary(row)) iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
            else if (mesh.is_boundary(col)) iter.valueRef() = 0;
        }

    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

    solver.compute(S);
    Eigen::VectorXd x = solver.solve(b);
    std::cout << "Size x :" << x.size() << std::endl;
    double error = 0.0;
    for (auto v: mesh.vertices()) {
//            std::cout << "x: " << x[v.idx()] << " Franke : " <<  franke3d( mesh.vertex(v) )<<std::endl;
        error += pow(x[v.idx()] - franke3d(mesh.vertex(v)), 2.);
    }

    std::cout << "DoF " << mesh.n_vertices() << std::endl;
    std::cout << "Franke RMSE error inner vertices: "
              << sqrt(error / (double) mesh.n_vertices()) << std::endl;
    return sqrt(error / (double) mesh.n_vertices());

}

