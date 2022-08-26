
#include "Franke_PoissonSystem_3D.h"
#include "LaplaceConstruction_3D.h"
#include "HarmonicBasis.h"
#include <igl/slice.h>
#include <Eigen/CholmodSupport>
#include <fstream>

//=============================================================================

enum LaplaceMethods
{
    Diamond = 0,
    Dual_Laplace = 2,
    Harmonic= 3,
    Sandwich = 4,
    AQAPoly_MG = 5,
    AQAPoly = 6
};

double franke3d(VolumeMesh::PointT vec)
{
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

double laplace_franke3d(VolumeMesh::PointT vec)
{
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
                            int cell_point, int degree, std::string meshname,
                            CoarseDimension dof)
{

    if (Laplace == AQAPoly)
    {
        return solve_3D_AQAPoly_Poisson(meshname, dof,degree);
    }else if(Laplace == AQAPoly_MG)
    {

        std::ofstream file_timings("proxy.csv");
        return solve_3D_AQAPoly_Poisson_mg(meshname, file_timings, dof, degree);
    }else if(Laplace == Harmonic){
        return solve_3D_Franke_harmonic(meshname);
    }else
    {
        Eigen::SparseMatrix<double> M, S, S_f;
        Eigen::VectorXd b(mesh.n_vertices());

            auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
            auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
            setup_3D_stiffness_matrix(mesh, S, Laplace, face_point, cell_point,degree, meshname);
            setup_3D_mass_matrix(mesh, M, Laplace, face_point, cell_point,degree, meshname);

        for (auto v : mesh.vertices())
        {
            b(v.idx()) = laplace_franke3d(mesh.vertex(v));
        }

        b = M * b;

        // Set the constraints at the locked vertices to the evluation of the Franke function
        for (auto v : mesh.vertices()){
            if (mesh.is_boundary(v))
            {
                // right-hand side: fix boundary values with franke function of the vertices
                b(v.idx()) = franke3d(mesh.vertex(v));
            }
        }

        // Adjust the right-hand-side to account for the locked nodes
        for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter )
            {
                OpenVolumeMesh::VertexHandle row =  OpenVolumeMesh::VertexHandle(iter.row());
                OpenVolumeMesh::VertexHandle col =  OpenVolumeMesh::VertexHandle(iter.col());
                if( !mesh.is_boundary(row) && mesh.is_boundary(col) ){
                    b[ iter.row() ] -= b[ iter.col() ] * iter.value();
                }
            }

//        std::cout << b << std::endl;
        // Adjust the system matrix to account for the locked nodes
        for( unsigned int i=0 ; i<S.outerSize() ; i++ ) for( Eigen::SparseMatrix< double >::InnerIterator iter(S,i) ; iter ; ++iter ){
                OpenVolumeMesh::VertexHandle row =  OpenVolumeMesh::VertexHandle(iter.row());
                OpenVolumeMesh::VertexHandle col =  OpenVolumeMesh::VertexHandle(iter.col());
                if( mesh.is_boundary(row) ) iter.valueRef() = iter.row()==iter.col() ? 1. : 0.;
                else if( mesh.is_boundary(col)) iter.valueRef() = 0;
        }


        Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;
//        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > solver;

        solver.compute(S);
        Eigen::VectorXd x = solver.solve( b );
        std::cout << "Size x :" <<  x.size() <<std::endl;
        double error = 0.0;
        for (auto v : mesh.vertices()){
//            std::cout << "x: " << x[v.idx()] << " Franke : " <<  franke3d( mesh.vertex(v) )<<std::endl;
            error += pow( x[v.idx()] - franke3d( mesh.vertex(v) ) , 2. );
        }

        std::cout << "DoF " << mesh.n_vertices() << std::endl;
        std::cout << "Franke RMSE error inner vertices: "
                  << sqrt(error / (double)mesh.n_vertices()) << std::endl;
        return sqrt(error / (double)mesh.n_vertices()) ;
    }
}

//-----------------------------------------------------------------------------

void solve_laplace_equation(VolumeMesh &mesh, int laplace, int face_point,
                            int cell_point, int degree)
{
    Eigen::SparseMatrix<double> M, S, S_f;
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(mesh.n_vertices(), 3);

    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    if (laplace == AQAPoly && degree == 2)
    {
        B = Eigen::MatrixXd::Zero(mesh.n_vertices() + mesh.n_edges(), 3);
    }

    setup_3D_stiffness_matrix(mesh, S, laplace, face_point, cell_point, degree);
    int nb = 0;

//-------------------------------------------------
//    if(laplace == Sandwich){
//        auto v_prop = mesh.request_vertex_property<bool>("v:boundary");
//
//        for(auto c: mesh.cells()){
//            for (auto v : mesh.cell_vertices(c))
//            {
//                if(mesh.is_boundary(c)){
//                    v_prop[v] =true;
//                }else{
//                    v_prop[v] =false;
//                }
//            }
//        }
//
//        for (auto v : mesh.vertices())
//        {
//            //count nr outer vertices
//            if (v_prop[v])
//            {
//                nb++;
//            }
//        }
//
//        Eigen::SparseMatrix<double> G, V, Div, P, Pc, Pf;
//
//        unsigned ins = 0;
//        unsigned out = 0;
//
//        Eigen::VectorXi in(mesh.n_vertices() - nb), bound(nb), x(3);
//        x << 0, 1, 2;
//
//        for (auto v : mesh.vertices())
//        {
//            // save indices of inner and outer vertices
//
//            if (!v_prop[v])
//            {
//                in(ins) = v.idx();
//                ins++;
//            }
//            else
//            {
//                bound(out) = v.idx();
//                out++;
//                // right-hand side: fix x coordinate of boundary vertices for the righthandsite
//                B(v.idx(), 0) = mesh.vertex(v)[0];
//                B(v.idx(), 1) = mesh.vertex(v)[1];
//                B(v.idx(), 2) = mesh.vertex(v)[2];
//            }
//        }
//
//        Eigen::SparseMatrix<double> L_in_in, L_in_b;
//        Eigen::MatrixXd b_in, b_out;
//
//        // slice S and b and bring boundary values on the righthandsite to solve only for inner vertices
//        igl::slice(S, in, in, L_in_in);
//        igl::slice(S, in, bound, L_in_b);
//        igl::slice(B, in, x, b_in);
//        igl::slice(B, bound, x, b_out);
//        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
//        solver.compute(L_in_in);
//        Eigen::MatrixXd X = solver.solve(b_in - L_in_b * b_out);
//        double error = 0.0;
//
//        int k = 0;
//        if (solver.info() != Eigen::Success)
//        {
//            std::cerr << "harmonic(): Could not solve linear system\n";
//        }
//        else
//        {
//            // copy solution
//            for (auto v : mesh.vertices())
//            {
//                if (!v_prop[v])
//                {
//                    error += pow(X(k, 0) - mesh.vertex(v)[0], 2) +
//                             pow(X(k, 1) - mesh.vertex(v)[1], 2) +
//                             pow(X(k, 2) - mesh.vertex(v)[2], 2);
//                    //                std::cout << "   Point computed: " << X.row(k) << "  Point mesh: " << mesh.vertex(v) <<std::endl;
//                    k++;
//                }
//            }
//            std::cout << "RMSE inner verticex positions : " << sqrt(error / (double)k) << std::endl;
//        }
//    }
////---------------------------------------------------
//    else
//    {

        for (auto v : mesh.vertices())
        {
            //count nr outer vertices
            if (mesh.is_boundary(v))
            {
                nb++;
            }
        }

        if (laplace == AQAPoly && degree == 2)
        {
            for (auto e : mesh.edges())
            {
                if (mesh.is_boundary(e))
                {
                    nb++;
                }
            }
        }
        Eigen::SparseMatrix<double> G, V, Div, P, Pc, Pf;

        unsigned ins = 0;
        unsigned out = 0;

        Eigen::VectorXi in(mesh.n_vertices() - nb), bound(nb), x(3);
        x << 0, 1, 2;

        if (laplace == AQAPoly && degree == 2)
        {
            in.resize(mesh.n_vertices() + mesh.n_edges() - nb);
        }

        for (auto v : mesh.vertices())
        {
            // save indices of inner and outer vertices

            if (!mesh.is_boundary(v))
            {
                in(ins) = v.idx();
                ins++;
            }
            else
            {
                bound(out) = v.idx();
                out++;
                // right-hand side: fix x coordinate of boundary vertices for the righthandsite
                B(v.idx(), 0) = mesh.vertex(v)[0];
                B(v.idx(), 1) = mesh.vertex(v)[1];
                B(v.idx(), 2) = mesh.vertex(v)[2];
            }
        }
        if (laplace == AQAPoly && degree == 2)
        {
            for (auto e : mesh.edges())
            {
                if (!mesh.is_boundary(e))
                {
                    in(ins) = mesh.n_vertices() + e.idx();
                    ins++;
                }
                else
                {
                    bound(out) = mesh.n_vertices() + e.idx();
                    out++;
                    VolumeMesh::PointT edge_mid =
                        0.5 * (mesh.vertex(mesh.edge_vertices(e)[0]) +
                               mesh.vertex(mesh.edge_vertices(e)[1]));
                    // right-hand side: fix boundary values with franke function of the vertices
                    B(mesh.n_vertices() + e.idx(), 0) = edge_mid[0];
                    B(mesh.n_vertices() + e.idx(), 1) = edge_mid[1];
                    B(mesh.n_vertices() + e.idx(), 2) = edge_mid[2];
                }
            }
        }

        Eigen::SparseMatrix<double> L_in_in, L_in_b;
        Eigen::MatrixXd b_in, b_out;

        // slice S and b and bring boundary values on the righthandsite to solve only for inner vertices
        igl::slice(S, in, in, L_in_in);
        igl::slice(S, in, bound, L_in_b);
        igl::slice(B, in, x, b_in);
        igl::slice(B, bound, x, b_out);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(L_in_in);
        Eigen::MatrixXd X = solver.solve(b_in - L_in_b * b_out);
        double error = 0.0;

        int k = 0;
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "harmonic(): Could not solve linear system\n";
        }
        else
        {
            // copy solution
            for (auto v : mesh.vertices())
            {
                if (!mesh.is_boundary(v))
                {
                    error += pow(X(k, 0) - mesh.vertex(v)[0], 2) +
                             pow(X(k, 1) - mesh.vertex(v)[1], 2) +
                             pow(X(k, 2) - mesh.vertex(v)[2], 2);
                    //                std::cout << "   Point computed: " << X.row(k) << "  Point mesh: " << mesh.vertex(v) <<std::endl;
                    k++;
                }
            }
            if (laplace == AQAPoly && degree == 2)
            {
                for (auto e : mesh.edges())
                {
                    if (!mesh.is_boundary(e))
                    {
                        VolumeMesh::PointT edge_mid =
                            0.5 * (mesh.vertex(mesh.edge_vertices(e)[0]) +
                                   mesh.vertex(mesh.edge_vertices(e)[1]));
                        // right-hand side: fix boundary values with franke function of the vertices
                        error += pow(X(k, 0) - edge_mid[0], 2) +
                                 pow(X(k, 1) - edge_mid[1], 2) +
                                 pow(X(k, 2) - edge_mid[2], 2);
                        //                    std::cout << "solved: " << X.row(k)
                        //                              << " edge vertex :" << edge_mid << std::endl;

                        k++;
                    }
                }
            }

            std::cout << "RMSE inner verticex positions : "
                      << sqrt(error / (double)k) << std::endl;
        }
//    }
}
