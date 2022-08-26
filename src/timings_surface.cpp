#include <iostream>
#include <pmp/Timer.h>
#include <unsupported/Eigen/SparseExtra>
#include "Surface/PolyLaplace.h"
#include "Surface/LaplaceConstruction.h"
#include "Surface/DisneyLaplace.h"
#include "Surface/SECLaplace.h"
#include "Surface/diffgeo.h"
#include <Eigen/CholmodSupport>


enum LaplaceMethods {
    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    IntrinsicDelaunay = 4,
    Disney = 5,
    SEC = 6,
    AQAPoly_Laplace = 7,
    quadratic_Triangle_Laplace = 8
};

enum InsertedPoint {
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2
};


using namespace std;


void takeTimingLGS(pmp::SurfaceMesh mesh, vector<double> &compute, vector<double> &solve, int laplace, int minpoint,
                   int solver_type) {

    auto points = mesh.vertex_property<pmp::Point>("v:point");
    Eigen::SparseMatrix<double> L, M, A;
    if (laplace == AlexaLaplace) {
        poly_laplace_lambda_ = 2.0;
        setup_stiffness_matrices(mesh, L, AlexaLaplace);
        setup_mass_matrices(mesh, M, AlexaLaplace);
    } else if (laplace == Disney) {
        disney_laplace_lambda_ = 1.0;
        setup_stiffness_matrices(mesh, L, Disney);
        setup_mass_matrices(mesh, M, Disney);
    } else if (laplace == IntrinsicDelaunay || laplace == CotanLaplace) {
        if (!mesh.is_triangle_mesh()) {
            SurfaceTriangulation tesselator(mesh);
            tesselator.triangulate(SurfaceTriangulation::Objective::MIN_AREA);
        }
        setup_stiffness_matrices(mesh, L, laplace, minpoint);
        setup_mass_matrices(mesh, M, laplace, minpoint);
    } else if (laplace == SEC) {
        sec_laplace_lambda_ = 2.0;
        setup_stiffness_matrices(mesh, L, SEC);
        setup_mass_matrices(mesh, M, SEC);
    } else {
        setup_stiffness_matrices(mesh, L, laplace, minpoint);
        setup_mass_matrices(mesh, M, laplace, minpoint);
    }

    Eigen::MatrixXd B(mesh.n_vertices(), 3);
    if (laplace == AQAPoly_Laplace || laplace == quadratic_Triangle_Laplace) {
        B.resize(mesh.n_vertices() + mesh.n_edges(), 3);
    }

    for (auto v : mesh.vertices()) {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }

    if (laplace == AQAPoly_Laplace || laplace == quadratic_Triangle_Laplace) {
        int nv = mesh.n_vertices();
        for (auto e : mesh.edges()) {
            Point mid_point = 0.5 * (mesh.position(mesh.vertex(e, 0)) + mesh.position(mesh.vertex(e, 1)));
            B(nv + e.idx(), 0) = mid_point[0];
            B(nv + e.idx(), 1) = mid_point[1];
            B(nv + e.idx(), 2) = mid_point[2];
        }
    }

    A = M - 0.1 * L;
    Eigen::MatrixXd M_B = M * B;


    if (solver_type == 0) {


        static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        pmp::Timer t;
        solver.compute(A);

        t.start();
        for (int i = 0; i < 20; i++) {
            solver.compute(A);
        }
        t.stop();

        std::cout << "Factorization of system matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
        compute.emplace_back(t.elapsed() / 20.0);

        Eigen::MatrixXd X = solver.solve(M_B);

        t.start();
        for (int i = 0; i < 20; i++) {
            Eigen::MatrixXd X = solver.solve(M_B);
        }
        t.stop();
        std::cout << "Time for solving CMCf : " << t.elapsed() / 20.0 << "ms" << std::endl;
        solve.emplace_back(t.elapsed() / 20.0);
    } else if (solver_type == 1) {
        //------------Cholmod----------------------

        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > cd;
        pmp::Timer t;
        cd.compute(A);

        t.start();
        for (int i = 0; i < 20; i++) {
            cd.compute(A);
        }
        t.stop();
        std::cout << "Factorization of system matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
        compute.emplace_back(t.elapsed() / 20.0);
        Eigen::MatrixXd X = cd.solve(M_B);

        t.start();
        for (int i = 0; i < 20; i++) {
            Eigen::MatrixXd X = cd.solve(M_B);
        }
        t.stop();
        std::cout << "Time for solving CMCf : " << t.elapsed() / 20.0 << "ms" << std::endl;
        solve.emplace_back(t.elapsed() / 20.0);

//    //-----------------------------------------
    } else {
        Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > cd;
        pmp::Timer t;
        cd.compute(A);

        t.start();
        for (int i = 0; i < 20; i++) {
            cd.compute(A);
        }
        t.stop();
        std::cout << "Factorization of system matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
        compute.emplace_back(t.elapsed() / 20.0);

        Eigen::MatrixXd X = cd.solve(M_B);

        t.start();
        for (int i = 0; i < 20; i++) {
            Eigen::MatrixXd X = cd.solve(M_B);
        }
        t.stop();
        std::cout << "Time for solving CMCf : " << t.elapsed() / 20.0 << "ms" << std::endl;
        solve.emplace_back(t.elapsed() / 20.0);
    }
}

double stiffness_sparsity(pmp::SurfaceMesh &mesh, unsigned int laplace, unsigned int min_point = 0) {

    if (laplace == SandwichLaplace) {
        std::cout << "Sandwich Laplace : \n";
    } else if (laplace == AlexaLaplace) {
        std::cout << "Alexa Laplace : \n";
    } else if (laplace == CotanLaplace) {
        std::cout << "Cotan Laplace : \n";
    } else if (laplace == Diamond) {
        std::cout << "Diamond Laplace : \n";
    } else if (laplace == Disney) {
        std::cout << "Disney Laplace : \n";
    } else if (laplace == SEC) {
        std::cout << " SEC Laplace : \n";
    } else if (laplace == AQAPoly_Laplace) {
        std::cout << " AQAPoly Laplace : \n";
    }
    Eigen::SparseMatrix<double> S, M;
    setup_stiffness_matrices(mesh, S, laplace, min_point);

//    double val = round(((double) sparsity(S) / (S.cols() * S.rows()) * 100.0) * 1000) / 1000;
    double val = S.nonZeros();
    double dim = S.rows() * S.cols();
    std::cout << "Entries: " << dim << " non Zeros: " << val << std::endl;

    return val ;
}

void compute_sparsity() {
    std::string filename_;
    filename_ = "./surface_sparsity.csv";
    std::ofstream file_sparse(filename_);
    std::vector<double> compute_alexa, compute_diamond, compute_sandwich, compute_disney, compute_SubdivTri, compute_AQAPoly;
    file_sparse << "Laplace,Hexas,Triangles,Concave,Quads\n";

    SurfaceMesh mesh;

    mesh.read("../data/fernando_meshes/sphere/hex_6.obj");
    compute_alexa.emplace_back(stiffness_sparsity(mesh, AlexaLaplace));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, AreaMinimizer));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, SandwichLaplace, AreaMinimizer));
    compute_disney.emplace_back(stiffness_sparsity(mesh, Disney));
    compute_AQAPoly.emplace_back(stiffness_sparsity(mesh, AQAPoly_Laplace, AreaMinimizer));
    insert_points(mesh,AreaMinimizer);
    compute_SubdivTri.emplace_back(stiffness_sparsity(mesh, quadratic_Triangle_Laplace));

    mesh.read("../data/fernando_meshes/sphere/triangle_6.obj");
    compute_alexa.emplace_back(stiffness_sparsity(mesh, AlexaLaplace));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, AreaMinimizer));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, SandwichLaplace, AreaMinimizer));
    compute_disney.emplace_back(stiffness_sparsity(mesh, Disney));
    compute_AQAPoly.emplace_back(stiffness_sparsity(mesh, AQAPoly_Laplace, AreaMinimizer));
    insert_points(mesh,AreaMinimizer);
    compute_SubdivTri.emplace_back(stiffness_sparsity(mesh, quadratic_Triangle_Laplace));

    mesh.read("../data/fernando_meshes/sphere/concave_6.obj");
    compute_alexa.emplace_back(stiffness_sparsity(mesh, AlexaLaplace));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, AreaMinimizer));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, SandwichLaplace, AreaMinimizer));
    compute_disney.emplace_back(stiffness_sparsity(mesh, Disney));
    compute_AQAPoly.emplace_back(stiffness_sparsity(mesh, AQAPoly_Laplace, AreaMinimizer));
    insert_points(mesh,AreaMinimizer);
    compute_SubdivTri.emplace_back(stiffness_sparsity(mesh, quadratic_Triangle_Laplace));


    mesh.read("../data/fernando_meshes/sphere/quad_6.obj");
    compute_alexa.emplace_back(stiffness_sparsity(mesh, AlexaLaplace));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, AreaMinimizer));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, SandwichLaplace, AreaMinimizer));
    compute_disney.emplace_back(stiffness_sparsity(mesh, Disney));
    compute_AQAPoly.emplace_back(stiffness_sparsity(mesh, AQAPoly_Laplace, AreaMinimizer));
    insert_points(mesh,AreaMinimizer);
    compute_SubdivTri.emplace_back(stiffness_sparsity(mesh, quadratic_Triangle_Laplace));

    file_sparse << "Diamond,";
    for (double i : compute_diamond) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "Sandwich,";
    for (double i : compute_sandwich) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "Alexa,";
    for (double i : compute_alexa) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "Disney,";
    for (double i : compute_disney) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "subdiv TriLaplace,";
    for (double i : compute_SubdivTri) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "AQAPoly,";
    for (double i : compute_AQAPoly) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse.close();
}

void testMeshes(vector<double> &compute, vector<double> &solve, int laplace, int solver, int minpoint = 0) {
    using namespace std;


    SurfaceMesh mesh;
    vector<vector<double>> V0, N, C;
    vector<vector<int>> F;
    std::cout << "===========================================" << std::endl;
    std::cout << "Fernando Meshes " << endl;
    std::cout << "===========================================" << std::endl;
    std::cout << " ----------------------------------------" << std::endl;
    std::cout << ": Fernando Hexa sphere " << endl;
    std::cout << " ----------------------------------------" << std::endl;
    mesh.read("../data/fernando_meshes/sphere/hex_6.obj");
    takeTimingLGS(mesh, compute, solve, laplace, minpoint, solver);
    std::cout << " ----------------------------------------" << std::endl;
    std::cout << ": Fernando Quad sphere " << endl;
    std::cout << " ----------------------------------------" << std::endl;
    mesh.read("../data/fernando_meshes/sphere/quad_6.obj");
    takeTimingLGS(mesh, compute, solve, laplace, minpoint, solver);

//    testOperator(mesh, laplace, minpoint);

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << ": Fernando Concave sphere " << endl;
    std::cout << " ----------------------------------------" << std::endl;
    mesh.read("../data/fernando_meshes/sphere/concave_6.obj");
    takeTimingLGS(mesh, compute, solve, laplace, minpoint, solver);
//    testOperator(mesh, laplace, minpoint);

    std::cout << " ----------------------------------------" << std::endl;
    std::cout << ": Fernando Triangle sphere " << endl;
    std::cout << " ----------------------------------------" << std::endl;
    mesh.read("../data/fernando_meshes/sphere/triangle_6.obj");
    takeTimingLGS(mesh, compute, solve, laplace, minpoint, solver);
//    testOperator(mesh, laplace, minpoint);


}

void test_timings(int solver) {
    std::string filename_, filename2_;
    if (solver == 0) {
        filename_ = "./surface_timings_compute_eigen_llt.csv";
        filename2_ = "./surface_timings_solve_eigen_llt.csv";
    } else if (solver == 1) {
        filename_ = "./surface_timings_compute_cholmod_supernodal.csv";
        filename2_ = "./surface_timings_solve_cholmod_supernodal.csv";
    } else {
        filename_ = "./surface_timings_compute_cholmod_llt.csv";
        filename2_ = "./surface_timings_solve_cholmod_llt.csv";
    }

    std::ofstream file_compute(filename_);
    std::ofstream file_solve(filename2_);
    std::vector<double> compute_alexa, compute_diamond, compute_sandwich, compute_disney,compute_AQAPoly;
    std::vector<double> solve_alexa, solve_diamond, solve_sandwich, solve_disney, solve_AQAPoly;
    file_compute << "Laplace,Hexas,Quads,Concave,Triangles\n";
    file_solve << "Laplace,Hexas,Quads,Concave,Triangles\n";

    std::cout << "Sandwich Laplace" << std::endl;
    testMeshes(compute_sandwich, solve_sandwich, SandwichLaplace, solver, AreaMinimizer);
    file_compute << "Sandwich,";
    for (int i = 0; i < compute_sandwich.size(); i++) {
        file_compute << compute_sandwich[i] << ",";
    }
    file_compute << "\n";
    file_solve << "Sandwich,";
    for (int i = 0; i < solve_sandwich.size(); i++) {
        file_solve << solve_sandwich[i] << ",";
    }
    file_solve << "\n";

    cout << " ==================================" << std::endl;
    std::cout << "Alexa Laplace" << std::endl;
    testMeshes(compute_alexa, solve_alexa, AlexaLaplace, solver);

    file_compute << "Alexa,";
    for (int i = 0; i < compute_alexa.size(); i++) {
        file_compute << compute_alexa[i] << ",";
    }
    file_compute << "\n";
    file_solve << "Alexa,";
    for (int i = 0; i < solve_alexa.size(); i++) {
        file_solve << solve_alexa[i] << ",";
    }
    file_solve << "\n";

    cout << " ==================================" << std::endl;
    std::cout << "Diamond Laplace" << std::endl;
    testMeshes(compute_diamond, solve_diamond, Diamond, solver, AreaMinimizer);

    file_compute << "Diamond,";
    for (int i = 0; i < compute_diamond.size(); i++) {
        file_compute << compute_diamond[i] << ",";
    }
    file_compute << "\n";
    file_solve << "Diamond,";
    for (int i = 0; i < solve_diamond.size(); i++) {
        file_solve << solve_diamond[i] << ",";
    }
    file_solve << "\n";

    cout << " ==================================" << std::endl;
    std::cout << "Disney Laplace" << std::endl;
    testMeshes(compute_disney, solve_disney, Disney, solver);

    file_compute << "Disney,";
    for (int i = 0; i < compute_disney.size(); i++) {
        file_compute << compute_disney[i] << ",";
    }
    file_compute << "\n";
    file_solve << "Disney,";
    for (int i = 0; i < solve_disney.size(); i++) {
        file_solve << solve_disney[i] << ",";
    }
    file_solve << "\n";

    std::cout << "AQAPoly" << std::endl;
    testMeshes(compute_AQAPoly, solve_AQAPoly, AQAPoly_Laplace, solver, AreaMinimizer);
    file_compute << "AQAPoly,";
    for (int i = 0; i < compute_AQAPoly.size(); i++) {
        file_compute << compute_AQAPoly[i] << ",";
    }
    file_compute << "\n";
    file_solve << "AQAPoly,";
    for (int i = 0; i < solve_AQAPoly.size(); i++) {
        file_solve << solve_AQAPoly[i] << ",";
    }
    file_solve << "\n";

    cout << " ==================================" << std::endl;
}

//=========================== X AXIS COMPUTATIONS==================================================
void write_DOF_mesh(pmp::SurfaceMesh mesh, std::ofstream &file_lin, std::ofstream &file_quad) {
    double dof_lin = mesh.n_vertices();
    double dof_quad = mesh.n_vertices() + mesh.n_edges();
    file_lin << dof_lin << ",";
    file_quad << dof_quad << ",";
}

void compute_DOF_x_axises(int lvl = 8) {
    std::string DOF_linear = "./DOF_linear.csv";
    std::string DOF_quadratic = "./DOF_quadratic.csv";
    std::string DOF_linear_subdiv = "./DOF_linear_subdiv.csv";
    std::string DOF_quadratic_subdiv = "./DOF_quadratic_subdiv.csv";

    std::ofstream file_DOF_lin(DOF_linear);
    std::ofstream file_DOF_quad(DOF_quadratic);
    std::ofstream file_DOF_quad_subdiv(DOF_quadratic_subdiv);
    std::ofstream file_DOF_lin_subdiv(DOF_linear_subdiv);

    pmp::SurfaceMesh mesh;

    //========SPHERES=============================
    file_DOF_lin << "Sphere_triangle,";
    file_DOF_quad << "Sphere_triangle,";
    file_DOF_lin_subdiv << "Sphere_triangle,";
    file_DOF_quad_subdiv << "Sphere_triangle,";

    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/sphere/triangle_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nSphere_quad,";
    file_DOF_quad << "\nSphere_quad,";
    file_DOF_lin_subdiv << "\nSphere_quad,";
    file_DOF_quad_subdiv << "\nSphere_quad,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/sphere/quad_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nSphere_hex,";
    file_DOF_quad << "\nSphere_hex,";
    file_DOF_lin_subdiv << "\nSphere_hex,";
    file_DOF_quad_subdiv << "\nSphere_hex,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/sphere/hex_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nSphere_concave,";
    file_DOF_quad << "\nSphere_concave,";
    file_DOF_lin_subdiv << "\nSphere_concave,";
    file_DOF_quad_subdiv << "\nSphere_concave,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/sphere/concave_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }

    //========PLANES=============================
    file_DOF_lin << "\nPlane_triangle,";
    file_DOF_quad << "\nPlane_triangle,";
    file_DOF_lin_subdiv << "\nPlane_triangle,";
    file_DOF_quad_subdiv << "\nPlane_triangle,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/grid/triangle_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nPlane_quad,";
    file_DOF_quad << "\nPlane_quad,";
    file_DOF_lin_subdiv << "\nPlane_quad,";
    file_DOF_quad_subdiv << "\nPlane_quad,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/grid/quad_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nPlane_Voronoi,";
    file_DOF_quad << "\nPlane_Voronoi,";
    file_DOF_lin_subdiv << "\nPlane_Voronoi,";
    file_DOF_quad_subdiv << "\nPlane_Voronoi,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/grid/voronoi_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);

    }
    file_DOF_lin << "\nPlane_concave,";
    file_DOF_quad << "\nPlane_concave,";
    file_DOF_lin_subdiv << "\nPlane_concave,";
    file_DOF_quad_subdiv << "\nPlane_concave,";
    for (int i = 1; i < lvl; i++) {
        std::string meshname = "../data/fernando_meshes/grid/concave_" + std::to_string(i) + ".obj";
        mesh.read(meshname);
        // DOF as x Axis
        write_DOF_mesh(mesh, file_DOF_lin, file_DOF_quad);
        insert_points(mesh,Centroid);
        write_DOF_mesh(mesh, file_DOF_lin_subdiv, file_DOF_quad_subdiv);
    }
}

//--------------------------Timings---------------------------------------------------------------------
void
compute_x_timings(int laplace, int solver, std::vector<double>& solve, std::vector<double>& compute, string meshname, int lvl=8) {
    pmp::SurfaceMesh mesh;
    if (meshname == "Plane_quad") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/quad_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            std::cout << "Mesh lvl "<< i << " with " << mesh.n_vertices()<< " vertices" << std::endl;
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Plane_triangle") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/triangle_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Plane_concave") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/concave_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Plane_Voronoi") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/voronoi_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Sphere_triangle") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/triangle_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Sphere_quad") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/quad_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Sphere_hex") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/hex_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    } else if (meshname == "Sphere_concave") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/concave_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            if(laplace == quadratic_Triangle_Laplace){
                insert_points(mesh,AreaMinimizer);
            }
            takeTimingLGS(mesh, compute, solve, laplace, AreaMinimizer,
                          solver);
        }
    }
}

void write_x_timings(int laplace, std::vector<double> timings, std::ofstream &file) {
    if (laplace == AQAPoly_Laplace) {
        file << "AQAPoly,";
    } else if (laplace == SandwichLaplace) {
        file << "[BHKB20],";

    } else if (laplace == AlexaLaplace) {
        file << "[AW11],";

    } else if (laplace == Diamond) {
        file << "[BBA21],";

    } else if (laplace == Disney) {
        file << "[dGBD20],";

    } else if (laplace == SEC) {
        file << "[dGDMD16],";

    } else if (laplace == quadratic_Triangle_Laplace) {
        file << "subdiv TriLaplace,";
    }
    for (int i = 0; i < timings.size(); i++) {
        file << timings[i] << ",";
    }
    file << "\n";
}

void Time_x_axises(int solver, std::string meshname) {
    std::string filename_, filename2_;
    if (solver == 0) {
        filename_ = "./surface_compute_llt_" + meshname + ".csv";
        filename2_ = "./surface_solve_llt_" + meshname + ".csv";
    } else if (solver == 1) {
        filename_ = "./surface_compute_cholmod_snodal_" + meshname + ".csv";
        filename2_ = "./surface_solve_cholmod_snodal_" + meshname + ".csv";
    } else {
        filename_ = "./surface_compute_cholmod_llt_" + meshname + ".csv";
        filename2_ = "./surface_solve_cholmod_llt_" + meshname + ".csv";
    }
    std::ofstream file_compute(filename_);
    std::ofstream file_solve(filename2_);

    std::vector<double> compute,solve ;
    file_compute << meshname << ",1,2,3,4,5,6,7\n";
    file_solve << meshname << ",1,2,3,4,5,6,7\n";
    int lvl = 7;

    // One iteration without saving to avoid starting noise, TODO: see why first timings are distored!
    compute_x_timings(SandwichLaplace, solver,  solve,  compute, meshname,1);
    solve.clear();
    compute.clear();

    compute_x_timings(SandwichLaplace, solver,  solve,  compute, meshname,lvl);
    write_x_timings(SandwichLaplace, solve, file_solve);
    write_x_timings(SandwichLaplace, compute,file_compute);

    solve.clear();
    compute.clear();

    compute_x_timings(AlexaLaplace, solver,  solve,  compute, meshname,lvl);
    std::cout <<"Solve size: " <<solve.size() << std::endl;
    write_x_timings(AlexaLaplace, solve, file_solve);
    write_x_timings(AlexaLaplace, compute,file_compute);

    solve.clear();
    compute.clear();


    compute_x_timings(Disney, solver,  solve,  compute, meshname,lvl);
    write_x_timings(Disney, solve, file_solve);
    write_x_timings(Disney, compute,file_compute);

    solve.clear();
    compute.clear();

    compute_x_timings(Diamond, solver,  solve,  compute, meshname,lvl);
    write_x_timings(Diamond, solve, file_solve);
    write_x_timings(Diamond, compute,file_compute);

    solve.clear();
    compute.clear();

    compute_x_timings(AQAPoly_Laplace, solver,  solve,  compute, meshname,lvl);
    write_x_timings(AQAPoly_Laplace, solve, file_solve);
    write_x_timings(AQAPoly_Laplace, compute,file_compute);


    solve.clear();
    compute.clear();

    compute_x_timings(quadratic_Triangle_Laplace, solver,  solve,  compute, meshname,lvl);
    write_x_timings(quadratic_Triangle_Laplace, solve, file_solve);
    write_x_timings(quadratic_Triangle_Laplace, compute,file_compute);

    solve.clear();
    compute.clear();


    cout << " ==================================" << std::endl;
}

void stiffness_timings(pmp::SurfaceMesh& mesh,int laplace, int minpoint, std::vector<double>& timing)
{
    pmp::Timer t;
    Eigen::SparseMatrix<double> S;
    setup_stiffness_matrices(mesh,S,laplace,minpoint);
    t.start();
    for (int i = 0; i < 20; i++) {
        setup_stiffness_matrices(mesh,S,laplace,minpoint);
    }
    t.stop();
    std::cout << "Computation time of stiffness matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
    timing.emplace_back(t.elapsed() / 20.0);
}

void compute_stiffness_timings(int laplace, int minpoint, std::vector<double>& timing, string meshname, int lvl) {
    pmp::SurfaceMesh mesh;
    if (meshname == "Quad_plane_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/quad_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);
        }
    } else if (meshname == "Triangle_plane_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/triangle_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Concave_plane_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/concave_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Voronoi_plane_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/grid/voronoi_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Triangle_sphere_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/triangle_" + std::to_string(i) + ".obj";
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Quad_sphere_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/quad_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Hex_sphere_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/hex_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);

        }
    } else if (meshname == "Concave_sphere_") {
        for (int i = 1; i < lvl; i++) {
            std::string data_path = "../data/fernando_meshes/sphere/concave_" + std::to_string(i) + ".obj";
            mesh.read(data_path);
            stiffness_timings(mesh,laplace,minpoint,timing);
        }
    }
}


void timings_stiffness_matrix_construction(std::string meshname){
    std::string filename_;
    filename_ = "./surface_stiffness_timings_"+meshname+".csv";
    std::ofstream file_(filename_);


    std::vector<double> timing ;
    int lvl = 8;
    file_ << meshname << ",1,2,3,4,5,6,7\n";

    // One test run without writing
    compute_stiffness_timings(SandwichLaplace, AreaMinimizer,timing,meshname,lvl);
    timing.clear();


    compute_stiffness_timings(AQAPoly_Laplace, Centroid,timing,meshname,lvl);
    write_x_timings(AQAPoly_Laplace, timing, file_);
    timing.clear();
    compute_stiffness_timings(SandwichLaplace, AreaMinimizer,timing,meshname,lvl);
    write_x_timings(SandwichLaplace, timing, file_);
    timing.clear();
    compute_stiffness_timings(Diamond, AreaMinimizer,timing,meshname,lvl);
    write_x_timings(Diamond, timing, file_);
    timing.clear();
    compute_stiffness_timings(AlexaLaplace, AreaMinimizer,timing,meshname,lvl);
    write_x_timings(AlexaLaplace, timing, file_);
    timing.clear();
    compute_stiffness_timings(Disney, AreaMinimizer,timing,meshname,lvl);
    write_x_timings(Disney, timing, file_);
    timing.clear();
}
//================================================================================================

int main(int argc, const char *argv[]) {
//
//    //Eigen LLT
//    test_timings(0);
//    //Cholmod Supernodal
//    test_timings(1);
//    //Cholmod LLT
//    test_timings(2);

      compute_sparsity();
//    compute_DOF_x_axises();
//    int solver = 2;
//    Time_x_axises(solver,"Plane_quad");
//    Time_x_axises(solver,"Sphere_quad");
//
//    Time_x_axises(solver,"Plane_triangle");
//    Time_x_axises(solver,"Sphere_triangle");
//
//    Time_x_axises(solver,"Plane_concave");
//    Time_x_axises(solver,"Sphere_concave");
//
//    Time_x_axises(solver,"Plane_Voronoi");
//    Time_x_axises(solver,"Sphere_hex");

}
