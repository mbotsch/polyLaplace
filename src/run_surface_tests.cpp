//=============================================================================
#include <iostream>
#include <fstream>
#include "Surface/SpectralProcessing.h"
#include "Surface/LaplaceConstruction.h"
#include "unsupported/Eigen/SparseExtra"

//=============================================================================

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    Disney = 5,
    SEC = 6,
    AQAPoly_Laplace = 7,
    quadratic_Triangle_Laplace = 8,
    Harmonic = 10

};

enum InsertedPoint
{
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2,
};

enum Function
{
    LinearPrecision = 0,
    Evalues = 1,
    Franke2d = 2
};

enum Sphere_Tesselation
{
    Triangle = 0,
    Quad = 1,
    Concave = 2,
    Hex = 3,
};

//=============================================================================
double inverse_mean_edgelenth(SurfaceMesh &mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {

        avgLen += mesh.edge_length(e);
    }

    return 1.0 / (avgLen / mesh.n_edges());
}

void normalize(SurfaceMesh &mesh)
{
    for (auto v : mesh.vertices())
    {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

void solve_eigenvalues_per_mesh(SurfaceMesh &mesh, std::ofstream &file,
                                std::string filename)
{
//    double rmse = solve_eigenvalue_problem(mesh, Diamond, AreaMinimizer, 1,
//                                           Vertices, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AlexaLaplace, AreaMinimizer, 1,
//                                    Vertices, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AQAPoly_Laplace, AreaMinimizer, 1,
//                                    Vertices, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, Disney, AreaMinimizer, 1, Vertices,
//                                    filename);
//    file << rmse << ",";
   double rmse = solve_eigenvalue_problem(mesh, Harmonic, AreaMinimizer, 1,Vertices, filename);
    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AQAPoly_Laplace, AreaMinimizer, 2,
//                                    Edges, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AQAPoly_Laplace, AreaMinimizer, 1,
//                                    Vertices, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AQAPoly_Laplace, AreaMinimizer, 2,
//                                    Vertices, filename);
//    file << rmse << ",";
//    rmse = solve_eigenvalue_problem(mesh, AQAPoly_Laplace, AreaMinimizer, 3,
//                                    Vertices, filename);
//    file << 0.0 << ",";

//    file << rmse << ",\n";
}
void write_diamond_stiffness()
{
    /*    for (int i = 1; i < lvl; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../../NaiveProlongation/voronoi_" +
                   std::to_string(i) + "c2.obj";
        std::string filename = "S_diamond_voronoi_"+ std::to_string(i)+".mtx";
        mesh.read(meshname);
        Eigen::SparseMatrix<double> S;
        setup_stiffness_matrices(mesh,S,Diamond,AreaMinimizer);
        Eigen::saveMarket(S,filename);
    }*/
    for (int i = 6; i < 8; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../data/fernando_meshes/grid/voronoi_" + std::to_string(i) + ".obj";
        std::string filename = "S_diamond_voronoi_" + std::to_string(i) + ".mtx";
        mesh.read(meshname);
        Eigen::SparseMatrix<double> S;
        setup_stiffness_matrices(mesh, S, Diamond, AreaMinimizer);
        Eigen::saveMarket(S, filename);
    }
}

void write_evalue_test_results(int lvl =7)
{
    std::string filename_;
    std::ofstream file(filename_);

    std::string filename;
    std::ofstream error_file("evalue_error_quad.csv");
    error_file
        << "Diamond,Alexa,Sandwich, Disney,Harmonic, Aqapoly\n";
//        << "Diamond,Alexa,Sandwich, Disney, Aqapoly ,AQAPoly v1,AQAPoly "
//                  "v2,AQAPoly v3,Subdiv Cotan,\n";
    for (int i = 2; i < lvl; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../data/fernando_meshes/sphere/quad_" +
                   std::to_string(i) + ".obj";
        filename = "quad_"+ std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file, filename);
//        insert_points(mesh, AreaMinimizer);
////        normalize(mesh);
//        double rmse =
//            solve_eigenvalue_problem(mesh, quadratic_Triangle_Laplace,
//                                     AreaMinimizer, 2, Edges, filename);
        error_file << "\n";
    }
    error_file.close();

    std::ofstream error_file_hex("evalue_error_hex.csv");
    error_file_hex
        << "Diamond,Alexa,Sandwich, Disney, Aqapoly ,Subdiv Cotan,\n";
//        << "Diamond,Alexa,Sandwich, Disney, Aqapoly,AQAPoly "
//                      "v1,AQAPoly v2,AQAPoly v3,Subdiv Cotan,\n";
    for (int i = 2; i < lvl; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../data/fernando_meshes/sphere/hex_" + std::to_string(i) +
                   ".obj";
        filename = "hex_"+ std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_hex, filename);
//        insert_points(mesh, AreaMinimizer);
////        normalize(mesh);
//
//        double rmse =
//            solve_eigenvalue_problem(mesh, quadratic_Triangle_Laplace,
//                                     AreaMinimizer, 2, Edges, filename);
        error_file_hex <<"\n";
    }

    error_file_hex.close();
    std::ofstream error_file_triangle("evalue_error_triangle.csv");
    error_file_triangle
        << "Diamond,Alexa,Sandwich, Disney, Aqapoly ,Subdiv Cotan,\n";
//        << "Diamond,Alexa,Sandwich, Disney, Aqapoly,AQAPoly "
//                           "v1,AQAPoly v2,AQAPoly v3,Subdiv Cotan,\n";
    for (int i = 2; i < lvl; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        int res = i;
        meshname = "../data/fernando_meshes/sphere/triangle_" +
                   std::to_string(res) + ".obj";
        filename = "triangle_"+ std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_triangle, filename);
//        insert_points(mesh, AreaMinimizer);
////        normalize(mesh);
//        double rmse =
//            solve_eigenvalue_problem(mesh, quadratic_Triangle_Laplace,
//                                     AreaMinimizer, 2, Edges, filename);
        error_file_triangle << "\n";
    }
    error_file_triangle.close();

    std::ofstream error_file_concave("evalue_error_concave.csv");
    error_file_concave
        << "Diamond,Alexa,Sandwich, Disney, Aqapoly ,Subdiv Cotan,\n";
//        << "Diamond,Alexa,Sandwich, Disney, Aqapoly,AQAPoly "
//                          "v1,AQAPoly v2,AQAPoly v3,Subdiv Cotan,\n";
    for (int i = 2; i < lvl; i++)
    {
        SurfaceMesh mesh;
        std::string meshname;
        int res = i;
        meshname = "../data/fernando_meshes/sphere/concave_" +
                   std::to_string(res) + ".obj";
        filename = "concave_"+ std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_concave, filename);
//        insert_points(mesh, AreaMinimizer);
////        normalize(mesh);
//        double rmse =
//            solve_eigenvalue_problem(mesh, quadratic_Triangle_Laplace,
//                                     AreaMinimizer, 2, Edges, filename);
        error_file_concave << "\n";
    }
    error_file_concave.close();
}

//----------------------------------------------------------------------------

int main()
{
    write_diamond_stiffness();

//    write_evalue_test_results(5);
}
//=============================================================================
