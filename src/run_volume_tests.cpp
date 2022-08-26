//=============================================================================
#include <iostream>
#include "Volume/Eigenmodes.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "VolumeMesh/VolumeSubdivision.h"
#include "Volume/HarmonicBasis.h"
#include "Volume/Diamond_3D.h"
#include "unsupported/Eigen/SparseExtra"
//=============================================================================

enum LaplaceMethods
{
    Diamond = 0,
    Dual_Laplace = 2,
    Harmonic = 3,
    Sandwich = 4,
    AQAPoly_MG = 5,
    AQAPoly = 6,
    refinedTetLaplace = 7
};
enum VolumePoints
{
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints
{
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

enum Function
{
    Franke3d = 1,
    Multigrid = 2,
    Coarse_Dim = 3
};

//=============================================================================
double inverse_mean_edgelenth(VolumeMesh &mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {

        avgLen += mesh.length(e);
    }

    return 1.0 / (avgLen / mesh.n_edges());
}

void write_eigenmodes_test_results()
{
    VolumeMesh mesh;
    //    std::string filename_ = "./errors_eigenmodes_3D.csv";
    //
    //    std::ofstream file(filename_);

    //    Eigen::VectorXd evalues;
    //    std::string filename;
    //    file << "Hexa [BBA21],Hexa [BHBK20],Hexa Harmonic,Hexa AQAPoly,Pyramids "
    //            "[BBA21],Pyramids [BHBK20],Pyramids Harmonic,Pyramids AQAPoly,Truncated "
    //            "[BBA21],Truncated [BHBK20],Truncated Harmonic,Truncated AQAPoly,"
    //            "Tetrahedral_Coarse [BBA21],Tetrahedral_Coarse [BHBK20],Tetrahedral_Coarse Harmonic,"
    //            "Tetrahedral_Coarse AQAPoly,"
    //            "Tetrahedral [BBA21],Tetrahedral [BHBK20],Tetrahedral Harmonic,Tetrahedral AQAPoly,"
    //         << std::endl;

    {
        std::string meshname = "../data/kong.ovm";
        mesh.read(meshname);
        Eigen::SparseMatrix<double> S;
        setup_3D_stiffness_matrix(mesh, S, Diamond, Quadratic_Areas_,
                                  Quadratic_Volume_);
        Eigen::saveMarket(S, "S_diamond_kong.mtx");
    }

    //    {
    //        std::string meshname =
    //            "../data/Volume/unit_balls/Hexahedral_sphere.ovm";
    //        mesh.read(meshname);
    //        filename = "Hexahedra";
    //        double rmse;
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Diamond,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    ////        file << rmse << ",";
    ////        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    ////                                        Quadratic_Areas_, Quadratic_Volume_,
    ////                                        filename, 1, Vertices);
    ////        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Harmonic,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    //        file << rmse << ",";
    ////        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    ////                                        Quadratic_Areas_, Quadratic_Volume_,
    ////                                        filename,2, Edges);
    ////        file << rmse << ",";
    //    }
    //    {
    //        std::string meshname = "../data/Volume/unit_balls/Pyramid_sphere.ovm";
    //        mesh.read(meshname);
    //        filename = "Pyramids";
    //        double rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Diamond,
    //                                               Quadratic_Areas_,
    //                                               Quadratic_Volume_, filename);
    ////        file << rmse << ",";
    ////        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    ////                                        Quadratic_Areas_, Quadratic_Volume_,
    ////                                        filename, 1, Vertices);
    ////        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Harmonic,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    //        file << rmse << ",";
    ////        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    ////                                        Face_Centroid, Cell_Centroid_, filename,
    ////                                        2, Edges);
    ////        file << rmse << ",";
    ////    }
    //    {
    //        std::string meshname = "../data/Volume/unit_balls/Truncated_sphere.ovm";
    //        mesh.read(meshname);
    //        filename = "Truncated";
    //        double rmse;
    //        double rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Diamond,
    //                                               Quadratic_Areas_,
    //                                               Quadratic_Volume_, filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename, 1, Vertices);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Harmonic,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Face_Centroid, Cell_Centroid_, filename,
    //                                        2, Edges);
    //        file << rmse << ",";
    //    }
    //    {
    //        std::string meshname =
    //            "../data/Volume/unit_balls/tetrahedra_sphere_irreg.ovm";
    //        mesh.read(meshname);
    //        filename = "Tetrahedra Coarse";
    //        double rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Diamond,
    //                                               Quadratic_Areas_,
    //                                               Quadratic_Volume_, filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename, 1, Vertices);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Harmonic,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Face_Centroid, Cell_Centroid_, filename,
    //                                        2, Edges);
    //        file << rmse << ",";
    //    }
    //    {
    //        std::string meshname =
    //            "../data/Volume/unit_balls/tetrahedra_sphere_reg.ovm";
    //        mesh.read(meshname);
    //        filename = "Tetrahedra";
    //        double rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Diamond,
    //                                               Quadratic_Areas_,
    //                                               Quadratic_Volume_, filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename, 1, Vertices);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, Harmonic,
    //                                        Quadratic_Areas_, Quadratic_Volume_,
    //                                        filename);
    //        file << rmse << ",";
    //        rmse = solve_eigenvalue_problem(mesh, meshname, evalues, AQAPoly,
    //                                        Face_Centroid, Cell_Centroid_, filename,
    //                                        2, Edges);
    //        file << rmse << ",";
    //    }
    //    file.close();
}

void write_data(VolumeMesh &mesh, std::string meshname, int laplace,
                std::ofstream &error_file, std::ofstream &timings_file,
                int function, int vcycle = 10)
{

    if (function == Franke3d)
    {
        if (laplace == Diamond)
        {
            error_file << solve_franke_poisson(mesh, Diamond, Quadratic_Areas_,
                                               Quadratic_Volume_, 1);
        }
        else if (laplace == Sandwich)
        {
            error_file << solve_franke_poisson(mesh, Sandwich, Quadratic_Areas_,
                                               Quadratic_Volume_, 1);
        }
        else if (laplace == Harmonic)
        {
            error_file << solve_franke_poisson(mesh, Harmonic, Quadratic_Areas_,
                                               Quadratic_Volume_, 1, meshname);
        }
        error_file << ",";
    }
    else
    {
        std::cout << " Function not implemented" << std::endl;
    }
}


void write_all_laplace_data(VolumeMesh &mesh, std::string meshname,
                            std::ofstream &file, std::ofstream &timings_file)
{
    write_data(mesh, meshname, Diamond, file, timings_file, Franke3d);
    write_data(mesh, meshname, Sandwich, file, timings_file, Franke3d);
    //    write_data(mesh, meshname, Harmonic, file, timings_file, Franke3d);
}

void write_text_headers(Function function, std::ofstream &file_error,
                        std::ofstream &file_timings)
{
    if (function == Franke3d)
    {
        //        file_error << "[BBA21],[BHKB20],Harmonic,AQAPoly,Tetrahedra, resolution, resolution tetrahedra"
        //                        << std::endl;
        //            file_error << "[BBA21],[BHKB20],AQAPoly,Tetrahedra, resolution, resolution tetrahedra"
        //           << std::endl;
        file_error << "[BBA21],[BHKB20],AQAPoly,Tetrahedra, dof liner,dof "
                      "quadratic, dof tetrahedra"
                   << std::endl;
    }
    else if (function == Multigrid)
    {
        file_error << "AQAPoly direct,Vcycles 3 it 3,Vcycles 3 it 10,Vcycles 3 "
                      "it 20,Vcycles 10 it 3,Vcycles 10 it 10,Vcycles 10 it "
                      "30,Vcycles 20 it 30,resolution"
                   << std::endl;
        file_timings
            << "Solver,Got system matrices,  Cholmod (Supernodal LLT) "
               "factorization,Solved direct system,Constructed "
               "multigrid system,Solved multigrid system,RMS,RMS Direct,DOF"
            << std::endl;
    }
    else if (function == Coarse_Dim)
    {
        file_error
            << "coarse Dim 0,coarse Dim 1,coarse Dim 2,coarse Dim 3,resolution"
            << std::endl;
        file_timings << "Solver,Got system matrices,  Cholmod (Supernodal LLT) "
                        "factorization,Solved direct system,Constructed "
                        "multigrid system,Solved multigrid system,RMS,DOF"
                     << std::endl;
    }
}

double write_3D_convergence_data_csv(Function function, int lvl_end = 6,
                                     int lvl_start = 2)
{

    VolumeMesh mesh;
    std::string filename_, timings_filename_;

    // --------------hex cubes----------------------------------

    if (function == Franke3d)
    {
        filename_ = "./errors_poisson_Franke3D_hex.csv";
        timings_filename_ = "./timings_mg_Franke3D_hex.csv";
    }
    else if (function == Multigrid)
    {
        filename_ = "./errors_mg_Franke3D_hex.csv";
        timings_filename_ = "./timings_mg_Franke3D_hex.csv";
    }
    else if (function == Coarse_Dim)
    {
        filename_ = "./errors_coarsedim_Franke3D_hex.csv";
        timings_filename_ = "./timings_mg_coarsedim_Franke3D_hex.csv";
    }

    std::ofstream file_franke_hex(filename_);
    std::ofstream file_timings_hex(timings_filename_);

    write_text_headers(function, file_franke_hex, file_timings_hex);

    for (int i = lvl_start; i < lvl_end; i++)
    {
        std::string meshname = "../data/Volume_data/cubes/cube_hexahedron_" +
                               std::to_string(i) + ".ovm";
        mesh.read(meshname);
        std::cout << meshname << std::endl;
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d)
        {
            write_all_laplace_data(mesh, meshname, file_franke_hex,
                                   file_timings_hex);
        }

                file_franke_hex << res << "," ;
                VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
                res = inverse_mean_edgelenth(mesh);
                file_franke_hex << res << std::endl;
//        -----------------Dof---------------------------
//        file_franke_hex << mesh.n_vertices() << ",";
//        file_franke_hex << mesh.n_vertices() + mesh.n_edges() << ",";
//        VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
//        file_franke_hex << mesh.n_vertices() + mesh.n_edges() << std::endl;
    }
    file_franke_hex.close();
    file_timings_hex.close();
    //--------------pyramid cubes----------------------------------

    if (function == Franke3d)
    {
        filename_ = "./errors_poisson_Franke3D_pyramid.csv";
        timings_filename_ = "./timings_mg_Franke3D_pyramid.csv";
    }
    else if (function == Multigrid)
    {
        filename_ = "./errors_mg_Franke3D_pyramid.csv";
        timings_filename_ = "./timings_mg_Franke3D_pyramid.csv";
    }
    else if (function == Coarse_Dim)
    {
        filename_ = "./errors_coarsedim_Franke3D_pyramid.csv";
        timings_filename_ = "./timings_mg_coarsedim_Franke3D_pyramid.csv";
    }
    std::ofstream file_franke_pyramid(filename_);
    std::ofstream file_timings_pyramid(timings_filename_);

    write_text_headers(function, file_franke_pyramid, file_timings_pyramid);

    for (int i = lvl_start; i < lvl_end; i++)
    {
        std::string meshname = "../data/Volume_data/cubes/cube_pyramids_" +
                               std::to_string(i) + ".ovm";
        std::cout << meshname << std::endl;

        mesh.read(meshname);
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d)
        {
            write_all_laplace_data(mesh, meshname, file_franke_pyramid,
                                   file_timings_pyramid);
        }
                file_franke_pyramid << res << ",";
                VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
                res = inverse_mean_edgelenth(mesh);
                file_franke_pyramid << res << std::endl;
//        -----------------Dof---------------------------
//        file_franke_pyramid << mesh.n_vertices() << ",";
//        file_franke_pyramid << mesh.n_vertices() + mesh.n_edges() << ",";
//        VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
//        file_franke_pyramid << mesh.n_vertices() + mesh.n_edges() << std::endl;
    }
    file_franke_pyramid.close();
    file_timings_pyramid.close();
    //--------------truncated cubes----------------------------------
    if (function == Franke3d)
    {
        filename_ = "./errors_poisson_Franke3D_truncated.csv";
        timings_filename_ = "./timings_mg_Franke3D_truncated.csv";
    }
    else if (function == Multigrid)
    {
        filename_ = "./errors_mg_Franke3D_truncated.csv";
        timings_filename_ = "./timings_mg_Franke3D_truncated.csv";
    }
    else if (function == Coarse_Dim)
    {
        filename_ = "./errors_coarsedim_Franke3D_truncated.csv";
        timings_filename_ = "./timings_mg_coarsedim_Franke3D_truncated.csv";
    }
    std::ofstream file_franke_truncated(filename_);
    std::ofstream file_timings_truncated(timings_filename_);

    write_text_headers(function, file_franke_truncated, file_timings_truncated);

    for (int i = lvl_start; i < lvl_end; i++)
    {
        std::string meshname = "../data/Volume_data/cubes/cube_truncated_" +
                               std::to_string(i) + ".ovm";
        std::cout << meshname << std::endl;

        mesh.read(meshname);
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d)
        {
            write_all_laplace_data(mesh, meshname, file_franke_truncated,
                                   file_timings_truncated);
        }

                file_franke_truncated << res << ",";
                VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
                res = inverse_mean_edgelenth(mesh);
                file_franke_truncated << res << std::endl;
        //-----------------Dof---------------------------
//        file_franke_truncated << mesh.n_vertices() << ",";
//        file_franke_truncated << mesh.n_vertices() + mesh.n_edges() << ",";
//        VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
//        file_franke_truncated << mesh.n_vertices() + mesh.n_edges()
//                              << std::endl;
    }
    file_franke_truncated.close();
    file_timings_truncated.close();
    //--------------tetrahedron cubes----------------------------------

    //    if (function == Franke3d)
    //    {
    //        filename_ = "./errors_poisson_Franke3D_tetrahedra.csv";
    //        timings_filename_ = "./timings_mg_coarsedim_Franke3D_tetrahedra.csv";
    //    }
    //    else if (function == Multigrid)
    //    {
    //        filename_ = "./errors_mg_Franke3D_tetrahedra.csv";
    //        timings_filename_ = "./timings_mg_Franke3D_tetrahedra.csv";
    //    }
    //    else if (function == Coarse_Dim)
    //    {
    //        filename_ = "./errors_coarsedim_Franke3D_tetrahedra.csv";
    //        timings_filename_ = "./timings_mg_Franke3D_tetrahedra.csv";
    //    }
    //    std::ofstream file_franke_tet(filename_);
    //    std::ofstream file_timings_tet(timings_filename_);
    //
    //    write_text_headers(function,file_franke_tet,file_timings_tet);
    //
    //    for (int i = lvl_start-1; i < lvl_end - 1; i++)
    //    {
    //        std::string meshname =
    //            "../data/Volume_data/cubes/cube_tet_" + std::to_string(i) + ".ovm";
    //        std::cout << meshname << std::endl;
    //
    //        mesh.read(meshname);
    //        double res = inverse_mean_edgelenth(mesh);
    //        if (function == Franke3d)
    //        {
    //            write_all_laplace_data(mesh, meshname, file_franke_tet,
    //                                   file_timings_tet);
    //        }
    //        else if (function == Multigrid)
    //        {
    //            write_all_mg_laplace_data(meshname, file_franke_tet,
    //                                      file_timings_tet);
    //        }
    //        else if (function == Coarse_Dim)
    //        {
    //            write_coarseDim_laplace_data(mesh, meshname, file_franke_tet,
    //                                         file_timings_tet);
    //        }
    //        file_franke_tet << res << std::endl;
    //    }
    //    file_franke_tet.close();
    //    file_timings_tet.close();

    //    //--------------voronoi cell cubes----------------------------------
    if (function == Franke3d)
    {
        filename_ = "./errors_poisson_Franke3D_voronoi.csv";
        timings_filename_ = "./timings_mg_Franke3D_voronoi.csv";
    }
    else if (function == Multigrid)
    {
        filename_ = "./errors_mg_Franke3D_voronoi.csv";
        timings_filename_ = "./timings_mg_Franke3D_voronoi.csv";
    }
    else if (function == Coarse_Dim)
    {
        filename_ = "./errors_coarsedim_Franke3D_voronoi.csv";
        timings_filename_ = "./timings_mg_coarsedim_Franke3D_voronoi.csv";
    }
    std::ofstream file_franke_voronoi(filename_);
    std::ofstream file_timings_voronoi(timings_filename_);

    if (function == Multigrid)
    {
        file_franke_voronoi
            << "AQAPoly direct,Vcycles 1,Vcycles 2,Vcycles 3,Vcycles 5,Vcycles "
               "10,Vcycles 30,Vcycles 60,resolution"
            << std::endl;
        file_timings_voronoi
            << "Solver,Got system matrices,  Cholmod (Supernodal LLT) "
               "factorization,Solved direct system,Constructed multigrid "
               "system,Solved multigrid system,RMS,DOF"
            << std::endl;
    }
    else
    {
        write_text_headers(function, file_franke_voronoi, file_timings_voronoi);
    }

    for (int i = lvl_start - 1; i < lvl_end - 1; i++)
    {
        std::string meshname = "../data/Volume_data/cubes/cube_voronoi_" +
                               std::to_string(i) + ".ovm";
        mesh.read(meshname);
        std::cout << meshname << std::endl;
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d)
        {
            write_all_laplace_data(mesh, meshname, file_franke_voronoi,
                                   file_timings_voronoi);
        }

                file_franke_voronoi << res << ",";
                VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
                res = inverse_mean_edgelenth(mesh);
                file_franke_voronoi << res << std::endl;
//        -----------------Dof---------------------------
//        file_franke_voronoi << mesh.n_vertices() << ",";
//        file_franke_voronoi << mesh.n_vertices() + mesh.n_edges() << ",";
//        VolumeSubdivision(mesh).tetrahedra(Quadratic_Areas_, Quadratic_Volume_);
//        file_franke_voronoi << mesh.n_vertices() + mesh.n_edges() << std::endl;
    }
    file_franke_voronoi.close();
    file_timings_voronoi.close();
}


void write_DOF()
{
    VolumeMesh mesh;
    std::string filename_cube = "./cube_x_DOF.csv";

    std::ofstream file(filename_cube);
    int lvl = 7;
    Eigen::VectorXd results;

    file << "Lvl,hex lin, hex quas,pyramid lin, pyramid quad,truncated lin, "
            "truncated quad, voronoi lin, voronoi quad, tetrahedra lin, "
            "tetrahedra quad,"
         << std::endl;
    for (int i = 2; i < lvl; i++)
    {

        file << i << ",";
        std::string meshname = "../data/Volume_data/cubes/cube_hexahedron_" +
                               std::to_string(i) + ".ovm";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file << Dof << ",";

        meshname = "../data/Volume_data/cubes/cube_pyramids_" +
                   std::to_string(i) + ".ovm";
        mesh.read(meshname);
        Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file << Dof << ",";

        meshname = "../data/Volume_data/cubes/cube_truncated_" +
                   std::to_string(i) + ".ovm";
        mesh.read(meshname);

        Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        meshname = "../data/Volume_data/cubes/cube_voronoi_" +
                   std::to_string(i - 1) + ".ovm";
        mesh.read(meshname);
        Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file << Dof << ",";

        meshname = "../data/Volume_data/cubes/cube_tet_" +
                   std::to_string(i - 1) + ".ovm";
        mesh.read(meshname);

        Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file << Dof << ",\n";
    }
    file.close();

    std::string filename_ball = "./ball_x_DOF.csv";

    std::ofstream file_2(filename_ball);
    file_2 << "Hex lin, Hex quas,Pyramid lin, Pyramid quad,Truncated lin, "
              "Truncated quad, Tetrahedral_Coarse lin, Tetrahedral_Coarse "
              "quad, Tetrahedral lin, Tetrahedral quad,"
           << std::endl;
    {
        std::string meshname =
            "../data/Volume_data/spheres/Hexahedral_sphere.ovm";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
    }
    {
        std::string meshname = "../data/Volume_data/spheres/Pyramid_sphere.ovm";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
    }
    {
        std::string meshname =
            "../data/Volume_data/spheres/Truncated_sphere.ovm";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
    }
    {
        std::string meshname =
            "../data/Volume_data/spheres/tetrahedra_sphere_2.mesh";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
    }
    {
        std::string meshname = "../data/Tetgen/spheres/sphere2.mesh";
        mesh.read(meshname);
        int Dof = 0;
        for (auto v : mesh.vertices())
        {
            if (!mesh.is_boundary(v))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
        for (auto e : mesh.edges())
        {
            if (!mesh.is_boundary(e))
            {
                Dof++;
            }
        }
        file_2 << Dof << ",";
    }
    file_2.close();
}

void write_harmonic_timings()
{
    std::string timings_filename2_;

    // --------------hex cube----------------------------------

    timings_filename2_ = "./timings_harmonic_selected.csv";

    std::ofstream file_timings(timings_filename2_);

    file_timings << "Voronoi Cube 2 " << std::endl;
    std::string meshname = "../data/Volume_data/cubes/cube_voronoi_2.ovm";
    {
        PolyhedralMesh mesh(meshname);
        Eigen::SparseMatrix<double> S;
        setup_3D_harmonic_stiffness_matrix(mesh, S, file_timings);
    }
    file_timings << "Hex Cube 2 " << std::endl;
    meshname = "../data/Volume_data/cubes/cube_hexahedron_2.ovm";
    {
        PolyhedralMesh mesh(meshname);
        Eigen::SparseMatrix<double> S;
        setup_3D_harmonic_stiffness_matrix(mesh, S, file_timings);
    }
    file_timings << "Voronoi Cube 3 " << std::endl;
    meshname = "../data/Volume_data/cubes/cube_voronoi_3.ovm";
    {
        PolyhedralMesh mesh(meshname);
        Eigen::SparseMatrix<double> S;
        setup_3D_harmonic_stiffness_matrix(mesh, S, file_timings);
    }
    file_timings << "Hex Cube 3" << std::endl;
    meshname = "../data/Volume_data/cubes/cube_hexahedron_3.ovm";
    {
        PolyhedralMesh mesh(meshname);
        Eigen::SparseMatrix<double> S;
        setup_3D_harmonic_stiffness_matrix(mesh, S, file_timings);
    }
    //    file_timings << "Dragon" << std::endl;
    //    meshname = "../data/Dragon_hex.ovm";
    //    {
    //        PolyhedralMesh mesh(meshname);
    //        Eigen::SparseMatrix<double> S;
    //        setup_3D_harmonic_stiffness_matrix(mesh, S, file_timings);
    //    }
}
//----------------------------------------------------------------------------

int main()
{
    //    write_eigenmodes_test_results();
    //    write_DOF();
    //          write_3D_mg_convergence_data_csv(Franke3d,6);
    write_3D_convergence_data_csv(Franke3d, 6, 2);
    //        write_3D_convergence_data_csv(Coarse_Dim,5);
    //        write_3D_convergence_data_csv(Multigrid,6,5);
    //    write_3D_mg();
    //    write_harmonic_timings();
}
//=============================================================================
