//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <fstream>
#include "Surface/PolyLaplace.h"
#include "Surface/Curvature.h"
#include "Surface/DisneyLaplace.h"
#include "Surface/Poisson_System.h"
#include "Surface/SpectralProcessing.h"
#include "Surface/HarmonicBasis2D.h"
#include "Surface/SmoothSubdivBasis.h"

//=============================================================================

using namespace pmp;

//=============================================================================

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    IntrinsicDelaunay = 4,
    Disney = 5,
    SEC = 6,
    AQAPoly = 7,
    quadratic_Triangle_Laplace = 8,
    AQAPoly_MG = 9,
    Harmonic = 10,
    Smoothsubdiv =11

};
enum InsertedPoint
{
    Centroid = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2
};

enum Function
{
    poisson_SH = 0,
    Franke2d = 2,
    curvature = 1,
    SH = 3,
    cond_number = 4
};

void normalize(SurfaceMesh &mesh)
{
    for (auto v : mesh.vertices())
    {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

void write_data(SurfaceMesh &mesh, int laplace, CoarseDimension coarseningType,
                int degree, std::ofstream &file, int function, int vcycles = 10,
                int l = 2, int m = 0)
{

    if (function == curvature)
    {
        Curvature analyzer(mesh, true);
        bool lumped = false;
        if (laplace == Diamond)
        {
            lumped = true;
        }
        file << analyzer.compute_curvature_error(laplace, AreaMinimizer, lumped,
                                                 degree, coarseningType);
        file << ",";
    }
    else if (function == SH)
    {
        bool lumped = true;
        if (laplace == Diamond)
        {
            lumped = false;
        }
        if (laplace == quadratic_Triangle_Laplace)
        {
            insert_points(mesh, AreaMinimizer);
            file << rmse_sh(mesh, laplace, AreaMinimizer, lumped, degree,
                            coarseningType);
        }
        else
        {
            file << rmse_sh(mesh, laplace, AreaMinimizer, lumped, degree,
                            coarseningType);
        }

        file << ",";
    }
    else if (function == poisson_SH || function == Franke2d)
    {
        if (laplace == AQAPoly_MG)
        {
            bool direct = true;
            file << solve_AQAPoly_Poisson_mg(
                mesh, coarseningType, degree, !direct,
                RELAXER_PARALLEL_GAUSS_SEIDEL, vcycles);
        }
        else if (laplace == Harmonic)
        {
            file << solve_2D_Franke_harmonic(mesh);
        }
        else if(laplace == Smoothsubdiv)
        {
            file << solve_poisson_non_dirichlet(mesh,1, 1);

        }
        else
        {
            if (laplace == quadratic_Triangle_Laplace)
            {
                if (function == poisson_SH)
                {
                    insert_points(mesh, AreaMinimizer);
                    double error = solve_poisson_system(
                        mesh, laplace, AreaMinimizer, function, degree,
                        coarseningType, l, m);
                    file << error;
                }
                else
                {
                    double error =
                        solve_poisson_system(mesh, AQAPoly, AreaMinimizer,
                                             function, 2, Refined_mesh, l, m);
                    file << error;
                }
            }
            else
            {
                double error =
                    solve_poisson_system(mesh, laplace, AreaMinimizer, function,
                                         degree, coarseningType, l, m);
                file << error;
            }
        }

        file << ",";
    }
    else if (function == cond_number)
    {
        if (laplace == AQAPoly)
        {
            std::vector<double> condnrs;
            double cond =
                AQAPoly_condition_nr(mesh, coarseningType, degree, condnrs);
            file << condnrs[0] << "," << condnrs[1];
        }

        else
        {
            file << condition_number(mesh, laplace, AreaMinimizer, degree,
                                     coarseningType);
        }
        file << ",";
    }
}


double inverse_mean_edgelenth(SurfaceMesh &mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {
        avgLen += mesh.edge_length(e);
    }

    return 1.0 / (avgLen / mesh.n_edges());
}
void write_mg_data(SurfaceMesh &mesh, std::ofstream &file, int function,
                   int l = 2, int m = 2)
{
    write_data(mesh, AQAPoly, Edges, 2, file, function, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 1, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 2, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 3, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 4, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 5, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 10, l, m);
    write_data(mesh, AQAPoly_MG, Edges, 2, file, function, 15, l, m);
}

void write_all_laplace_data(SurfaceMesh &mesh, std::ofstream &file,
                            int function, int l = 2, int m = 2)
{
    poly_laplace_lambda_ = 2.0;
    write_data(mesh, AlexaLaplace, Vertices, 1, file, function, 10, l, m);
    disney_laplace_lambda_ = 1.0;
    write_data(mesh, Disney, Vertices, 1, file, function, 10, l, m);
    write_data(mesh, SandwichLaplace, Vertices, 1, file, function, 10, l, m);
    write_data(mesh, Diamond, Vertices, 1, file, function, 10, l, m);
    subdiv_lvl = 1;
    write_data(mesh, Smoothsubdiv, Vertices, 1, file, function, 10, l, m);
    subdiv_lvl = 4;
    write_data(mesh, Smoothsubdiv, Vertices, 1, file, function, 10, l, m);
//    subdiv_lvl = 6;
//    write_data(mesh, Smoothsubdiv, Vertices, 1, file, function, 10, l, m);

//    write_data(mesh, AQAPoly, Edges, 2, file, function, 10, l, m);
//    write_data(mesh, AQAPoly, Refined_mesh, 2, file, function, 10, l, m);

    //    write_data(mesh, SEC, Vertices, 1, file, function, 10, l, m);
}
void write_cond_numbers(SurfaceMesh &mesh, std::ofstream &file, int function,
                        int l = 2, int m = 2)
{
    write_data(mesh, AQAPoly, Vertices, 1, file, function, 10, l, m);
    write_data(mesh, AQAPoly, Edges, 2, file, function, 10, l, m);
    write_data(mesh, AQAPoly, Edges, 3, file, function, 10, l, m);
    write_data(mesh, AQAPoly, Edges, 4, file, function, 10, l, m);
}

void write_text_headers(std::ofstream &file_error, int function = Franke2d)
{

    //    file_error << "[AW11],[dGBD20],[BHKB20],[BBA21],Butterfly,AQAPoly,resolution" << std::endl;
    if (function == Franke2d)
    {
        file_error
//            << "[AW11],[dGBD20],[BHKB20],[BBA21],AQAPoly,Refined Triangles,Dof "
//               "linear, Dof quadratic, Dof quadratic subdiv"
//            << std::endl;
            << "[AW11],[dGBD20],[BHKB20],[BBA21],Smoothsubdiv lvl 1,Smoothsubdiv lvl 4,Dof "
               "linear, Dof quadratic, Dof quadratic subdiv"
            << std::endl;
    }
    else if (function == cond_number)
    {
        file_error << "Vertex deg 1, relative, Edge deg 2, relative, Edge deg "
                      "3, relative,Edge deg 4, relative, linear Dof, quad Dof"
                   << std::endl;
    }
}

double write_convergence_data_csv(Function function, int lvl = 7, int l = 2,
                                  int m = 0, int start_lvl = 1)
{

    SurfaceMesh mesh;
    std::string filename_;
    std::string sh = "_Y_" + std::to_string(l) + std::to_string(m);

    if ((function != Franke2d) && (function != cond_number))
    {
        //--------------quad spheres----------------------------------
        if (function == poisson_SH)
        {
            filename_ = "./errors_poisson_SH" + sh + "_quad.csv";
        }
        else if (function == curvature)
        {
            filename_ = "./errors_curvature_quad.csv";
        }
        else if (function == SH)
        {
            filename_ = "./errors_SH_band_recreation_quad.csv";
        }
        std::ofstream file(filename_);
        write_text_headers(file);
        //        file << "[AW11],[dGBD20],[BHKB20],[BBA21],AQAPoly,SEC 2.0,Subdiv "
        //                "Cotan,resolution, resolution subdiv"
        //             << std::endl;
        //            << "[AW11],[dGBD20],[BHKB20],[BBA21],AQAPoly,AQAPoly v1,AQAPoly "
        //               "v2,AQAPoly v3,Subdiv Cotan, resolution, resolution subdiv"
        //            << std::endl;

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/sphere/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file, function, l, m);
            file << res << ",";
            res = inverse_mean_edgelenth(mesh);
            file << res << std::endl;
        }
        file.close();

        //--------------concave spheres----------------------------------
        if (function == poisson_SH)
        {
            filename_ = "./errors_poisson_SH" + sh + "_concave.csv";
        }
        else if (function == curvature)
        {
            filename_ = "./errors_curvature_concave.csv";
        }
        else if (function == SH)
        {
            filename_ = "./errors_SH_band_recreation_concave.csv";
        }
        std::ofstream file_concave(filename_);
        write_text_headers(file_concave);

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/sphere/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_concave, function, l, m);
            file_concave << res << ",";
            res = inverse_mean_edgelenth(mesh);
            file_concave << res << std::endl;
        }
        file_concave.close();

        //--------------hex spheres----------------------------------
        if (function == poisson_SH)
        {
            filename_ = "./errors_poisson_SH" + sh + "_hex.csv";
        }
        else if (function == curvature)
        {
            filename_ = "./errors_curvature_hex.csv";
        }
        else if (function == SH)
        {
            filename_ = "./errors_SH_band_recreation_hex.csv";
        }
        std::ofstream file_hex(filename_);
        write_text_headers(file_hex);

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/sphere/hex_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_hex, function, l, m);
            file_hex << res << ",";
            res = inverse_mean_edgelenth(mesh);
            file_hex << res << std::endl;
        }
        file_hex.close();

        //--------------triangle spheres----------------------------------
        if (function == poisson_SH)
        {
            filename_ = "./errors_poisson_SH" + sh + "_triangle.csv";
        }
        else if (function == curvature)
        {
            filename_ = "./errors_curvature_triangle.csv";
        }
        else if (function == SH)
        {
            filename_ = "./errors_SH_band_recreation_triangle.csv";
        }
        std::ofstream file_tri(filename_);
        write_text_headers(file_tri);

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/sphere/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_tri, function, l, m);
            file_tri << res << ",";
            res = inverse_mean_edgelenth(mesh);
            file_tri << res << std::endl;
        }
        file_tri.close();
    }
    else
    {
        //--------------quad planes----------------------------------
        if (function == Franke2d)
        {
            filename_ = "./errors_poisson_Franke2D_quad.csv";
        }
        else
        {
            filename_ = "./condition_number_plane_quad.csv";
        }
        std::ofstream file_franke_quad(filename_);
        write_text_headers(file_franke_quad, function);

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/grid/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            if (function == Franke2d)
            {
                write_all_laplace_data(mesh, file_franke_quad, function);
            }
            else
            {
                write_cond_numbers(mesh, file_franke_quad, function);
            }
            //            file_franke_quad << res << ",";
            //            insert_points(mesh, AreaMinimizer);
            //            res = inverse_mean_edgelenth(mesh);
            //            file_franke_quad << res << std::endl;

            //--------------------Dof--------------
            file_franke_quad << mesh.n_vertices() << ",";
            //            file_franke_quad << mesh.n_vertices()+mesh.n_edges() << ",";
            //            insert_points(mesh, AreaMinimizer);
            file_franke_quad << mesh.n_vertices() + mesh.n_edges() << std::endl;
        }
        file_franke_quad.close();

        //--------------concave planes----------------------------------
        if (function == Franke2d)
        {
            filename_ = "./errors_poisson_Franke2D_concave.csv";
        }
        else
        {
            filename_ = "./condition_number_plane_concave.csv";
        }

        std::ofstream file_franke_concave(filename_);
        write_text_headers(file_franke_concave, function);
        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/grid/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            if (function == Franke2d)
            {
                write_all_laplace_data(mesh, file_franke_concave, function);
            }
            else
            {
                write_cond_numbers(mesh, file_franke_concave, function);
            }
            //            file_franke_concave << res << ",";
            //            insert_points(mesh, AreaMinimizer);
            //            res = inverse_mean_edgelenth(mesh);
            //            file_franke_concave << res << std::endl;
            //--------------------Dof--------------
            file_franke_concave << mesh.n_vertices() << ",";
            //            file_franke_concave << mesh.n_vertices()+mesh.n_edges() << ",";
            //            insert_points(mesh, AreaMinimizer);
            file_franke_concave << mesh.n_vertices() + mesh.n_edges()
                                << std::endl;
        }
        file_franke_concave.close();

        //--------------voronoi planes----------------------------------

        if (function == Franke2d)
        {
            filename_ = "./errors_poisson_Franke2D_Voronoi.csv";
        }
        else
        {
            filename_ = "./condition_number_plane_voronoi.csv";
        }

        std::ofstream file_franke_voronoi(filename_);
        write_text_headers(file_franke_voronoi, function);

        int voronoi_lvl = lvl;
        if (lvl > 6)
            voronoi_lvl = 6;
        for (int i = start_lvl; i < voronoi_lvl; i++)
        {
            std::string meshname =
                "../data/fernando_meshes/grid/clean/voronoi_" +
                std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            if (function == Franke2d)
            {
                write_all_laplace_data(mesh, file_franke_voronoi, function);
            }
            else
            {
                write_cond_numbers(mesh, file_franke_voronoi, function);
            }
            //            file_franke_voronoi << res << ",";
            //            insert_points(mesh, AreaMinimizer);
            //            res = inverse_mean_edgelenth(mesh);
            //            file_franke_voronoi << res << std::endl;

            //--------------------Dof--------------
            file_franke_voronoi << mesh.n_vertices() << ",";
            //            file_franke_voronoi << mesh.n_vertices()+mesh.n_edges() << ",";
            //            insert_points(mesh, AreaMinimizer);
            file_franke_voronoi << mesh.n_vertices() + mesh.n_edges()
                                << std::endl;
        }
        file_franke_voronoi.close();

        //        --------------triangle planes----------------------------------

        if (function == Franke2d)
        {
            filename_ = "./errors_poisson_Franke2D_triangle.csv";
        }
        else
        {
            filename_ = "./condition_number_plane_triangle.csv";
        }
        std::ofstream file_tri(filename_);
        write_text_headers(file_tri, function);

        for (int i = start_lvl; i < lvl; i++)
        {
            std::string meshname = "../data/fernando_meshes/grid/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            if (function == Franke2d)
            {
                write_all_laplace_data(mesh, file_tri, function);
            }
            else
            {
                write_cond_numbers(mesh, file_tri, function);
            }
            //            file_tri << res << ",";
            //            insert_points(mesh, AreaMinimizer);
            //            res = inverse_mean_edgelenth(mesh);
            //            file_tri << res << std::endl;
            //--------------------Dof--------------
            file_tri << mesh.n_vertices() << ",";
            //            file_franke_voronoi << mesh.n_vertices()+mesh.n_edges() << ",";
            //            insert_points(mesh, AreaMinimizer);
            file_tri << mesh.n_vertices() + mesh.n_edges() << std::endl;
        }
        file_tri.close();
        //
        //
        //        if (function == Franke2d)
        //        {
        //            filename_ = "./errors_poisson_Franke2D_triangle_reg.csv";
        //        }
        //        else
        //        {
        //            filename_ = "./condition_number_plane_triangle_reg.csv";
        //        }
        //        std::ofstream file_tri2(filename_);
        //        write_text_headers(file_tri2);
        //
        //
        //        for (int i = start_lvl; i < lvl; i++)
        //        {
        //            std::string meshname = "../data/fernando_meshes/grid/triangle_reg" +
        //                                   std::to_string(i) + ".off";
        //            mesh.read(meshname);
        //            double res = inverse_mean_edgelenth(mesh);
        //            if(function == Franke2d){
        //                write_all_laplace_data(mesh, file_tri2, function);
        //            }
        //            else{
        //                write_v0_degree_data(mesh, file_tri2, function);
        //            }
        //            file_tri2 << res << std::endl;
        //        }
        //        file_tri2.close();
    }
}

void write_nonstarshaped_convergence(int lvl)
{
    SurfaceMesh mesh;
    std::string filename_;
    std::ofstream file_franke_voronoi("error_2D_Franke_NonStarshaped.csv");
    file_franke_voronoi << "Linear, Quadratic,res" << std::endl;
    int c = 1;
    for (int i = 1; i < lvl; i++)
    {

        //    int i=5;
        std::string meshname = "../data/UMeshes/u" + std::to_string(c) + ".obj";
        mesh.read(meshname);
        //        filename_ = "../../NaiveProlongation/voronoi_" +
        //                    std::to_string(i) + "c.trip";
        //        filename_ = "../../NaiveProlongation/P" +
        //                    std::to_string(i) + ".mtx";
        //        double error = naive_Franke(mesh,filename_ );
        double res = inverse_mean_edgelenth(mesh);
        //        file_franke_voronoi << error <<",";

        file_franke_voronoi << solve_2D_AQAPoly_Poisson(mesh, Vertices, 1)
                            << ",";
        file_franke_voronoi << solve_2D_AQAPoly_Poisson(mesh, Edges, 2) << ",";
        file_franke_voronoi << res << "\n";
        c *= 2;
    }
    file_franke_voronoi.close();
}
void write_naive_convergence(int lvl)
{
    SurfaceMesh mesh;
    std::string filename_;
    std::ofstream file_franke_voronoi("error_2D_Franke_NaiveProlontaion.csv");
    file_franke_voronoi << "Naive Prolongation, Linear, Quadratic,res"
                        << std::endl;
    for (int i = 1; i < lvl; i++)
    {
        //    int i=5;
        std::string meshname =
            "../../NaiveProlongation/voronoi_" + std::to_string(i) + "c2.obj";
        mesh.read(meshname);
        //        filename_ = "../../NaiveProlongation/voronoi_" +
        //                    std::to_string(i) + "c.trip";
        filename_ = "../../NaiveProlongation/P" + std::to_string(i) + ".mtx";
        double error = naive_Franke(mesh, filename_);
        double res = inverse_mean_edgelenth(mesh);
        file_franke_voronoi << error << ",";

        file_franke_voronoi << solve_2D_AQAPoly_Poisson(mesh, Vertices, 1)
                            << ",";
        file_franke_voronoi << solve_2D_AQAPoly_Poisson(mesh, Edges, 2) << ",";
        file_franke_voronoi << res << "\n";
    }
    file_franke_voronoi.close();
}

//=============================================================================
int main()
{

    //    write_convergence_data_csv(poisson_SH, 6, 2, 0);
    //    write_convergence_data_csv(curvature);
    //      write_convergence_data_csv(cond_number, 3,0,0,2);
    //        write_naive_convergence(5);
    //        write_nonstarshaped_convergence(6);
        write_convergence_data_csv(Franke2d, 5);
//    write_convergence_data_csv(cond_number, 5);

    //
    //          write_mg_convergence_data_csv(Franke2d,6);
    //    write_convergence_data_csv(SH);
}
