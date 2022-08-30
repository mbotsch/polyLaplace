//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <fstream>
#include "Surface/[AW11]Laplace.h"
#include "Surface/Curvature.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/Poisson_System.h"
#include "Surface/SpectralProcessing.h"
#include "Surface/HarmonicBasis2D.h"

//=============================================================================

using namespace pmp;

//=============================================================================

enum LaplaceMethods {

    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 3,
    deGoesLaplace = 4,
    Harmonic = 5
};
enum InsertedPoint {
    Centroid = 0,
    AreaMinimizer = 2
};

enum Function {
    poisson_SH = 0,
    Franke2d = 2,
    curvature = 1,
    SH = 3,
};

void normalize(SurfaceMesh &mesh) {
    for (auto v: mesh.vertices()) {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

void write_data(SurfaceMesh &mesh, int laplace, std::ofstream &file, int function,
                int l = 2, int m = 0) {
    if (function == curvature) {
        Curvature analyzer(mesh, true);
        bool lumped = false;
        if (laplace == Diamond) {
            lumped = true;
        }
        file << analyzer.compute_curvature_error(laplace, AreaMinimizer, lumped);
        file << ",";
    } else if (function == SH) {
        bool lumped = true;
        if (laplace == Diamond) {
            lumped = false;
        } else {
            file << rmse_sh(mesh, laplace, AreaMinimizer, lumped);
        }

        file << ",";
    } else if (function == poisson_SH || function == Franke2d) {
        if (laplace == Harmonic) {
            file << solve_2D_Franke_harmonic(mesh);
        } else {

            double error =
                    solve_poisson_system(mesh, laplace, AreaMinimizer, function,
                                         l, m);
            file << error;
        }

        file << ",";
    }

}

double inverse_mean_edgelenth(SurfaceMesh &mesh) {
    double avgLen = 0.;

    for (auto e: mesh.edges()) {
        avgLen += mesh.edge_length(e);
    }

    return 1.0 / (avgLen / mesh.n_edges());
}

void write_all_laplace_data(SurfaceMesh &mesh, std::ofstream &file,
                            int function, int l = 2, int m = 2) {
    poly_laplace_lambda_ = 2.0;
    write_data(mesh, AlexaWardetzkyLaplace, file, function, l, m);
    disney_laplace_lambda_ = 1.0;
    write_data(mesh, deGoesLaplace, file, function, l, m);
    write_data(mesh, PolySimpleLaplace, file, function, l, m);
    write_data(mesh, Diamond, file, function, l, m);
}

void write_text_headers(std::ofstream &file_error, int function = Franke2d) {

        file_error
                << "[AW11],[dGBD20],[BHKB20],[BBA21],MEL"
                << std::endl;

}

double write_convergence_data_csv(Function function, int lvl = 7, int l = 2,
                                  int m = 0, int start_lvl = 1) {

    SurfaceMesh mesh;
    std::string filename_;
    std::string sh = "_Y_" + std::to_string(l) + std::to_string(m);

    if ((function != Franke2d)) {
        //--------------quad spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_quad.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_quad.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_quad.csv";
        }
        std::ofstream file(filename_);
        write_text_headers(file);
        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/sphere/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file, function, l, m);
            file << res << std::endl;
        }
        file.close();

        //--------------concave spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_concave.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_concave.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_concave.csv";
        }
        std::ofstream file_concave(filename_);
        write_text_headers(file_concave);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/sphere/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_concave, function, l, m);
            file_concave << res << std::endl;
        }
        file_concave.close();

        //--------------hex spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_hex.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_hex.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_hex.csv";
        }
        std::ofstream file_hex(filename_);
        write_text_headers(file_hex);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/sphere/hex_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_hex, function, l, m);
            file_hex << res << std::endl;
        }
        file_hex.close();

        //--------------triangle spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_triangle.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_triangle.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_triangle.csv";
        }
        std::ofstream file_tri(filename_);
        write_text_headers(file_tri);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/sphere/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_tri, function, l, m);
            file_tri << res << std::endl;
        }
        file_tri.close();
    } else {
        //--------------quad planes----------------------------------
        filename_ = "./errors_poisson_Franke2D_quad.csv";

        std::ofstream file_franke_quad(filename_);
        write_text_headers(file_franke_quad, function);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/grid/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_quad, function);
            //--------------------Inverse Mean Edge length--------------
            file_franke_quad << res << "\n";
        }
        file_franke_quad.close();

        //--------------concave planes----------------------------------
        filename_ = "./errors_poisson_Franke2D_concave.csv";

        std::ofstream file_franke_concave(filename_);
        write_text_headers(file_franke_concave, function);
        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/grid/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_concave, function);
            file_franke_concave << res << "\n";
        }
        file_franke_concave.close();

        //--------------voronoi planes----------------------------------
        filename_ = "./errors_poisson_Franke2D_Voronoi.csv";

        std::ofstream file_franke_voronoi(filename_);
        write_text_headers(file_franke_voronoi, function);

        int voronoi_lvl = lvl;
        if (lvl > 6)
            voronoi_lvl = 6;
        for (int i = start_lvl; i < voronoi_lvl; i++) {
            std::string meshname =
                    "../data/fernando_meshes/grid/clean/voronoi_" +
                    std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_voronoi, function);
            file_franke_voronoi << res << "\n";
        }
        file_franke_voronoi.close();

        //        --------------triangle planes----------------------------------

        filename_ = "./errors_poisson_Franke2D_triangle.csv";

        std::ofstream file_tri(filename_);
        write_text_headers(file_tri, function);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/fernando_meshes/grid/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);

            write_all_laplace_data(mesh, file_tri, function);
            file_tri << res << "\n";
            file_tri.close();
        }
    }
}

//=============================================================================
int main() {

    //    write_convergence_data_csv(poisson_SH, 6, 2, 0);
    //    write_convergence_data_csv(curvature);
    write_convergence_data_csv(Franke2d, 5);
    //    write_convergence_data_csv(SH);
}
