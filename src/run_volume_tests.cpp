//=============================================================================
#include <iostream>
#include "Volume/Eigenmodes.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "VolumeMesh/VolumeSubdivision.h"
#include "Volume/HarmonicBasis.h"
//=============================================================================

enum LaplaceMethods {
    Diamond = 0,
    Harmonic = 1,
    PolySimpleLaplace = 2,
};

enum VolumePoints {
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints {
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

enum Function {
    Franke3d = 1,
};

//=============================================================================
double inverse_mean_edgelenth(VolumeMesh &mesh) {
    double avgLen = 0.;

    for (auto e: mesh.edges()) {

        avgLen += mesh.length(e);
    }

    return 1.0 / (avgLen / (double) mesh.n_edges());
}

void write_eigenmodes_test_results() {
    VolumeMesh mesh;
    std::string filename_ = "./errors_eigenmodes_3D.csv";

    std::ofstream file(filename_);

    Eigen::VectorXd evalues;
    std::string filename;
    file << "Hexa [BBA21],Hexa [BHBK20],Hexa [MKB08],Pyramids "
            "[BBA21],Pyramids [BHBK20],Pyramids [MKB08],Truncated "
            "[BBA21],Truncated [BHBK20], Truncated [MKB08],"
            "Tetrahedral_Coarse [BBA21],Tetrahedral_Coarse [BHBK20], Tetrahedral_Coarse [MKB08],"
            "Tetrahedral [BBA21],Tetrahedral [BHBK20],Tetrahedral [MKB08],"
         << std::endl;

    {
        std::string meshname =
                "../data/volume_meshes/unit_balls/Hexahedral_sphere.ovm";
        mesh.read(meshname);
        filename = "Hexahedra";
        double rmse;
        rmse = solve_eigenvalue_problem(mesh,  evalues, Diamond,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse = solve_eigenvalue_problem(mesh,  evalues, PolySimpleLaplace,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse =solve_harmonic_eigenvalue_problem(meshname,evalues,filename);
        file << rmse << ",";
    }
    {
        std::string meshname = "../data/volume_meshes/unit_balls/Pyramid_sphere.ovm";
        mesh.read(meshname);
        filename = "Pyramids";
        double rmse = solve_eigenvalue_problem(mesh, evalues, Diamond,
                                               Quadratic_Areas_,
                                               Quadratic_Volume_, filename);
        file << rmse << ",";
        rmse = solve_eigenvalue_problem(mesh,  evalues, PolySimpleLaplace,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse =solve_harmonic_eigenvalue_problem(meshname,evalues,filename);
        file << rmse << ",";
    }
    {
        std::string meshname = "../data/volume_meshes/unit_balls/Truncated_sphere.ovm";
        mesh.read(meshname);
        double rmse;
        filename = "Truncated";
        rmse = solve_eigenvalue_problem(mesh, evalues, Diamond,
                                               Quadratic_Areas_,
                                               Quadratic_Volume_, filename);
        file << rmse << ",";
        rmse = solve_eigenvalue_problem(mesh, evalues, PolySimpleLaplace,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse =solve_harmonic_eigenvalue_problem(meshname,evalues,filename);
        file << rmse << ",";
    }
    {
        std::string meshname =
                "../data/volume_meshes/unit_balls/tetrahedra_sphere_irreg.ovm";
        mesh.read(meshname);
        filename = "Tetrahedra Coarse";
        double rmse = solve_eigenvalue_problem(mesh,evalues, Diamond,
                                               Quadratic_Areas_,
                                               Quadratic_Volume_, filename);
        file << rmse << ",";
        rmse = solve_eigenvalue_problem(mesh, evalues, PolySimpleLaplace,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse =solve_harmonic_eigenvalue_problem(meshname,evalues,filename);
        file << rmse << ",";
    }
    {
        std::string meshname =
                "../data/volume_meshes/unit_balls/tetrahedra_sphere_reg.ovm";
        mesh.read(meshname);
        filename = "Tetrahedra";
        double rmse = solve_eigenvalue_problem(mesh,  evalues, Diamond,
                                               Quadratic_Areas_,
                                               Quadratic_Volume_, filename);
        file << rmse << ",";
        rmse = solve_eigenvalue_problem(mesh,  evalues, PolySimpleLaplace,
                                        Quadratic_Areas_, Quadratic_Volume_,
                                        filename);
        file << rmse << ",";
        rmse =solve_harmonic_eigenvalue_problem(meshname,evalues,filename);
        file << rmse << ",";
    }
    file.close();
}

void write_data(VolumeMesh &mesh,const std::string& meshname, int laplace,
                std::ofstream &error_file, int function) {

    if (function == Franke3d) {
        if (laplace == Diamond) {
            error_file << solve_franke_poisson(mesh, Diamond, Quadratic_Areas_,
                                               Quadratic_Volume_);
        } else if (laplace == PolySimpleLaplace) {
            error_file << solve_franke_poisson(mesh, PolySimpleLaplace, Quadratic_Areas_,
                                               Quadratic_Volume_);
        } else if(laplace == Harmonic){
            error_file << solve_3D_Franke_harmonic(meshname);
        }
        error_file << ",";
    } else {
        std::cout << " Function not implemented" << std::endl;
    }
}


void write_all_laplace_data(VolumeMesh &mesh, const std::string& meshname,
                            std::ofstream &file) {
    write_data(mesh, meshname, Diamond, file, Franke3d);
    write_data(mesh, meshname, PolySimpleLaplace, file, Franke3d);
    write_data(mesh, meshname, Harmonic, file, Franke3d);
}

void write_text_headers(Function function, std::ofstream &file_error) {
    if (function == Franke3d) {
        file_error << "[BBA21],[BHKB20],Harmonic, MEL" << std::endl;
    }
}

void write_3D_convergence_data_csv(Function function, int lvl_end = 6,
                                     int lvl_start = 2) {

    VolumeMesh mesh;
    std::string filename_, timings_filename_;

    // --------------hex cubes----------------------------------

    if (function == Franke3d) {
        filename_ = "./errors_poisson_Franke3D_hex.csv";
    }

    std::ofstream file_franke_hex(filename_);
    write_text_headers(function, file_franke_hex);

    for (int i = lvl_start; i < lvl_end; i++) {
        std::string meshname = "../data/volume_meshes/cubes/cube_hexahedron_" +
                               std::to_string(i) + ".ovm";
        mesh.read(meshname);
        std::cout << meshname << std::endl;
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d) {
            write_all_laplace_data(mesh,meshname, file_franke_hex);
        }
        file_franke_hex << res << std::endl;
    }
    file_franke_hex.close();

    //--------------pyramid cubes----------------------------------
    if (function == Franke3d) {
        filename_ = "./errors_poisson_Franke3D_pyramid.csv";
    }
    std::ofstream file_franke_pyramid(filename_);
    write_text_headers(function, file_franke_pyramid);

    for (int i = lvl_start; i < lvl_end; i++) {
        std::string meshname = "../data/volume_meshes/cubes/cube_pyramids_" +
                               std::to_string(i) + ".ovm";
        std::cout << meshname << std::endl;

        mesh.read(meshname);
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d) {
            write_all_laplace_data(mesh, meshname,file_franke_pyramid);
        }
        file_franke_pyramid << res << std::endl;
    }
    file_franke_pyramid.close();
    //--------------truncated cubes----------------------------------
    if (function == Franke3d) {
        filename_ = "./errors_poisson_Franke3D_truncated.csv";
    }
    std::ofstream file_franke_truncated(filename_);

    write_text_headers(function, file_franke_truncated);

    for (int i = lvl_start; i < lvl_end; i++) {
        std::string meshname = "../data/volume_meshes/cubes/cube_truncated_" +
                               std::to_string(i) + ".ovm";
        std::cout << meshname << std::endl;
        mesh.read(meshname);
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d) {
            write_all_laplace_data(mesh, meshname,file_franke_truncated);
        }
        file_franke_truncated << res << ",\n";
    }
    file_franke_truncated.close();

    //--------------voronoi cell cubes----------------------------------
    if (function == Franke3d) {
        filename_ = "./errors_poisson_Franke3D_voronoi.csv";
    }
    std::ofstream file_franke_voronoi(filename_);
    write_text_headers(function, file_franke_voronoi);
    for (int i = lvl_start - 1; i < lvl_end - 1; i++) {
        std::string meshname = "../data/volume_meshes/cubes/cube_voronoi_" +
                               std::to_string(i) + ".ovm";
        mesh.read(meshname);
        std::cout << meshname << std::endl;
        double res = inverse_mean_edgelenth(mesh);
        if (function == Franke3d) {
            write_all_laplace_data(mesh,meshname, file_franke_voronoi);
        }
        file_franke_voronoi << res << ",\n";
    }
    file_franke_voronoi.close();
}

//----------------------------------------------------------------------------

int main() {

    write_3D_convergence_data_csv(Franke3d, 6);
    write_eigenmodes_test_results();
}
//=============================================================================
