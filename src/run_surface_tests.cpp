//=============================================================================

#include <pmp/SurfaceMesh.h>
#include <fstream>
#include <cfloat>
#include "Surface/[AW11]Laplace.h"
#include "Surface/Curvature.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/Poisson_System.h"
#include "Surface/SpectralProcessing.h"
#include "Surface/GeodesicsInHeat.h"
#include "Surface/HarmonicBasis2D.h"

//=============================================================================

using namespace pmp;

//=============================================================================

enum LaplaceMethods {
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
    Harmonic = 4
};

enum InsertedPoint {
    Centroid = 0,
    AreaMinimizer = 2
};

enum Function {
    poisson_SH = 0,
    Geodesics = 1,
    Franke2d = 2,
    SH = 3,
    curvature = 4,
    Eigenvalues = 5
};

void normalize(SurfaceMesh &mesh) {
    for (auto v: mesh.vertices()) {
        mesh.position(v) = normalize(mesh.position(v));
    }
}

void write_data(SurfaceMesh &mesh, int laplace, std::ofstream &file, int function, bool spheredist,
                int l = 2, int m = 0) {
    if (function == curvature) {
        Curvature analyzer(mesh, true);
        bool lumped = true;
        if (laplace == Diamond) {
            lumped = false;
        }
        file << analyzer.compute_curvature_error(laplace, AreaMinimizer, lumped);
        file << ",";
    } else if (function == SH) {
        bool lumped = true;
        file << rmse_sh(mesh, laplace, AreaMinimizer, lumped);
        file << ",";
    } else if (function == poisson_SH || function == Franke2d) {
        double error =
                solve_poisson_system(mesh, laplace, AreaMinimizer, function,
                                     l, m);
        file << error;
        file << ",";
    } else if (function == Geodesics) {
        if (laplace == Harmonic) {
            file << 0.0;
            file << ",";
        } else {
            GeodesicsInHeat heat(mesh, laplace, AreaMinimizer,
                                 spheredist, !spheredist, MaxDiagonal);
            Eigen::VectorXd dist, geodist;
            int n = 0;
            heat.compute_geodesics();
            if (spheredist) {
                // pick vertex 0 for spheres
                n = 0;
            } else {
                //Center of plane
                Point p(0.0, 0.0, 0.0);
                double min_norm = 555555555555;
                for (auto v: mesh.vertices()) {
                    if (norm(mesh.position(v) - p) < min_norm) {
                        min_norm = norm(mesh.position(v) - p);
                        n = v.idx();
                    }
                }
            }
            double error = heat.getDistance(n, dist, geodist);
            file << error;
            file << ",";
        }
    }
}

double inverse_mean_edgelenth(SurfaceMesh &mesh) {
    double avgLen = 0.;

    for (auto e: mesh.edges()) {
        avgLen += mesh.edge_length(e);
    }

    return 1.0 / (avgLen / (double) mesh.n_edges());
}

void write_all_laplace_data(SurfaceMesh &mesh, std::ofstream &file,
                            int function, bool spheredist, int l = 2, int m = 2) {
    poly_laplace_lambda_ = 2.0;
    write_data(mesh, AlexaWardetzkyLaplace, file, function, spheredist, l, m);
    poly_laplace_lambda_ = 1.0;
    write_data(mesh, AlexaWardetzkyLaplace, file, function, spheredist, l, m);
    poly_laplace_lambda_ = 0.5;
    write_data(mesh, AlexaWardetzkyLaplace, file, function, spheredist, l, m);
    poly_laplace_lambda_ = 0.1;
    write_data(mesh, AlexaWardetzkyLaplace, file, function, spheredist, l, m);
    deGoes_laplace_lambda_ = 2.0;
    write_data(mesh, deGoesLaplace, file, function, spheredist, l, m);
    deGoes_laplace_lambda_ = 1.0;
    write_data(mesh, deGoesLaplace, file, function, spheredist, l, m);
    deGoes_laplace_lambda_ = 0.5;
    write_data(mesh, deGoesLaplace, file, function, spheredist, l, m);
    deGoes_laplace_lambda_ = 0.1;
    write_data(mesh, deGoesLaplace, file, function, spheredist, l, m);
    write_data(mesh, PolySimpleLaplace, file, function, spheredist, l, m);
    write_data(mesh, Diamond, file, function, spheredist, l, m);
    write_data(mesh, Harmonic, file, function, spheredist, l, m);

}

void write_text_headers(std::ofstream &file_error) {

    file_error
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08],MEL"
            << std::endl;

}

void write_kernelProbe_convergence_test_csv(int lvl = 6) {

    //--------------voronoi planes----------------------------------
    std::string filename_ = "./errors_harmonicKernel_Franke2D_Voronoi.csv";
    std::ofstream file_franke_voronoi(filename_);
    SurfaceMesh mesh;
    file_franke_voronoi
    << "nk=2 np=6,nk=2 np=8,nk=2 np=10,nk=5 np=15,nk=5 np=20,nk=5 np=25,nk=10 np=30,nk=10 np=40,nk=10 np=50,nk=20 np=60,nk=20 np=80,nk=20 np=100,MEL"
            << std::endl;
    int voronoi_lvl = lvl;
    if (lvl > 6)
        voronoi_lvl = 6;
    Eigen::SparseMatrix<double> S, M;
    for (int j = 1; j < voronoi_lvl; j++) {
        std::string meshname =
                "../data/surface_meshes/grid/clean/voronoi_" +
                std::to_string(j) + ".obj";
        mesh.read(meshname);
        std::cout << meshname << std::endl;

        int nv = (int) mesh.n_vertices();
        Eigen::VectorXd b(nv), B(nv);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        Eigen::MatrixXd X;
        for (auto v: mesh.vertices()) {
            b(v.idx()) = laplace_of_poisson_function(mesh.position(v), Franke2d);
        }
        int n = 4;
        int nKernel,nProbes;
        for (int nk = 0; nk < n; nk++) {
            for (int np = 0; np < 3; np++) {
                double error = 0.0;
                int kernel[4] = {2,5,10,20};
                int probes[3] = {3,4,5};
                nKernel = kernel[nk];
                nProbes = nKernel*probes[np];
                buildStiffnessAndMass2d(mesh, S, M, nKernel, nProbes);
                lump_matrix(M);
                B = M * b;
                S *= -1.0;
                // Set the constraints at the locked vertices to the evluation of the Franke function
                for (auto v: mesh.vertices()) {
                    if (mesh.is_boundary(v)) {
                        // right-hand side: fix boundary values with franke function of the vertices
                        B(v.idx()) = poisson_function(mesh.position(v), Franke2d);
                    }
                }

                // Adjust the right-hand-side to account for the locked nodes
                for (unsigned int i = 0; i < S.outerSize(); i++)
                    for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i);
                         iter; ++iter) {
                        Vertex row = pmp::Vertex(iter.row());
                        Vertex col = pmp::Vertex(iter.col());
                        if (!mesh.is_boundary(row) && mesh.is_boundary(col)) {
                            B[iter.row()] -= B[iter.col()] * iter.value();
                        }
                    }

                // Adjust the system matrix to account for the locked nodes
                for (unsigned int i = 0; i < S.outerSize(); i++)
                    for (Eigen::SparseMatrix<double>::InnerIterator iter(S, i);
                         iter; ++iter) {
                        Vertex row = pmp::Vertex(iter.row());
                        Vertex col = pmp::Vertex(iter.col());
                        if (mesh.is_boundary(row))
                            iter.valueRef() = iter.row() == iter.col() ? 1. : 0.;
                        else if (mesh.is_boundary(col))
                            iter.valueRef() = 0;
                    }

                solver.compute(S);
                Eigen::VectorXd x = solver.solve(B);
                for (auto v: mesh.vertices()) {
                    error += pow(
                            x[v.idx()] - poisson_function(mesh.position(v), Franke2d),
                            2.);
                }

                std::cout << "Franke RMSE error inner vertices: "
                          << sqrt(error / (double) mesh.n_vertices()) << std::endl;
                file_franke_voronoi << sqrt(error / (double) mesh.n_vertices()) << ",";
            }
        }
        double res = inverse_mean_edgelenth(mesh);
        file_franke_voronoi << res << "\n";
    }
    file_franke_voronoi.close();


    //--------------concave spheres----------------------------------
    filename_ = "./errors_harmonicKernel_SH32_hex.csv";
    std::ofstream file_concave(filename_);
    file_concave
            << "nk=2 np=6,nk=2 np=8,nk=2 np=10,nk=5 np=15,nk=5 np=20,nk=5 np=25,nk=10 np=30,nk=10 np=40,nk=10 np=50,nk=20 np=60,nk=20 np=80,nk=20 np=100,MEL"
            << std::endl;
    for (int j = 1; j < lvl; j++) {
        std::string meshname = "../data/surface_meshes/sphere/hex_" + std::to_string(j) + ".obj";
        mesh.read(meshname);
        std::cout << meshname << std::endl;

        int nv = (int) mesh.n_vertices();
        Eigen::VectorXd b(nv);
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        Eigen::MatrixXd X;
        int l = 3;
        int m = 2;
        for (auto v: mesh.vertices()) {
            b(v.idx()) = sphericalHarmonic(mesh.position(v), l, m);
        }
        int n = 4;
        int nKernel, nProbes;
        for (int nk = 0; nk < n; nk++) {
            for (int np = 0; np < 3; np++) {
                double error = 0.0;
                int kernel[4] = {2, 5, 10, 20};
                int probes[3] = {3, 4, 5};
                nKernel = kernel[nk];
                nProbes = nKernel * probes[np];
                buildStiffnessAndMass2d(mesh, S, M, nKernel, nProbes);
                lump_matrix(M);
                S *= -1.0;
                solver.analyzePattern(M);
                solver.factorize(M);
                X = solver.solve(S * b);
                if (solver.info() != Eigen::Success) {
                    std::cout << "Issue: " << solver.info() << std::endl;
                    std::cerr << "Could not solve linear system\n";
                }
                double eval = -l * (l + 1);
                error = (b - 1.0 / eval * X).transpose() * M *
                        (b - 1.0 / eval * X);
                error = sqrt(error / double(nv));
                file_concave << sqrt(error / (double) mesh.n_vertices()) << ",";
                std::cout << "error spherical harmonics " << "Y_" << l << "^" << m << ": " << error
                          << std::endl;
            }
        }
        double res = inverse_mean_edgelenth(mesh);
        file_concave << res << "\n";
    }
    file_concave.close();
}

void write_convergence_data_csv(Function function, int lvl = 7, int start_lvl = 1, int l = 2,
                                int m = 0) {

    SurfaceMesh mesh;
    std::string filename_;
    std::string sh = "_Y_" + std::to_string(l) + std::to_string(m);

    if (function == Franke2d || function == Geodesics) {
        bool sphere_dist = false;
        //--------------quad planes----------------------------------
        if (function == Franke2d) {
            filename_ = "./errors_poisson_Franke2D_quad.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_quad_plane.csv";

        }

        std::ofstream file_franke_quad(filename_);
        write_text_headers(file_franke_quad);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/grid/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_quad, function, sphere_dist);
            //--------------------Inverse Mean Edge length--------------
            file_franke_quad << res << "\n";
        }
        file_franke_quad.close();

        //--------------concave planes----------------------------------
        if (function == Franke2d) {
            filename_ = "./errors_poisson_Franke2D_concave.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_concave_plane.csv";
        }
        std::ofstream file_franke_concave(filename_);
        write_text_headers(file_franke_concave);
        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/grid/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_concave, function, sphere_dist);
            file_franke_concave << res << "\n";
        }
        file_franke_concave.close();

        //--------------voronoi planes----------------------------------
        if (function == Franke2d) {
            filename_ = "./errors_poisson_Franke2D_Voronoi.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_Voronoi_plane.csv";
        }
        std::ofstream file_franke_voronoi(filename_);
        write_text_headers(file_franke_voronoi);
        int voronoi_lvl = lvl;
        if (lvl > 6)
            voronoi_lvl = 6;
        for (int i = start_lvl; i < voronoi_lvl; i++) {
            std::string meshname =
                    "../data/surface_meshes/grid/clean/voronoi_" +
                    std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;

            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_franke_voronoi, function, sphere_dist);
            file_franke_voronoi << res << "\n";
        }
        file_franke_voronoi.close();

        //  --------------triangle planes----------------------------------


        if (function == Franke2d) {
            filename_ = "./errors_poisson_Franke2D_triangle.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_triangle_plane.csv";
        }
        std::ofstream file_tri(filename_);
        write_text_headers(file_tri);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/grid/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;

            double res = inverse_mean_edgelenth(mesh);

            write_all_laplace_data(mesh, file_tri, function, sphere_dist);
            file_tri << res << "\n";
        }
        file_tri.close();

    }
    if (function != Franke2d) {
        bool sphere_dist = true;
        //--------------quad spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_quad.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_quad.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_quad.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_quad_sphere.csv";
        }
        std::ofstream file(filename_);
        write_text_headers(file);
        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/sphere/quad_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file, function, sphere_dist, l, m);
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
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_concave_sphere.csv";
        }
        std::ofstream file_concave(filename_);
        write_text_headers(file_concave);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/sphere/concave_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_concave, function, sphere_dist, l, m);
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
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_hex_sphere.csv";
        }
        std::ofstream file_hex(filename_);
        write_text_headers(file_hex);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/sphere/hex_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_hex, function, sphere_dist, l, m);
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
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_triangle_sphere.csv";
        }
        std::ofstream file_tri(filename_);
        write_text_headers(file_tri);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/sphere/triangle_" +
                                   std::to_string(i) + ".obj";
            mesh.read(meshname);
            std::cout << meshname << std::endl;

            normalize(mesh);
            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_tri, function, sphere_dist, l, m);
            file_tri << res << std::endl;
        }
        file_tri.close();

        //--------------noisy hex spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_noisyHex.csv";
        } else if (function == curvature) {
            filename_ = "./errors_curvature_noisyHex.csv";
        } else if (function == SH) {
            filename_ = "./errors_SH_band_recreation_noisyHex.csv";
        } else if (function == Geodesics) {
            filename_ = "./errors_Geodesics_noisyHex_sphere.csv";
        }
        std::ofstream file_noise(filename_);
        write_text_headers(file_noise);

        for (int i = start_lvl; i < lvl; i++) {
            std::string meshname = "../data/surface_meshes/sphere/hex_noise_" +
                                   std::to_string(i) + ".off";
            mesh.read(meshname);
            std::cout << meshname << std::endl;

            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_noise, function, sphere_dist, l, m);
            file_noise << res << std::endl;
        }
        file_noise.close();

        //--------------noisy quad spheres----------------------------------
        if (function == poisson_SH) {
            filename_ = "./errors_poisson_SH" + sh + "_noisyQuad.csv";
        }
        std::ofstream file_noise_quad(filename_);
        write_text_headers(file_noise_quad);

        for (int i = start_lvl; i < 5; i++) {
            std::string meshname = "../data/surface_meshes/sphere/noiseQuad_" +
                                   std::to_string(i) + ".off";
            mesh.read(meshname);
            std::cout << meshname << std::endl;

            double res = inverse_mean_edgelenth(mesh);
            write_all_laplace_data(mesh, file_noise_quad, function, sphere_dist, l, m);
            file_noise_quad << res << std::endl;
        }
        file_noise_quad.close();
    }
}

void solve_eigenvalues_per_mesh(SurfaceMesh &mesh, std::ofstream &file,
                                const std::string &filename) {
    double rmse;
    poly_laplace_lambda_ = 2.0;
    rmse = solve_eigenvalue_problem(mesh, AlexaWardetzkyLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    poly_laplace_lambda_ = 1.0;
    rmse = solve_eigenvalue_problem(mesh, AlexaWardetzkyLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    poly_laplace_lambda_ = 0.5;
    rmse = solve_eigenvalue_problem(mesh, AlexaWardetzkyLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    poly_laplace_lambda_ = 0.1;
    rmse = solve_eigenvalue_problem(mesh, AlexaWardetzkyLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    deGoes_laplace_lambda_ = 2.0;
    rmse = solve_eigenvalue_problem(mesh, deGoesLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    deGoes_laplace_lambda_ = 1.0;
    rmse = solve_eigenvalue_problem(mesh, deGoesLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    deGoes_laplace_lambda_ = 0.5;
    rmse = solve_eigenvalue_problem(mesh, deGoesLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    deGoes_laplace_lambda_ = 0.1;
    rmse = solve_eigenvalue_problem(mesh, deGoesLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    rmse = solve_eigenvalue_problem(mesh, PolySimpleLaplace, AreaMinimizer, filename);
    file << rmse << ",";
    rmse = solve_eigenvalue_problem(mesh, Diamond, AreaMinimizer, filename);
    file << rmse << ",";
    rmse = solve_eigenvalue_problem(mesh, Harmonic, AreaMinimizer, filename);
    file << rmse << ",";
}

void evaluate_conditionNr_per_mesh(SurfaceMesh &mesh, std::ofstream &file) {
    double result;
    Eigen::Vector3d condnr;
    poly_laplace_lambda_ = 2.0;
    result = condition_number(mesh, AlexaWardetzkyLaplace, AreaMinimizer, condnr);
    file << result << ",";
    poly_laplace_lambda_ = 1.0;
    result = condition_number(mesh, AlexaWardetzkyLaplace, AreaMinimizer, condnr);
    file << result << ",";
    poly_laplace_lambda_ = 0.5;
    result = condition_number(mesh, AlexaWardetzkyLaplace, AreaMinimizer, condnr);
    file << result << ",";
    poly_laplace_lambda_ = 0.1;
    result = condition_number(mesh, AlexaWardetzkyLaplace, AreaMinimizer, condnr);
    file << result << ",";
    deGoes_laplace_lambda_ = 2.0;
    result = condition_number(mesh, deGoesLaplace, AreaMinimizer, condnr);
    file << result << ",";
    deGoes_laplace_lambda_ = 1.0;
    result = condition_number(mesh, deGoesLaplace, AreaMinimizer, condnr);
    file << result << ",";
    deGoes_laplace_lambda_ = 0.5;
    result = condition_number(mesh, deGoesLaplace, AreaMinimizer, condnr);
    file << result << ",";
    deGoes_laplace_lambda_ = 0.1;
    result = condition_number(mesh, deGoesLaplace, AreaMinimizer, condnr);
    file << result << ",";
    result = condition_number(mesh, PolySimpleLaplace, AreaMinimizer, condnr);
    file << result << ",";
    result = condition_number(mesh, Diamond, AreaMinimizer, condnr);
    file << result << ",";
    result = condition_number(mesh, Harmonic, AreaMinimizer, condnr);
    file << result << ",";
}

void write_evalue_test_results(int lvl = 7) {
    std::string filename_;
    std::ofstream file(filename_);

    std::string filename;
    std::ofstream error_file("evalue_error_quad.csv");
    error_file
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08]\n";
    for (int i = 2; i < lvl; i++) {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../data/surface_meshes/sphere/quad_" +
                   std::to_string(i) + ".obj";
        filename = "quad_" + std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file, filename);
        error_file << "\n";
    }
    error_file.close();

    std::ofstream error_file_hex("evalue_error_hex.csv");
    error_file_hex
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08]\n";
    for (int i = 2; i < lvl; i++) {
        SurfaceMesh mesh;
        std::string meshname;
        meshname = "../data/surface_meshes/sphere/hex_" + std::to_string(i) +
                   ".obj";
        filename = "hex_" + std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_hex, filename);
        error_file_hex << "\n";
    }

    error_file_hex.close();
    std::ofstream error_file_triangle("evalue_error_triangle.csv");
    error_file_triangle
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08]\n";
    for (int i = 2; i < lvl; i++) {
        SurfaceMesh mesh;
        std::string meshname;
        int res = i;
        meshname = "../data/surface_meshes/sphere/triangle_" +
                   std::to_string(res) + ".obj";
        filename = "triangle_" + std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_triangle, filename);
        error_file_triangle << "\n";
    }
    error_file_triangle.close();

    std::ofstream error_file_concave("evalue_error_concave.csv");
    error_file_concave
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08]\n";
    for (int i = 2; i < lvl; i++) {
        SurfaceMesh mesh;
        std::string meshname;
        int res = i;
        meshname = "../data/surface_meshes/sphere/concave_" +
                   std::to_string(res) + ".obj";
        filename = "concave_" + std::to_string(i);
        mesh.read(meshname);
        normalize(mesh);
        solve_eigenvalues_per_mesh(mesh, error_file_concave, filename);
        error_file_concave << "\n";
    }
    error_file_concave.close();
}


void write_conditionNR_test_results() {
    int i = 3;
    std::string filename_;
    std::ofstream file(filename_);

    std::string filename;
    std::ofstream error_file("conditionNr_error_quad.csv");
    error_file
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08],Inverse Edge lenght\n";
    SurfaceMesh mesh;
    std::string meshname;
    meshname = "../data/surface_meshes/grid/quad_" +
               std::to_string(i) + ".obj";
    mesh.read(meshname);
    Edge origin_edge;
    auto min = DBL_MAX;
    for (auto e: mesh.edges()) {
        Vertex v1 = mesh.vertex(e, 1);
        Vertex v0 = mesh.vertex(e, 0);
        Point p = 0.5 * (mesh.position(v0) + mesh.position(v1));
        Point origin(0.0, 0.0, 0.0);
        if (norm(origin - p) < min) {
            origin_edge = e;
            min = norm(origin - p);
        }
    }
    Vertex v1 = mesh.vertex(origin_edge, 1);
    Vertex v0 = mesh.vertex(origin_edge, 0);
    std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;

    Point dir = (mesh.position(v1) - mesh.position(v0));
    double el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
    evaluate_conditionNr_per_mesh(mesh, error_file);
    error_file << el << "\n";
    while (norm((mesh.position(v1) - mesh.position(v0)) - 0.1 * dir) > 0.000001) {
        mesh.position(v0) += 0.1 * dir;
        std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;
        std::cout << "------------------------------" << std::endl;
        el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
        std::cout << mesh.position(v0) << std::endl;
        evaluate_conditionNr_per_mesh(mesh, error_file);
        error_file << el << "\n";
    }
    error_file.close();

    std::ofstream error_file_voronoi("conditionNr_error_voronoi.csv");
    error_file_voronoi
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08],Inverse Edge lenght\n";

    meshname = "../data/surface_meshes/grid/voronoi_" + std::to_string(i) +
               ".obj";
    mesh.read(meshname);
    min = DBL_MAX;
    for (auto e: mesh.edges()) {
        v1 = mesh.vertex(e, 1);
        v0 = mesh.vertex(e, 0);
        Point p = 0.5 * (mesh.position(v0) + mesh.position(v1));
        Point origin(0.0, 0.0, 0.0);
        if (norm(origin - p) < min) {
            origin_edge = e;
            min = norm(origin - p);
        }
    }
    v1 = mesh.vertex(origin_edge, 1);
    v0 = mesh.vertex(origin_edge, 0);
    std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;

    dir = (mesh.position(v1) - mesh.position(v0));
    el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
    evaluate_conditionNr_per_mesh(mesh, error_file_voronoi);
    error_file_voronoi << el << "\n";
    while (norm((mesh.position(v1) - mesh.position(v0)) - 0.1 * dir) > 0.000001) {
        mesh.position(v0) += 0.1 * dir;
        std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;
        std::cout << "------------------------------" << std::endl;
        el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
        std::cout << mesh.position(v0) << std::endl;
        evaluate_conditionNr_per_mesh(mesh, error_file_voronoi);
        error_file_voronoi << el << "\n";
    }

    error_file_voronoi.close();
    std::ofstream error_file_triangle("conditionNr_error_triangle.csv");
    error_file_triangle
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08],Inverse Edge lenght\n";

    meshname = "../data/surface_meshes/grid/triangle_" +
               std::to_string(i) + ".obj";
    filename = "triangle_" + std::to_string(i);
    mesh.read(meshname);
    min = DBL_MAX;
    for (auto e: mesh.edges()) {
        v1 = mesh.vertex(e, 1);
        v0 = mesh.vertex(e, 0);
        Point p = 0.5 * (mesh.position(v0) + mesh.position(v1));
        Point origin(0.0, 0.0, 0.0);
        if (norm(origin - p) < min) {
            origin_edge = e;
            min = norm(origin - p);
        }
    }
    v1 = mesh.vertex(origin_edge, 1);
    v0 = mesh.vertex(origin_edge, 0);
    std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;

    dir = (mesh.position(v1) - mesh.position(v0));
    el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
    evaluate_conditionNr_per_mesh(mesh, error_file_triangle);
    error_file_triangle << el << "\n";
    while (norm((mesh.position(v1) - mesh.position(v0)) - 0.1 * dir) > 0.000001) {
        mesh.position(v0) += 0.1 * dir;
        std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;
        std::cout << "------------------------------" << std::endl;
        el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
        std::cout << mesh.position(v0) << std::endl;
        evaluate_conditionNr_per_mesh(mesh, error_file_triangle);
        error_file_triangle << el << "\n";
    }

    error_file_triangle.close();

    std::ofstream error_file_concave("conditionNr_error_concave.csv");
    error_file_concave
            << "[AW11] l=2,[AW11] l=1,[AW11] l=0.5,[AW11] l=0.1,[dGBD20] l=2,[dGBD20] l=1,[dGBD20] l=0.5,[dGBD20] l=0.1,[BHKB20],[BBA21],[MKB08],Inverse Edge lenght\n";

    meshname = "../data/surface_meshes/grid/concave_" +
               std::to_string(i) + ".obj";
    mesh.read(meshname);
    min = DBL_MAX;
    for (auto e: mesh.edges()) {
        v1 = mesh.vertex(e, 1);
        v0 = mesh.vertex(e, 0);
        Point p = 0.5 * (mesh.position(v0) + mesh.position(v1));
        Point origin(0.0, 0.0, 0.0);
        if (norm(origin - p) < min) {
            origin_edge = e;
            min = norm(origin - p);
        }
    }
    v1 = mesh.vertex(origin_edge, 1);
    v0 = mesh.vertex(origin_edge, 0);
    std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;

    dir = (mesh.position(v1) - mesh.position(v0));
    el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
    evaluate_conditionNr_per_mesh(mesh, error_file_concave);
    error_file_concave << el << "\n";
    while (norm((mesh.position(v1) - mesh.position(v0)) - 0.1 * dir) > 0.000001) {
        mesh.position(v0) += 0.1 * dir;
        std::cout << "edge length" << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;
        std::cout << "------------------------------" << std::endl;
        el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
        std::cout << mesh.position(v0) << std::endl;
        evaluate_conditionNr_per_mesh(mesh, error_file_concave);
        error_file_concave << el << "\n";
    }
    error_file_concave.close();


}

//=============================================================================
int main() {
    // choose between AlexaWardetzky laplace implementations
    Herholz_version_ = true;
    write_convergence_data_csv(Geodesics, 6);
    write_convergence_data_csv(curvature, 6);
    write_convergence_data_csv(Franke2d, 6);
    write_convergence_data_csv(poisson_SH, 6, 1, 4,2 );
    write_convergence_data_csv(poisson_SH, 6, 1, 3,-2);
    write_convergence_data_csv(poisson_SH, 6, 1, 3, 2);
    write_convergence_data_csv(poisson_SH, 6, 1, 4,-2);
    write_convergence_data_csv(SH, 6);
    write_evalue_test_results(6);
    write_conditionNR_test_results();
    write_kernelProbe_convergence_test_csv(6);
}
