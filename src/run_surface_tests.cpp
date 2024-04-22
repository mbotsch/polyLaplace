//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================//=============================================================================
#include "../common_util.h"
#include <pmp/surface_mesh.h>
#include <pmp/bounding_box.h>
#include <pmp/algorithms/utilities.h>
#include <pmp/io/io.h>
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

enum Function
{
    poisson_SH = 0,
    Geodesics = 1,
    Franke2d = 2,
    SH = 3,
    curvature = 4,
    Eigenvalues = 5
};

LaplaceMode2D methods[] = {AlexaWardetzkyLaplace, deGoesLaplace,
                         PolySimpleLaplace, Diamond};

float AWLaplace[] = {2.0, 1.0, 0.5, 0.1};
float dGLaplace[] = {2.0, 1.0, 0.5, 0.1};

InsertedPoint vvm[] = {AreaMinimizer, TraceMinimizer};

std::string getModeName(LaplaceMode2D laplace)
{
    switch (laplace)
    {
        case PolySimpleLaplace:
            return "[BHKB20]";
        case AlexaWardetzkyLaplace:
            return "[AW11]";
        case Diamond:
            return "[BBA21]";
        case deGoesLaplace:
            return "[dGBD20]";
        case Harmonic:
            return "[MKB08]";
        default:
            return "";
    }
}

std::string getVertexName(InsertedPoint virt_vert_m)
{
    switch (virt_vert_m)
    {
        case Centroid_:
            return "Centroid";
        case AreaMinimizer:
            return "SqAreaMin";
        case TraceMinimizer:
            return "TraceMin";
        default:
            return "";
    }
}

// Different Mesh Types
enum class PlanarMeshes
{
    Triangle = 0,
    Quad = 1,
    Voronoi = 2,
    Concave = 3,
    Diamond_Quad = 4,
    Compressed_Hex = 5,
    Stretched_Hex = 6,
    Irregular_Polygon = 7,
};

enum class SphericalMeshes
{
    Triangle = 0,
    Quad = 1,
    Hex = 2,
    Voronoi = 3,
    Concave = 4,
    Noisy_Hex = 5,
    Noisy_Hex_NP = 6,
    Cos_Hex = 7,
    Acos_Hex = 8,
    Sin_Hex = 9,
};

DiffusionStep steps[] = {MaxDiagonal};

std::string getDiffusionStepName(DiffusionStep step)
{
    switch (step)
    {
        case MeanEdge:
            return "MeanEdge";
        case MaxEdge:
            return "MaxEdge";
        case MaxDiagonal:
            return "MaxDiagonal";
        default:
            return "";
    }
}

std::string getPlanarMeshName(PlanarMeshes pm)
{
    switch (pm)
    {
        case PlanarMeshes::Triangle:
            return "triangle";
        case PlanarMeshes::Quad:
            return "quad";
        case PlanarMeshes::Voronoi:
            return "voronoi";
        case PlanarMeshes::Concave:
            return "concave";
        case PlanarMeshes::Diamond_Quad:
            return "dia_quad";
        case PlanarMeshes::Compressed_Hex:
            return "comp_hex";
        case PlanarMeshes::Stretched_Hex:
            return "str_hex";
        case PlanarMeshes::Irregular_Polygon:
            return "irr_poly";
        default:
            return "";
    }
}

std::string getSphericalMeshName(SphericalMeshes sm)
{
    switch (sm)
    {
        case SphericalMeshes::Triangle:
            return "triangle";
        case SphericalMeshes::Quad:
            return "quad";
        case SphericalMeshes::Hex:
            return "hex";
        case SphericalMeshes::Voronoi:
            return "voronoi";
        case SphericalMeshes::Concave:
            return "concave";
        case SphericalMeshes::Noisy_Hex:
            return "hex_noise";
        case SphericalMeshes::Noisy_Hex_NP:
            return "hex2_noise";
        case SphericalMeshes::Cos_Hex:
            return "hex_cos";
        case SphericalMeshes::Acos_Hex:
            return "hex_acos";
        case SphericalMeshes::Sin_Hex:
            return "hex_sin";
        default:
            return "";
    }
}

std::string getPlanarMeshNameCSV(PlanarMeshes pm)
{
    switch (pm)
    {
        case PlanarMeshes::Triangle:
            return "Triangle";
        case PlanarMeshes::Quad:
            return "Quad";
        case PlanarMeshes::Voronoi:
            return "Voronoi";
        case PlanarMeshes::Concave:
            return "Concave";
        case PlanarMeshes::Diamond_Quad:
            return "QuadDiamond";
        case PlanarMeshes::Compressed_Hex:
            return "HexCompressed";
        case PlanarMeshes::Stretched_Hex:
            return "HexStretched";
        case PlanarMeshes::Irregular_Polygon:
            return "PolyIrregular";
        default:
            return "";
    }
}

std::string getSphericalMeshNameCSV(SphericalMeshes sm)
{
    switch (sm)
    {
        case SphericalMeshes::Triangle:
            return "Triangle";
        case SphericalMeshes::Quad:
            return "Quad";
        case SphericalMeshes::Hex:
            return "Hexagonal";
        case SphericalMeshes::Voronoi:
            return "Voronoi";
        case SphericalMeshes::Concave:
            return "Concave";
        case SphericalMeshes::Noisy_Hex:
            return "HexNoisy";
        case SphericalMeshes::Noisy_Hex_NP:
            return "HexNoisy2";
        case SphericalMeshes::Cos_Hex:
            return "HexAnisotropic";
        case SphericalMeshes::Acos_Hex:
            return "HexAnisotropic2";
        case SphericalMeshes::Sin_Hex:
            return "HexAnisotropic3";
        default:
            return "";
    }
}

void normalize(SurfaceMesh& mesh)
{
    for (auto v : mesh.vertices())
    {
        mesh.position(v) = normalize(mesh.position(v));
    }
}


void write_data(SurfaceMesh& mesh, int laplace, int virtualVertexMode,
                std::ofstream& file, int function,
                bool spheredist, DiffusionStep step,
                const std::vector<int>& indices, int l = 2, int m = 0)
{

    if (function == curvature)
    {
        Curvature analyzer(mesh, true);
        bool lumped = true;
        if (laplace == Diamond)
        {
            lumped = false;
        }
        file << analyzer.compute_curvature_error(laplace, virtualVertexMode, lumped);
        file << ",";
    }
    else if (function == SH)
    {
        bool lumped = true;
        file << rmse_sh(mesh, laplace, virtualVertexMode, lumped);
        file << ",";
    }
    else if (function == poisson_SH || function == Franke2d)
    {
        double error = solve_poisson_system(
            mesh, laplace, virtualVertexMode, function, l, m);
        file << error;
        file << ",";
    }
    else if (function == Geodesics)
    {
        double condition_number = 0;
        if (laplace == Harmonic)
        {
            file << 0.0;
            file << ",";
        }
        else
        {
            GeodesicsInHeat heat(mesh, laplace, virtualVertexMode, spheredist, !spheredist, step);
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();

            double error = 0;
            for (int n : indices)
            {
                error += heat.getDistance(n, dist, geodist, false);
            }
            error /= (float)indices.size();
            file << error;
            file << ",";


            std::string meshtype;
            if (spheredist)
            {
                meshtype = " sphere ";
            }
            else
            {
                meshtype = " plane ";
            }

            if (laplace == AlexaWardetzkyLaplace)
            {
                std::cout << "Distance deviation" + meshtype +
                                 "(AlexaWardetzky Laplace, l="
                          << poly_laplace_lambda_ << "): " << error << std::endl;
            }
            else if (laplace == deGoesLaplace)
            {
                std::cout << "Distance deviation" + meshtype + "(deGoes Laplace, l="
                          << deGoes_laplace_lambda_ << "): " << error << std::endl;
            }
            else
            {
                if (laplace == PolySimpleLaplace)
                {
                    std::cout << "Distance deviation" + meshtype + "(Polysimple ";
                }
                else if (laplace == Diamond)
                {
                    std::cout << "Distance deviation" + meshtype + "(Diamond ";
                }

                if (virtualVertexMode == AreaMinimizer)
                {
                    std::cout << "squared area minimizer): " << error << std::endl;
                }
                else if (virtualVertexMode == Centroid_)
                {
                    std::cout << "centroid): " << error << std::endl;
                }
                else
                {
                    std::cout << "trace minimizer): " << error << std::endl;
                }
            }
        }
    }
}

double inverse_mean_edge_length(SurfaceMesh& mesh)
{
    double avgLen = 0.;

    for (auto e : mesh.edges())
    {
        auto v0 = mesh.vertex(e, 0);
        auto v1 = mesh.vertex(e, 1);
        avgLen += distance(mesh.position(v0), mesh.position(v1));
    }

    return 1.0 / (avgLen / (double)mesh.n_edges());
}

void write_all_laplace_data(SurfaceMesh& mesh, std::ofstream& file,
                            int function, bool spheredist,
                            DiffusionStep st, int l = 2, int m = 2)
{
    //int num_samples = 100;
    std::vector<int> indices;
    if (function == Geodesics) {
        /*
		for (int i = 0; i < num_samples; ++i) {
			indices.push_back(rand() % mesh.vertices_size());
		}
		*/

        BoundingBox bb_ = bounds(mesh);
        bool is2D = bb_.min()[2] - bb_.max()[2] < 1e-6;
        Eigen::Vector3d ext = bb_.max() - bb_.min();
        int n = is2D ? 16 : 8;

        pmp::Point sample(0,0,0);
        for (int i = 0; i < n; i++) {
            sample[0] = ((i+0.5) / n) * ext[0] + bb_.min()[0];
            for (int j = 0; j < n; j++) {
                sample[1] = ((j+0.5) / n) * ext[1] + bb_.min()[1];

                if(is2D) {
                    double minDist = ext.norm() * 2;
                    int argMin;
                    for (auto v : mesh.vertices()) {
                        if (norm(mesh.position(v) - sample) < minDist){
                            minDist = norm(mesh.position(v) - sample);
                            argMin = v.idx();
                        }
                    }
                    indices.push_back(argMin);
                } else {
                    for (int k = 0; k < n; k++) {
                        sample[2] = ((k+0.5) / n) * ext[2] + bb_.min()[2];

                        double minDist = ext.norm() * 2;
                        int argMin;
                        for (auto v : mesh.vertices()) {
                            if (norm(mesh.position(v) - sample) < minDist){
                                minDist = norm(mesh.position(v) - sample);
                                argMin = v.idx();
                            }
                        }
                        indices.push_back(argMin);
                    }
                }
            }
        }

    }

    for (LaplaceMode2D mode : methods)
    {
        switch (mode)
        {
            case PolySimpleLaplace:
            case Diamond:
                for (auto & i : vvm)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << getVertexName(i) << std::endl;
                    write_data(mesh, mode, i, file, function,
                               spheredist, st, indices, l,
                               m);
                }
                break;
            case AlexaWardetzkyLaplace:
                for (float lambda : AWLaplace)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << lambda << std::endl;
                    poly_laplace_lambda_ = lambda;
                    write_data(mesh, mode, AreaMinimizer, file,
                               function, spheredist, st,
                               indices, l, m);
                }
                break;
            case deGoesLaplace:
                for (float lambda : dGLaplace)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << lambda << std::endl;
                    deGoes_laplace_lambda_ = lambda;
                    write_data(mesh, mode, AreaMinimizer, file,
                               function, spheredist, st,
                               indices, l, m);
                }
                break;
            default:
                std::cout << "=============================" << std::endl;
                std::cout << getModeName(mode) << std::endl;
                write_data(mesh, mode, AreaMinimizer, file,
                           function, spheredist, st,
                           indices, l, m);
        }
    }
}

void evaluate_conditionNr_per_mesh(SurfaceMesh& mesh, std::ofstream& file,
                                   bool generalized)
{
    double result;
    Eigen::Vector3d condnr;
    for (LaplaceMode2D mode : methods)
    {
        switch (mode)
        {
            case PolySimpleLaplace:
            case Diamond:
                for (auto & i : vvm)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << getVertexName(i) << std::endl;
                    result = condition_number(mesh, mode, i, condnr, generalized);
                    file << result << ",";
                    std::cout << "Condition Number: " << result << std::endl;
                }
                break;
            case AlexaWardetzkyLaplace:
                for (float lambda : AWLaplace)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << lambda << std::endl;
                    poly_laplace_lambda_ = lambda;
                    result =
                        condition_number(mesh, mode, AreaMinimizer, condnr, generalized);
                    file << result << ",";
                    std::cout << "Condition Number: " << result << std::endl;
                }
                break;
            case deGoesLaplace:
                for (float lambda : dGLaplace)
                {
                    std::cout << "=============================" << std::endl;
                    std::cout << getModeName(mode) << ": hyperparameter " << lambda << std::endl;
                    deGoes_laplace_lambda_ = lambda;
                    result =
                        condition_number(mesh, mode, AreaMinimizer, condnr, generalized);
                    file << result << ",";
                    std::cout << "Condition Number: " << result << std::endl;
                }
                break;
            default:
                std::cout << "=============================" << std::endl;
                std::cout << getModeName(mode) << std::endl;
                result =
                    condition_number(mesh, mode, AreaMinimizer, condnr, generalized);
                file << result << ",";
                std::cout << "Condition Number of stiffness matrix: " << result << std::endl;
                break;
        }
    }
}

void write_text_headers(std::ofstream& file_error,
                        const std::string& xval = ",MEL")
{
    std::ostringstream string_stream;
    for (LaplaceMode2D mode : methods)
    {
        switch (mode)
        {
            case AlexaWardetzkyLaplace:
                for (double lambda : AWLaplace)
                {
                    string_stream << getModeName(mode) << " l=" << lambda
                                  << ",";
                }
                break;
            case deGoesLaplace:
                for (double lambda : dGLaplace)
                {
                    string_stream << getModeName(mode) << " l=" << lambda
                                  << ",";
                }
                break;
            case PolySimpleLaplace:
            case Diamond:
                for (auto & i : vvm)
                {
                    string_stream << getModeName(mode)
                                  << "m=" << getVertexName(i) << ",";
                }
                break;
            default:
                string_stream << getModeName(mode) << ",";
                break;
        }
    }
    file_error << string_stream.str().substr(0, string_stream.str().size() - 1)
               << xval << std::endl;
}


void write_convergence_data_csv(Function function,
                                const std::vector<PlanarMeshes>& pm,
                                const std::vector<SphericalMeshes>& sm,
                                int lvl = 7, int start_lvl = 1, int l = 2,
                                int m = 0)
{
    SurfaceMesh mesh;
    std::string filename_;

    std::map<Function, std::string> filename_prefix{
        {Franke2d, "Franke"},
        {Geodesics, "Geodesics"},
        {poisson_SH, "SphericalHarmonics"},
        {curvature, "Curvature"},
        {SH, "SphericalHarmonicsRec"}};

    if (function == Franke2d)
    {
        bool sphere_dist = false;

        for (PlanarMeshes p : pm)
        {
            filename_ = "./Plane_" + filename_prefix[function] + "_" +
                        getPlanarMeshNameCSV(p) + ".csv";
            std::ofstream file(filename_);
            write_text_headers(file);

            std::string type_name = getPlanarMeshName(p);
            int adjlvl = lvl;
            if (p == PlanarMeshes::Voronoi)
            {
                type_name = "clean/" +  type_name;
                adjlvl = std::min(lvl, 6);
            }

            for (int i = start_lvl; i < adjlvl; i++)
            {
                std::string meshname = "../data/surface_meshes/grid/" +
                                       type_name + "_" + std::to_string(i) +
                                       ".obj";
                read(mesh, meshname);
                std::cout << meshname << std::endl;
                double res = inverse_mean_edge_length(mesh);
                write_all_laplace_data(mesh, file, function, sphere_dist,
                                       MaxDiagonal);
                //--------------------Inverse Mean Edge length--------------
                file << res << std::endl;
            }
            file.close();
        }
    }

    if (function == Geodesics)
    {
        for (DiffusionStep step : steps)
        {
            bool sphere_dist = false;
            for (PlanarMeshes p : pm)
            {
                filename_ = "./Plane_" + filename_prefix[function] +
                            getDiffusionStepName(step) + "_" +
                            getPlanarMeshNameCSV(p) + ".csv";
                std::ofstream file(filename_);
                write_text_headers(file);

                std::string type_name = getPlanarMeshName(p);
                int adjlvl = lvl;
                if (p == PlanarMeshes::Voronoi)
                {
                    type_name = "clean/" +  type_name;
                    adjlvl = std::min(lvl, 6);
                }

                for (int i = start_lvl; i < adjlvl; i++)
                {
                    std::string meshname = "../data/surface_meshes/grid/" +
                                           type_name + "_" + std::to_string(i) +
                                           ".obj";
                    read(mesh, meshname);
                    std::cout << meshname << std::endl;
                    double res = inverse_mean_edge_length(mesh);

                    std::ofstream iter_file;
                    write_all_laplace_data(mesh, file, function, sphere_dist,
                                           step);
                    //--------------------Inverse Mean Edge length--------------
                    file << res << std::endl;
                }
                file.close();
            }

            sphere_dist = true;

            for (SphericalMeshes s : sm)
            {
                filename_ = "./Sphere_" + filename_prefix[function] +
                            getDiffusionStepName(step) + "_" +
                            getSphericalMeshNameCSV(s) + ".csv";

                std::ofstream file(filename_);
                write_text_headers(file);

                std::string type_name = getSphericalMeshName(s);

                if (s == SphericalMeshes::Cos_Hex ||
                    s == SphericalMeshes::Acos_Hex ||
                    s == SphericalMeshes::Sin_Hex)
                {
                    type_name = "anispheres/" +  type_name;
                }
                for (int i = start_lvl; i < lvl; i++)
                {
                    std::string meshname = "../data/surface_meshes/sphere/" +
                                           type_name + "_" + std::to_string(i) +
                                           ".obj";
                    read(mesh, meshname);
                    std::cout << meshname << std::endl;
                    if (s != SphericalMeshes::Noisy_Hex)
                    {
                        normalize(mesh);
                    }
                    double res = inverse_mean_edge_length(mesh);
                    std::ofstream iter_file;
                    write_all_laplace_data(mesh, file, function, sphere_dist,
                                           step);
                    file << res << std::endl;
                }
                file.close();
            }
        }
    }

    if (function != Franke2d && function != Geodesics)
    {
        bool sphere_dist = true;

        for (SphericalMeshes s : sm)
        {
            filename_ = "./Sphere_" + filename_prefix[function] + "_" +
                        getSphericalMeshNameCSV(s) + ".csv";

            std::ofstream file(filename_);
            write_text_headers(file);

            std::string type_name = getSphericalMeshName(s);

            if (s == SphericalMeshes::Cos_Hex ||
                s == SphericalMeshes::Acos_Hex || s == SphericalMeshes::Sin_Hex)
            {
                type_name = "anispheres/" +  type_name;
            }
            for (int i = start_lvl; i < lvl; i++)
            {
                std::string meshname = "../data/surface_meshes/sphere/" +
                                       type_name + "_" + std::to_string(i) +
                                       ".obj";
                read(mesh, meshname);
                std::cout << meshname << std::endl;
                if (s != SphericalMeshes::Noisy_Hex)
                {
                    normalize(mesh);
                }
                double res = inverse_mean_edge_length(mesh);
                write_all_laplace_data(mesh, file, function, sphere_dist,
                                       MaxDiagonal, l, m);
                file << res << std::endl;
            }
            file.close();
        }

        if (function == poisson_SH)
        {
            filename_ = "./Sphere_" + filename_prefix[function] + "_" +
                        getSphericalMeshNameCSV(SphericalMeshes::Noisy_Hex_NP) +
                        ".csv";

            std::ofstream error_file(filename_);
            write_text_headers(error_file, ",Mean plane distance");
            for (int i = start_lvl; i < 9; i++)
            {
                std::string type_name =
                    getSphericalMeshName(SphericalMeshes::Noisy_Hex_NP);
                std::string meshname = "../data/surface_meshes/sphere/" +
                                       type_name + "_" + std::to_string(i) +
                                       ".obj";
                read(mesh, meshname);
                std::cout << meshname << std::endl;
                auto vpoint = mesh.get_vertex_property<Point>("v:point");
                double dist_sum = 0.0;
                for (auto f : mesh.faces())
                {
                    // fit plane to face
                    Eigen::MatrixXd poly(mesh.valence(f), 3);
                    int j = 0;
                    for (auto v : mesh.vertices(f))
                    {
                        Point p = vpoint[v];
                        poly(j, 0) = p[0];
                        poly(j, 1) = p[1];
                        poly(j, 2) = p[2];
                        j++;
                    }
                    Eigen::Vector3d n, o;
                    fit_plane_to_polygon(poly, n, o);
                    // compute mean distance to plane
                    double dist = 0.0;
                    for (auto v : mesh.vertices(f))
                    {
                        Eigen::Vector3d p(vpoint[v][0], vpoint[v][1],
                                          vpoint[v][2]);
                        dist += abs(n.dot(p - o));
                    }
                    dist /= (double)mesh.valence(f);
                    dist_sum += dist;
                }
                write_all_laplace_data(mesh, error_file, function, sphere_dist,
                                       MaxDiagonal, l, m);
                error_file << dist_sum / (double)mesh.n_faces() << std::endl;
            }
            error_file.close();
        }
    }
}


void write_conditionNR_refined(int start_lvl, int lvl, bool generalized,
                               const std::vector<PlanarMeshes>& pm,
                               const std::vector<SphericalMeshes>& sm)
{
    std::string stream_prefix = (generalized) ? "Generalized" : "Stiffness";

    //PLANAR MESHES

    for (PlanarMeshes p : pm)
    {
        std::ofstream error_file("./Plane_" + stream_prefix + "Refinement_" +
                                 getPlanarMeshNameCSV(p) + ".csv");
        write_text_headers(error_file, ",Inverse Edge Length");

        SurfaceMesh mesh;
        std::string meshname;
        std::string type_name = getPlanarMeshName(p);
        int adjlvl = lvl;
        if (p == PlanarMeshes::Voronoi)
        {
            type_name = "clean/" +  type_name;
            adjlvl = std::min(lvl, 6);
        }

        for (int i = start_lvl; i < adjlvl; i++)
        {
            meshname = "../data/surface_meshes/grid/" + type_name + "_" +
                       std::to_string(i) + ".obj";
            read(mesh, meshname);
            double res = inverse_mean_edge_length(mesh);
            evaluate_conditionNr_per_mesh(mesh, error_file, generalized);
            error_file << res << "\n";
        }
        error_file.close();
    }

    //SPHERICAL MESHES

    for (SphericalMeshes s : sm)
    {
        std::ofstream error_file("./Sphere_" + stream_prefix + "Refinement_" +
                                 getSphericalMeshNameCSV(s) + ".csv");
        write_text_headers(error_file, ",Inverse Edge Length");

        SurfaceMesh mesh;
        std::string meshname;

        std::string type_name = getSphericalMeshName(s);

        if (s == SphericalMeshes::Cos_Hex || s == SphericalMeshes::Acos_Hex ||
            s == SphericalMeshes::Sin_Hex)
        {
            type_name = "anispheres/" +  type_name;
        }

        for (int i = start_lvl; i < lvl; i++)
        {
            meshname = "../data/surface_meshes/sphere/" + type_name + "_" +
                       std::to_string(i) + ".obj";
            read(mesh, meshname);
            normalize(mesh);
            double res = inverse_mean_edge_length(mesh);
            evaluate_conditionNr_per_mesh(mesh, error_file, generalized);
            error_file << res << "\n";
        }
        error_file.close();
    }
}

void write_edgecollapse_per_mesh(SurfaceMesh& mesh, std::ofstream& error_file,
                                 bool generalized, bool sphere)
{
    Edge origin_edge;
    auto min = DBL_MAX;
    Point origin(0.0, 0.0, 0.0);
    if (sphere)
    {
        normalize(mesh);
        origin[0] += 1;
    }

    for (auto e : mesh.edges())
    {
        Vertex v1 = mesh.vertex(e, 1);
        Vertex v0 = mesh.vertex(e, 0);
        Point pos = 0.5 * (mesh.position(v0) + mesh.position(v1));

        if (norm(origin - pos) < min)
        {
            origin_edge = e;
            min = norm(origin - pos);
        }
    }
    Vertex v1 = mesh.vertex(origin_edge, 1);
    Vertex v0 = mesh.vertex(origin_edge, 0);
    std::cout << "edge length " << norm(mesh.position(v1) - mesh.position(v0))
              << std::endl;

    Point dir = (mesh.position(v1) - mesh.position(v0));
    double el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
    evaluate_conditionNr_per_mesh(mesh, error_file, generalized);
    error_file << el << "\n";
    while (norm((mesh.position(v1) - mesh.position(v0)) - 0.1 * dir) > 0.000001)
    {
        mesh.position(v0) += 0.1 * dir;
        std::cout << "edge length "
                  << norm(mesh.position(v1) - mesh.position(v0)) << std::endl;
        std::cout << "------------------------------" << std::endl;
        el = 1.0 / norm(mesh.position(v1) - mesh.position(v0));
        std::cout << mesh.position(v0) << std::endl;
        evaluate_conditionNr_per_mesh(mesh, error_file, generalized);
        error_file << el << "\n";
    }
    error_file.close();
}

void write_conditionNR_edgeshortening(bool generalized,
                                      const std::vector<PlanarMeshes>& pm,
                                      const std::vector<SphericalMeshes>& sm)
{
    int i = 3;

    SurfaceMesh mesh;
    std::string meshname;

    std::string stream_prefix = (generalized) ? "Generalized" : "Stiffness";

    // PLANAR MESHES

    for (PlanarMeshes p : pm)
    {
        std::ofstream error_file("./Plane_" + stream_prefix + "EdgeCollapse_" +
                                 getPlanarMeshNameCSV(p) + ".csv");

        write_text_headers(error_file, ",Inverse Edge length");
        std::string type_name = getPlanarMeshName(p);
        if (p == PlanarMeshes::Voronoi)
        {
            type_name = "clean/" +  type_name;
        }
        meshname = "../data/surface_meshes/grid/" + type_name + "_" +
                   std::to_string(i) + ".obj";
        read(mesh, meshname);
        write_edgecollapse_per_mesh(mesh, error_file, generalized, false);
    }

    // SPHERICAL MESHES

    for (SphericalMeshes s : sm)
    {
        std::ofstream error_file("./Sphere_" + stream_prefix + "EdgeCollapse_" +
                                 getSphericalMeshNameCSV(s) + ".csv");
        write_text_headers(error_file);
        std::string type_name = getSphericalMeshName(s);

        if (s == SphericalMeshes::Cos_Hex || s == SphericalMeshes::Acos_Hex ||
            s == SphericalMeshes::Sin_Hex)
        {
            type_name = "anispheres/" +  type_name;
        }

        meshname =
            "../data/surface_meshes/sphere/" + type_name + "_" +
            std::to_string(i + ((s == SphericalMeshes::Voronoi) ? 2 : 0)) +
            ".obj";
        read(mesh, meshname);
        write_edgecollapse_per_mesh(mesh, error_file, generalized, true);
    }
}

//=============================================================================
int main()
{
    std::vector<PlanarMeshes> pm_eval = {
        //PlanarMeshes::Triangle,      PlanarMeshes::Quad,
        PlanarMeshes::Stretched_Hex, PlanarMeshes::Compressed_Hex,
        PlanarMeshes::Concave,       PlanarMeshes::Voronoi,
        //PlanarMeshes::Diamond_Quad,  PlanarMeshes::Irregular_Polygon
    };

    std::vector<SphericalMeshes> sm_eval = {
        //SphericalMeshes::Triangle,  SphericalMeshes::Quad,
        SphericalMeshes::Hex,       SphericalMeshes::Cos_Hex,
        SphericalMeshes::Concave,   SphericalMeshes::Voronoi,
        //SphericalMeshes::Noisy_Hex, SphericalMeshes::Acos_Hex,
        //SphericalMeshes::Sin_Hex
    };

    std::vector<PlanarMeshes> pm_cond = {
        //PlanarMeshes::Triangle,      PlanarMeshes::Quad,
        PlanarMeshes::Stretched_Hex, PlanarMeshes::Compressed_Hex,
        PlanarMeshes::Concave,       PlanarMeshes::Voronoi,
        //PlanarMeshes::Diamond_Quad,  PlanarMeshes::Irregular_Polygon
    };

    std::vector<SphericalMeshes> sm_cond = {
        //SphericalMeshes::Triangle,  SphericalMeshes::Quad,
        SphericalMeshes::Hex,       SphericalMeshes::Cos_Hex,
        SphericalMeshes::Concave,   SphericalMeshes::Voronoi,
        //SphericalMeshes::Noisy_Hex, SphericalMeshes::Acos_Hex,
        //SphericalMeshes::Sin_Hex
    };

    std::vector<PlanarMeshes> pm_cond_es = {PlanarMeshes::Quad};

    std::vector<SphericalMeshes> sm_cond_es = {SphericalMeshes::Quad};

    // choose between AlexaWardetzky laplace implementations
    Herholz_version_ = true;

    write_convergence_data_csv(Geodesics, pm_eval, sm_eval, 6);
    //write_convergence_data_csv(curvature, pm_eval, sm_eval, 6);

    write_convergence_data_csv(Franke2d, pm_eval, sm_eval, 6);
    //write_convergence_data_csv(poisson_SH, pm_eval, sm_eval, 6, 1, 4, 2);
    //write_convergence_data_csv(poisson_SH, pm_eval, sm_eval, 6, 1, 3, -2);
    write_convergence_data_csv(poisson_SH, pm_eval, sm_eval, 6, 1, 3, 2);
    //write_convergence_data_csv(poisson_SH, pm_eval, sm_eval, 6, 1, 4, -2);
    //write_convergence_data_csv(SH, pm_eval, sm_eval, 6);

    write_conditionNR_edgeshortening(true, pm_cond_es, sm_cond_es);
    write_conditionNR_edgeshortening(false, pm_cond_es, sm_cond_es);

    write_conditionNR_refined(1, 5, true, pm_cond, sm_cond);
    write_conditionNR_refined(1, 6, false, pm_cond, sm_cond);

}
