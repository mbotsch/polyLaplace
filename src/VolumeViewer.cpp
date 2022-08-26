//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include <imgui.h>
#include <VolumeSubdivision.h>
#include <values.h>
#include "VolumeViewer.h"
#include "Volume/Eigenmodes.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "Volume/diffgeo_3D.h"
#include "Volume/AQAPoly_Laplacian_3D.h"
#include "VolumeMeshIO.h"
#include "Volume/LaplaceConstruction_3D.h"

//=============================================================================

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

bool VolumeViewer::load_mesh(const char *filename)
{
    if (VolumeMeshViewer::load_mesh(filename))
    {
        filename_ = filename;
        return true;
    }
    return false;
}

bool VolumeViewer::write_histogram(const char *filename)
{
    if (VolumeMeshViewer::write_histogram(filename))
    {
        filename_ = filename;
        return true;
    }
    return false;
}

void VolumeViewer::keyboard(int key, int code, int action, int mod)
{
    if (action != GLFW_PRESS) // only react on key press events
        return;

    switch (key)
    {

        case GLFW_KEY_B: {
            auto mesh = static_cast<VolumeMesh>(mesh_);
//
//            mat3 rotation(0.0);
//            for (int i = 0; i < 3; ++i)
//            {
//                for (int j = 0; j < 3; ++j)
//                {
//                    rotation(i, j) = modelview_matrix_(i, j);
//                }
//            }
//
//            {
//                auto box = mesh.get_bounding_sphere();
//
//                Vec3d cur = Vec3d(0.0);
//                Vec3d max = Vec3d(std::numeric_limits<double>::min());
//                Vec3d min = Vec3d(std::numeric_limits<double>::max());
//                for (auto v : mesh.vertices())
//                {
//                    auto vec = mesh.vertex(v);
//                    vec -= box.first;
//                    auto pmp_vec = pmp::vec3(vec[0], vec[1], vec[2]);
//                    pmp_vec = rotation * pmp_vec;
//                    vec = Vec3d(pmp_vec[0], pmp_vec[1], pmp_vec[2]);
//
//                    cur = vec;
//                    for (int i = 0; i < 3; ++i)
//                    {
//                        if (cur[i] > max[i])
//                        {
//                            max[i] = cur[i];
//                        }
//                        if (cur[i] < min[i])
//                        {
//                            min[i] = cur[i];
//                        }
//                    }
//                    mesh.set_vertex(v, vec);
//                }
//
//                Vec3d diag = max - min;
//                double maxD = *std::max_element(diag.data(), diag.data() + 3);
//
//                for (auto v : mesh.vertices())
//                {
//                    auto vec = mesh.vertex(v);
//                    vec /= maxD;
//                    mesh.set_vertex(v, vec);
//                }
//            }

            std::string name = filename_;
            int last_slash = name.find_last_of('/');
            int last_dot = name.find_last_of('.');
            name = name.substr(last_slash + 1, last_dot - last_slash - 1);

            file_manager_.writeFile(name + ".ovm", mesh);
            VolumeMeshIO output(name + ".HYBRID");
            output.write(mesh);
            break;
        }

        default: {
            VolumeMeshViewer::keyboard(key, code, action, mod);
            break;
        }
    }
}
//----------------------------------------------------------------------------

void VolumeViewer::process_imgui()
{
    VolumeMeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    static int laplace = 0;
    ImGui::RadioButton("Diamond", &laplace, 0);
    ImGui::RadioButton("Primal Laplace", &laplace, 1);
    ImGui::RadioButton("Dual Laplace", &laplace, 2);
    ImGui::RadioButton("Harmonic Basis", &laplace, 3);
    ImGui::RadioButton("Sandwich Laplace", &laplace, 4);
    ImGui::RadioButton("AQAPoly", &laplace, 6);
    ImGui::Text("Degree of Basis functions");

    ImGui::Spacing();
    static int degree = 2;
    ImGui::RadioButton("Linear", &degree, 1);
    ImGui::RadioButton("Quadratic", &degree, 2);
    ImGui::RadioButton("Cubic", &degree, 3);

    degree_ = degree;
    ImGui::Spacing();
    ImGui::Text("DoF to keep");

    ImGui::Spacing();
    static int coarsening = Edges;
    ImGui::RadioButton("Vertices", &coarsening, Vertices);
    ImGui::RadioButton("Edges", &coarsening, Edges);
    ImGui::RadioButton("Virtual Edges", &coarsening, Virtual_Edges);

    coarseningType_ = (CoarseDimension)coarsening;
    ImGui::Spacing();
    laplace_matrix = (unsigned)laplace;
    ImGui::Spacing();
    ImGui::Spacing();

    static int face_point = 0;
    ImGui::RadioButton("Area minimizer", &face_point, Quadratic_Areas_);
    ImGui::RadioButton("Face Centroid", &face_point, Face_Centroid);
    face_point_ = (unsigned)face_point;
    ImGui::Spacing();
    ImGui::Spacing();

    static int cell_point = 0;
    ImGui::RadioButton("Volume minimizer", &cell_point, Quadratic_Volume_);
    ImGui::RadioButton("Cell Centroid", &cell_point, Cell_Centroid_);
    cell_point_ = (unsigned)cell_point;
//    if(ImGui::Button("Print simpel Prolongation")){
//        Eigen::SparseMatrix<double>P;
//        setup_simple_3D_prolongation( mesh_,  P);
//    }
    if (ImGui::Button("Kugelize"))
    {
        kugelize(mesh_);
        update_mesh();
        double min_x = MAXFLOAT;
        double max_x = MINFLOAT;
        VolumeMesh::PointT min(0, 0, 0), max(0, 0, 0);
        for (auto v : mesh_.vertices())
        {
            VolumeMesh::PointT P = mesh_.vertex(v);
            if (P[0] < min_x)
            {
                min = P;
                min_x = P[0];
            }
            if (P[0] > max_x)
            {
                max_x = P[0];
                max = P;
            }
        }
        std::cout << " new Radius: " << (max_x - min_x) / 2.0 << std::endl;
        std::cout << "new center: " << 0.5 * (min + max) << std::endl;
    }
    if(ImGui::Button("Mean edge length")){
        double mean = 0.0;
        for(auto e : mesh_.edges()){
            mean += mesh_.length(e);
        }
        std::cout << "mean edge length: " << mean/mesh_.n_edges() << std::endl;
        std::cout << "Quad DOF: " << mesh_.n_edges()+mesh_.n_vertices() << std::endl;
        std::cout << "Linear DOF: " << mesh_.n_vertices() << std::endl;
    }
    if (ImGui::Button("RMSE Eigenvalues Sphere"))
    {
        Eigen::VectorXd evalues;
        solve_eigenvalue_problem(mesh_,filename_, evalues, laplace_matrix, face_point_,
                                 cell_point_,"",degree_,coarseningType_);
    }
    if (ImGui::Button("RMSE Franke Poisson System"))
    {
        if(laplace_matrix == 6)
        {
            std::ofstream timings_file("test.csv");
            if (degree_ == 1)
            {
                std::cout << "Linear " << std::endl;
//                solve_3D_AQAPoly_Poisson(filename_, Vertices,1);
                solve_3D_AQAPoly_Poisson_mg(filename_, timings_file, Vertices, 1,
                                            true);
            }
            else if (degree_ == 2)
            {
                std::cout << "Quadratic " << std::endl;
//                solve_3D_AQAPoly_Poisson(filename_, Edges,2);
                solve_3D_AQAPoly_Poisson_mg(filename_, timings_file, Edges, 2,
                                            true);

            }
        }else{
            solve_franke_poisson(mesh_, laplace_matrix, face_point_, cell_point_,2,filename_,coarseningType_);
        }
    }
    if (ImGui::Button("RMSE Franke Poisson System MG"))
    {
        if(laplace_matrix == 6)
        {
            std::ofstream timings_file("test.csv");
            if (degree_ == 1)
            {
                std::cout << "Linear " << std::endl;
                solve_3D_AQAPoly_Poisson_mg(filename_, timings_file, Vertices,
                                            1, false);
            }
            else if (degree_ == 2)
            {
                std::cout << "Quadratic " << std::endl;
                solve_3D_AQAPoly_Poisson_mg(filename_, timings_file, Edges, 2,
                                            false);
            }
        }else{
            solve_franke_poisson(mesh_, laplace_matrix, face_point_, cell_point_,2,filename_,coarseningType_);
        }
    }
    if (ImGui::Button("Linear Precision"))
    {
        if(laplace_matrix == 6){
            if (degree_ == 1)
            {
                solve_3D_AQAPoly_Poisson(filename_, Vertices, 1, 1);
            }else if(degree_ == 2)
            {
                solve_3D_AQAPoly_Poisson(filename_, Edges, 2, 1);
            }
        }else{
            solve_laplace_equation(mesh_, laplace_matrix, face_point_, cell_point_);
        }
    }
    if (ImGui::Button("Read Elasticity Sate"))
    {
        std::string filename = "./elastic_state_position.txt";
        std::vector<double> positions;
        std::ifstream myfile;
        myfile.open(filename);
        double num = 0.0;
        //keep storing values from the text file so long as data exists:
        while (myfile >> num) {
            positions.push_back(num);
        }
        int i = 0;
        for(auto v: mesh_.vertices()){
            std::cout << v << std::endl;
            VolumeMesh::PointT p (positions[i],positions[i+1],positions[i+2]);
            VolumeMesh::PointT r = mesh_.vertex(v);
            mesh_.set_vertex(v,r+p);
            i+=3;
        }
        update_mesh();
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Polyhedra!"))
    {
        if (ImGui::Button("virtual Points"))
        {
            VolumeSubdivision(mesh_).tetrahedra(face_point_, cell_point_);
            update_mesh();
        }

        if (ImGui::Button("irregular pyrmaids"))
        {
            VolumeSubdivision(mesh_).irregular_mesh(5);
            update_mesh();
        }

        if (ImGui::Button("full truncation"))
        {
            VolumeSubdivision(mesh_).full_truncation();
            update_mesh();
        }

        if (ImGui::Button("Quads"))
        {
            VolumeSubdivision(mesh_).quads();
            update_mesh();
        }

        if (ImGui::Button("Linear Subdivision"))
        {
            VolumeSubdivision(mesh_).linear_subdivision();
            update_mesh();
        }
    }
}

void VolumeViewer::mouse(int button, int action, int mods)
{
    VolumeMeshViewer::mouse(button, action, mods);
}
