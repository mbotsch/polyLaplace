//=============================================================================
// Copyright 2024 Astrid Bunge, Sven Wagner, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "../common_util.h"
#include <pmp/visualization/mesh_viewer.h>
#include "Surface/Smoothing.h"
#include <imgui.h>
#include "Surface/[AW11]Laplace.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/GeodesicsInHeat.h"
#include "Surface/Curvature.h"
#include "Surface/Parameterization.h"
#include "Surface/PolySimpleLaplace.h"
#include "Surface/PolySmoothing.h"
#include "Surface/SpectralProcessing.h"
using namespace pmp;

//=============================================================================

class Viewer : public MeshViewer
{
public:
    Viewer(const char* title, int width, int height);
    void load_mesh(const char* filename) override;

protected:
    void process_imgui() override;
    void draw(const std::string& _draw_mode) override;
    void color_code_condition_numbers(int laplace, int min_point);

private:
    Smoothing smoother_;
    double min_cond = -1, max_cond = -1;
    bool show_uv_layout_;
    int mesh_index_ = 0;
};


void Viewer::load_mesh(const char* filename)
{
    min_cond = -1;
    max_cond = -1;
    MeshViewer::load_mesh(filename);
}

Viewer::Viewer(const char* title, int width, int height)
    : MeshViewer(title, width, height), smoother_(mesh_)
{
    set_draw_mode("Hidden Line");
    crease_angle_ = 0.0;
    show_uv_layout_ = false;
}

void Viewer::process_imgui()
{
    if (ImGui::CollapsingHeader("Load Mesh", ImGuiTreeNodeFlags_DefaultOpen))
    {
        const char* listbox_items[] = {"dual_Bunny.off", "fandisk.obj",
                                       "fertility.obj", "gear.obj",
                                       "horse.obj", "rockerArm.obj"};
        int listbox_item_current = mesh_index_;
        ImGui::ListBox(" ", &listbox_item_current, listbox_items,
                       IM_ARRAYSIZE(listbox_items), 5);
        if (listbox_item_current != mesh_index_ || mesh_.n_vertices() == 0)
        {
            std::stringstream ss;
            ss << DATA_PATH << listbox_items[listbox_item_current];

            load_mesh(ss.str().c_str());

            mesh_index_ = listbox_item_current;
        }
        ImGui::Spacing();
        ImGui::Spacing();
    }
    // initialize settings
    static int laplace = 0;
    static int min_point = 2;

    MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Polygon Laplace",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::RadioButton("Alexa & Wardetzky Laplace", &laplace, 1);
        ImGui::RadioButton("deGoes Laplace", &laplace, 3);
        ImGui::RadioButton("Virtual Refinement Laplace", &laplace, 0);
        ImGui::RadioButton("Diamond Laplace", &laplace, 2);
        ImGui::RadioButton("Harmonic Laplace", &laplace, 4);

        ImGui::Spacing();
        if (laplace == 0 || laplace == 2)
        {
            ImGui::Text("Choose your minimizing point ");

            ImGui::Spacing();

            ImGui::RadioButton("Naïve (Centroid)", &min_point, 0);
            ImGui::RadioButton("Simple (Area Minimizer)", &min_point, 2);
            ImGui::RadioButton("Robust (Trace Minimizer)", &min_point, 3);

            ImGui::Spacing();
        }
        if (laplace == 1)
        {
            ImGui::PushItemWidth(100);
            ImGui::SliderFloat("Lambda", &poly_laplace_lambda_, 0.01, 3.0);
            ImGui::PopItemWidth();
        }
        else if (laplace == 3)
        {
            ImGui::PushItemWidth(100);
            ImGui::SliderFloat("Lambda", &deGoes_laplace_lambda_, 0.01, 3.0);
            ImGui::PopItemWidth();
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Make it robust",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Indent(10);

        if (ImGui::Button("Mesh Optimization"))
        {
            SmoothingConfigs oConf(25, false, false, false, false);
            PolySmoothing polySmoothing(mesh_, oConf);
            polySmoothing.optimize(5);
            update_mesh();
        }

        if (ImGui::Button("Color Code Condition Number"))
        {
            color_code_condition_numbers(laplace, min_point);
            renderer_.set_specular(0);
            update_mesh();
            set_draw_mode("Hidden Line");
        }

        ImGui::Unindent(10);
    }

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Applications", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Indent(10);

        if (ImGui::Button("Mean Curvature"))
        {
            Curvature analyzer(mesh_, false);
            analyzer.visualize_curvature(laplace, min_point, true);
            renderer_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
            show_uv_layout_ = false;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        static float timestep = 0.01;
        float lb = 0.001;
        float ub = 1.0;
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("TimeStep", &timestep, lb, ub);
        ImGui::PopItemWidth();

        if (ImGui::Button("Smoothing"))
        {
            Scalar dt = timestep * radius_ * radius_;
            try
            {
                smoother_.implicit_smoothing(dt, laplace, min_point, true);
            }
            catch (const SolverException& e)
            {
                std::cerr << e.what() << std::endl;
            }

            update_mesh();
            set_draw_mode("Hidden Line");
            show_uv_layout_ = false;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Geodesics"))
        {
            GeodesicsInHeat heat(mesh_, laplace, min_point, false, false,
                                 DiffusionStep(2));
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance(0, dist, geodist);

            update_mesh();
            renderer_.use_checkerboard_texture();
            set_draw_mode("Texture");
            show_uv_layout_ = false;
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();

        if (ImGui::Button("Parametrization"))
        {
            try
            {
                Parameterization param(mesh_);
                param.harmonic(laplace, min_point);
            }
            catch (const std::exception& e)
            {
                std::cerr << e.what() << std::endl;
            }

            update_mesh();
            renderer_.use_checkerboard_texture();
            set_draw_mode("Texture");
            show_uv_layout_ = true;
        }

        ImGui::Unindent(10);
    }
}

void Viewer::color_code_condition_numbers(int laplace, int min_point)
{
    auto face_color = mesh_.face_property<Color>("f:color");
    auto face_cond = mesh_.face_property<double>("f:condition");

    // compute global condition number
    Eigen::Vector3d values;
    double cond = condition_number(mesh_, laplace, min_point, values, false);
    std::cout << "Condition Number: " << cond << std::endl;

    // compute per-face condition numbers
    Eigen::MatrixXd Si;
    Eigen::VectorXd w;
    Eigen::Vector3d p;
    Eigen::MatrixXd poly;
    for (Face f : mesh_.faces())
    {
        // collect polygon vertices
        get_polygon_from_face(mesh_, f, poly);

        // compute weights for the polygon
        if (min_point == Centroid_)
        {
            int val = (int)poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else if (min_point == AreaMinimizer)
        {
            find_area_minimizer_weights(poly, w);
        }
        else
        {
            find_trace_minimizer_weights(poly, w);
        }
        Eigen::Vector3d min;

        // compute virtual vertex
        min = poly.transpose() * w;

        // get per-face Laplace matrix
        localCotanMatrix(poly, min, w, Si);

        // compute per-face condition number
        Eigen::SelfAdjointEigenSolver<Eigen::SparseMatrix<double>> eigs(Si);
        face_cond[f] =
            eigs.eigenvalues()[Si.rows() - 1] / eigs.eigenvalues()[1];
    }

    // compute max/min condition number for color coding
    if (max_cond == -1 && min_cond == -1)
    {
        std::vector<double> cond_numbers;
        for (auto f : mesh_.faces())
        {
            cond_numbers.push_back(face_cond[f]);
        }
        std::ranges::sort(cond_numbers);
        max_cond = cond_numbers[int(0.99 * mesh_.n_faces())];
        min_cond = cond_numbers[0];
    }

    // turn condition number into color
    for (auto f : mesh_.faces())
    {
        auto good_col = Color(0.39, 0.74, 1); // Turquoise (good)
        auto ok_col = Color(1, 0.74, 0);      // Orange (okay)
        auto bad_col = Color(1, 0.0, 1);      // Purple (bad)

        double blend = fmin(
            1.0, fmax(0.0, (face_cond[f] - min_cond) / (max_cond - min_cond)));
        face_color[f] = (blend < 0.5)
                            ? (1 - blend) * good_col + blend * ok_col
                            : (1 - blend) * ok_col + blend * bad_col;
    }
}

void Viewer::draw(const std::string& draw_mode)
{
    // draw the mesh
    renderer_.draw(projection_matrix_, modelview_matrix_, draw_mode);

    // draw uv layout
    if (draw_mode == "Texture" && show_uv_layout_)
    {
        // clear depth buffer
        glClear(GL_DEPTH_BUFFER_BIT);

        // setup viewport
        GLint size = std::min(width(), height()) / 4;
        glViewport(width() - size - 1, height() - size - 1, size, size);

        // setup matrices
        mat4 P = ortho_matrix(0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 1.0f);
        mat4 M = mat4::identity();

        // draw mesh once more
        renderer_.draw(P, M, "Texture Layout");

        // reset viewport
        glViewport(0, 0, width(), height());
    }
}

int main(int argc, char** argv)
{
    Viewer window("Polygon Laplace Demo", 800, 600);
    if (argc == 2)
        window.load_mesh(argv[1]);
    return window.run();
}
