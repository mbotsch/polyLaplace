//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "../common_util.h"
#include <pmp/visualization/MeshViewer.h>
#include "Surface/Smoothing.h"
#include <imgui.h>
#include "Surface/[AW11]Laplace.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/GeodesicsInHeat.h"
#include "Surface/Curvature.h"
#include "Surface/Parameterization.h"
using namespace pmp;

class Viewer : public MeshViewer
{
public:
    Viewer(const char* title, int width, int height);

protected:
    void process_imgui() override;
    void draw(const std::string& _draw_mode) override;

private:
    Smoothing smoother_;
    bool show_uv_layout_;
};

Viewer::Viewer(const char* title, int width, int height)
    : MeshViewer(title, width, height), smoother_(mesh_)
{
    set_draw_mode("Hidden Line");
    crease_angle_ = 0.0;
    show_uv_layout_ = false;
}

void Viewer::process_imgui()
{
    // initialize settings
    static int laplace = 0;
    static int min_point = 2;

    MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Laplace", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::RadioButton("Alexa & Wardetzky Laplace", &laplace, 1);
        ImGui::RadioButton("deGoes Laplace", &laplace, 3);
        ImGui::RadioButton("Polysimple/-robust Laplace", &laplace, 0);
        ImGui::RadioButton("Diamond", &laplace, 2);
        ImGui::RadioButton("Harmonic", &laplace, 4);

        ImGui::Spacing();
        if(laplace == 0 || laplace == 2)
        {
            ImGui::Text("Choose your minimizing Point ");

            ImGui::Spacing();

            ImGui::RadioButton("Centroid", &min_point, 0);
            ImGui::RadioButton("Area Minimizer & Min Norm Weights", &min_point,
                               2);
            ImGui::RadioButton("Trace Minimizer & Discr. Harm. Weights",
                               &min_point, 3);

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


    if (ImGui::CollapsingHeader("Applications", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Indent(10);

        if (ImGui::Button("Mean Curvature"))
        {
            Curvature analyzer(mesh_, false);
            analyzer.visualize_curvature(laplace, min_point, true);
            mesh_.use_cold_warm_texture();
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

            double condition_number;
            heat.compute_geodesics(condition_number);
            heat.getDistance(0, dist, geodist);

            update_mesh();
            mesh_.use_checkerboard_texture();
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
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
            show_uv_layout_ = true;
        }

        ImGui::Unindent(10);
    }
}


void Viewer::draw(const std::string& draw_mode)
{
    // draw the mesh
    mesh_.draw(projection_matrix_, modelview_matrix_, draw_mode);

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
        mesh_.draw(P, M, "Texture Layout");

        // reset viewport
        glViewport(0, 0, width(), height());
    }
}


int main(int argc, char** argv)
{
#ifndef __EMSCRIPTEN__
    Viewer window("Laplace Demo", 800, 600);
    if (argc == 2)
        window.load_mesh(argv[1]);
    return window.run();
#else
    Viewer window("Smoothing", 800, 600);
    window.load_mesh(argc == 2 ? argv[1] : "input.off");
    return window.run();
#endif
}
