// Copyright 2011-2019 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include <pmp/visualization/MeshViewer.h>
#include "Surface/Smoothing.h"
#include <imgui.h>
#include "Surface/[AW11]Laplace.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/GeodesicsInHeat.h"
#include "Surface/Curvature.h"

enum InsertedPoint
{
    Centroid_ = 0,
    AreaMinimizer = 2,
};

enum LaplaceMethods {
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3
};
using namespace pmp;

class Viewer : public MeshViewer
{
public:
    Viewer(const char* title, int width, int height);

protected:
    void process_imgui() override;

private:
    Smoothing smoother_;
};

Viewer::Viewer(const char* title, int width, int height)
    : MeshViewer(title, width, height), smoother_(mesh_)
{
    crease_angle_ = 180.0;
}

void Viewer::process_imgui()
{
    MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();
    ImGui::PushItemWidth(100);
    ImGui::PushItemWidth(100);
    ImGui::SliderFloat("Hyperparameter Alexa & Wardetzky Laplace",
                       &poly_laplace_lambda_, 0.01, 3.0);
    ImGui::PushItemWidth(100);
    ImGui::SliderFloat("Hyperparameter deGoes Laplace",
                       &deGoes_laplace_lambda_, 0.01, 3.0);
    ImGui::PopItemWidth();

    ImGui::Spacing();

    ImGui::Text("Choose your Laplacian");

    ImGui::Spacing();

    static int laplace = 0;
    ImGui::RadioButton("Polysimple Laplace", &laplace, 0);
    ImGui::RadioButton("Alexa & Wardetzky Laplace", &laplace, 1);
    ImGui::RadioButton("Diamond", &laplace, 2);
    ImGui::RadioButton("deGoes Laplace", &laplace, 3);

    ImGui::Spacing();

    ImGui::Text("Choose your minimizing Point ");

    ImGui::Spacing();

    static int min_point = 2;
    ImGui::RadioButton("Centroid", &min_point, 0);
    ImGui::RadioButton("Area Minimizer", &min_point, 2);
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Smoothing", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static int iterations = 10;
        ImGui::PushItemWidth(100);
        ImGui::SliderInt("Iterations", &iterations, 1, 100);
        ImGui::PopItemWidth();

        static float timestep = 0.001;
        float lb =  0.001;
        float ub =  1.0;
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("TimeStep", &timestep, lb, ub);
        ImGui::PopItemWidth();

        if (ImGui::Button("Implicit Smoothing"))
        {
            // does the mesh have a boundary?
            bool has_boundary = false;
            for (auto v : mesh_.vertices())
                if (mesh_.is_boundary(v))
                    has_boundary = true;

            // only re-scale if we don't have a (fixed) boundary
            bool rescale = !has_boundary;

            Scalar dt =  timestep * radius_ * radius_;
            try
            {
                smoother_.implicit_smoothing(dt, laplace, min_point,rescale);
            }
            catch (const SolverException& e)
            {
                std::cerr << e.what() << std::endl;
                return;
            }
            update_mesh();
        }
    }
    if (ImGui::CollapsingHeader("Curvature", ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("Mean Curvature"))
        {
            Curvature analyzer(mesh_, false);
            analyzer.visualize_curvature(laplace, min_point, true);
            mesh_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
    if (ImGui::CollapsingHeader("Geodesic in Heat", ImGuiTreeNodeFlags_DefaultOpen)) {
        static int ts = 0;
        ImGui::Text("Choose your diffusion time step ");

        ImGui::RadioButton("Mean edge length", &ts, 0);
        ImGui::RadioButton("Max edge length", &ts, 1);
        ImGui::RadioButton("Max diagonal length", &ts, 2);

        if (ImGui::Button("Compute Distances Vertex 0")) {

            GeodesicsInHeat heat(mesh_, laplace, min_point,
                                 false, false, DiffusionStep(ts));
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance(0, dist, geodist);

            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");

            mesh_.use_checkerboard_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
}

int main(int argc, char** argv)
{
#ifndef __EMSCRIPTEN__
    Viewer window("Smoothing", 800, 600);
    if (argc == 2)
        window.load_mesh(argv[1]);
    return window.run();
#else
    Viewer window("Smoothing", 800, 600);
    window.load_mesh(argc == 2 ? argv[1] : "input.off");
    return window.run();
#endif
}
