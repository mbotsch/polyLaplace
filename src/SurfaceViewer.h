//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/mesh_viewer.h>
#include "Surface/SpectralProcessing.h"
#include "Surface/GeodesicsInHeat.h"

//=============================================================================

class Viewer : public pmp::MeshViewer
{
public: // public methods
    Viewer(const char* title, int width, int height)
        : MeshViewer(title, width, height),
          compare_sphere(false),
          compare_cube(false)
    {
        set_draw_mode("Hidden Line");
    }
    void load_mesh(const char* filename) override;

protected:
    // overload GUI callbacks
    void keyboard(int key, int code, int action, int mod) override;
    void process_imgui() override;
    void update_mesh() override;
    void draw(const std::string& _draw_mode) override;

    void mouse(int button, int action, int mods) override;

    void dualize();
    void close_holes();
    void Centroid();
    void insert_points(unsigned int min_point);
    static void calc_colors(int minpoint, SurfaceMesh& mesh);

private: // private methods
    int laplace_matrix, min_point_;
    double cond_maxi = -1, cond_mini = -1;
    bool compare_sphere;
    bool compare_cube;
    DiffusionStep time_step_;
    SurfaceMesh original = mesh_;
};

//=============================================================================
