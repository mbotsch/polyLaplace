//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/MeshViewer.h>
#include "Surface/Smoothing.h"
#include "Surface/SpectralProcessing.h"

//=============================================================================

class Viewer : public pmp::MeshViewer
{
public: // public methods
    Viewer(const char* title, int width, int height)
        : MeshViewer(title, width, height),
          smooth_(mesh_),
          compare_sphere(false),
          compare_cube(false)
    {
        set_draw_mode("Hidden Line");
    }
    virtual void load_mesh(const char* filename) override;
protected:
    // overload GUI callbacks
    virtual void keyboard(int key, int code, int action, int mod) override;
    virtual void process_imgui() override;
    virtual void update_mesh() override;
    virtual void draw(const std::string& _draw_mode) override;

    void mouse(int button, int action, int mods) override;

    void dualize();
    void close_holes();
    void open_holes();
    void Centroid();
    void insert_points(unsigned int min_point);
    void mirror_mesh();

private: // private methods
    Smoothing smooth_;

    unsigned int laplace_matrix, min_point_, degree_;
    CoarseDimension coarseningType_;
    bool compare_sphere;
    bool compare_cube;
    bool time_step_;
    Eigen::MatrixXd basis_values;
    SurfaceMeshGL original = mesh_;

    std::vector<Face> holes_;
};

//=============================================================================
