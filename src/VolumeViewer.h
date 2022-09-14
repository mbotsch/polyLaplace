//=============================================================================
// Copyright 2021 Astrid Bunge, Mario Botsch, Marc Alexa.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <VolumeMeshViewer.h>
//=============================================================================

class VolumeViewer : public VolumeMeshViewer
{
public: // public methods
    VolumeViewer(const char* title, int width, int height)
        : VolumeMeshViewer(title, width, height)
    {
        set_draw_mode("Smooth Shading");
    }
    bool load_mesh(const char* filename) override;


protected:
    // overload GUI callbacks
    void keyboard(int key, int code, int action, int mod) override;
    void process_imgui() override;
    void mouse(int button, int action, int mods) override;

    int face_point_;
    int cell_point_;
    int laplace_matrix;
    std::string filename_;
    // private methods
};
