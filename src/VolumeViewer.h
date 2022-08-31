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
    virtual bool load_mesh(const char* filename) override;

    bool write_histogram(const char* filename);

protected:
    // overload GUI callbacks
    virtual void keyboard(int key, int code, int action, int mod) override;
    virtual void process_imgui() override;
    virtual void mouse(int button, int action, int mods) override;

    int face_point_;
    int cell_point_;
    int laplace_matrix;
    std::string filename_;
    // private methods
};
