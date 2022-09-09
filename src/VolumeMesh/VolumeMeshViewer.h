//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include "VolumeMeshGL.h"
#include <pmp/visualization/TrackballViewer.h>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

//=============================================================================

class VolumeMeshViewer : public pmp::TrackballViewer
{
public:
    //! constructor
    VolumeMeshViewer(const char* title, int width, int height,
                     bool showgui = true);

    //! destructor
    ~VolumeMeshViewer();

    //! load a mesh from file \c filename
    virtual bool load_mesh(const char* filename);

    //! update mesh normals and all buffers for OpenGL rendering.  call this
    //! function whenever you change either the vertex positions or the
    //! triangulation of the mesh
    virtual void update_mesh();

    //! draw the scene in different draw modes
    virtual void draw(const std::string& draw_mode) override;

    //! handle ImGUI interface
    virtual void process_imgui() override;

    //! this function handles keyboard events
    virtual void keyboard(int key, int code, int action, int mod) override;

    //! this function hadle mouse events
    virtual void mouse(int button, int action, int mods) override;

    //! draw the scene
    virtual void display() override;

    //! pick vertex next to the the mouse click
    VHandle pick_vertex(double x, double y);

protected:
    VolumeMeshGL mesh_;    //!< the volume mesh
    std::string filename_; //!< the current file that stores the mesh
    OpenVolumeMesh::IO::FileManager
        file_manager_; //!< filemanager for reading mesh

    float embedding_;            //!< value for the embedding of vertices
    float near_clipping_ = 0.0f; //!< near clipping plane for openGL
    float alpha_;                //!< alpha value for opacity
};
