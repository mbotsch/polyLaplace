//=============================================================================
// Copyright (C) 2013-2019 The pmp-library developers
//
// This file is part of the Polygon Mesh Processing Library.
// Distributed under the terms of the MIT license, see LICENSE.txt for details.
//
// SPDX-License-Identifier: MIT
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/TrackballViewer.h>
#include "MySurfaceMeshGL.h"

//=============================================================================

namespace pmp {

    class MyViewer : public pmp::TrackballViewer {
    public:
        //! constructor
        MyViewer(const char *title, int width, int height, bool showgui = true);

        //! destructor
        virtual ~MyViewer();

        //! load a mesh from file \c filename
        virtual bool load_mesh(const char *filename);

        //! load a matcap texture from file \c filename
        bool load_matcap(const char *filename);

        //! load a texture from file \c filename
        bool load_texture(const char *filename, GLint format = GL_RGB,
                          GLint min_filter = GL_LINEAR_MIPMAP_LINEAR,
                          GLint mag_filter = GL_LINEAR,
                          GLint wrap = GL_CLAMP_TO_EDGE);

        //! update mesh normals and all buffers for OpenGL rendering.  call this
        //! function whenever you change either the vertex positions or the
        //! triangulation of the mesh
        virtual void update_mesh();

        //! draw the scene in different draw modes
        virtual void draw(const std::string &draw_mode) override;

        //! handle ImGUI interface
        virtual void process_imgui() override;

        //! this function handles keyboard events
        virtual void keyboard(int key, int code, int action, int mod) override;

        //! get vertex closest to 3D position Distributed under the mouse cursor
        pmp::Vertex pick_vertex(int x, int y);

    protected:
        pmp::MySurfaceMeshGL mesh_;   //!< the mesh
        std::string filename_; //!< the current file
        float crease_angle_;


    private:
        /// which energy to use?
        Energy energy_ = CROSS_EDGE_DIFF;
        TextureMode textureMode_  = ColdWarmTexture;
    };
}
//=============================================================================
