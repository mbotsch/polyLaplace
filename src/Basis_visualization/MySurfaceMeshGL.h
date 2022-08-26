// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#pragma once

#include <limits>

#include <pmp/SurfaceMesh.h>
#include <pmp/visualization/GL.h>
#include <pmp/visualization/Shader.h>
#include <pmp/MatVec.h>

#include "MyShader.h"

namespace pmp {
enum Energy
{
    LINEAR_UNIFORM,
    CROSS_EDGE_DIFF
};
enum TextureMode
{
    ColdWarmTexture,
    CheckerboardTexture,
    MatCapTexture,
    OtherTexture,
    CoffeeTexture
};
//! Class for rendering surface meshes using OpenGL
//! \ingroup visualization
class MySurfaceMeshGL : public SurfaceMesh
{
public:
    //! Constructor
    MySurfaceMeshGL();

    //! default destructor
    ~MySurfaceMeshGL();

    //! get front color
    const vec3 &front_color() const { return front_color_; }

    //! set front color
    void set_front_color(const vec3 &color) { front_color_ = color; }

    //! get back color
    const vec3 &back_color() const { return back_color_; }

    //! set back color
    void set_back_color(const vec3 &color) { back_color_ = color; }

    //! get ambient reflection coefficient
    float ambient() const { return ambient_; }

    //! set ambient reflection coefficient
    void set_ambient(float a) { ambient_ = a; }

    //! get diffuse reflection coefficient
    float diffuse() const { return diffuse_; }

    //! set diffuse reflection coefficient
    void set_diffuse(float d) { diffuse_ = d; }

    //! get specular reflection coefficient
    float specular() const { return specular_; }

    //! set specular reflection coefficient
    void set_specular(float s) { specular_ = s; }

    //! get specular shininess coefficient
    float shininess() const { return shininess_; }

    //! set specular shininess coefficient
    void set_shininess(float s) { shininess_ = s; }

    //! get alpha value for transparent rendering
    float alpha() const { return alpha_; }

    //! set alpha value for transparent rendering
    void set_alpha(float a) { alpha_ = a; }

    //! get crease angle (in degrees) for visualization of sharp edges
    Scalar crease_angle() const { return crease_angle_; }

    //! set crease angle (in degrees) for visualization of sharp edges
    void set_crease_angle(Scalar ca);

    //! get basis function for visualization of single basis functions
    int basis_function() const { return basis_function_; }

    //! set basis function for visualization of single basis functions
    void set_basis_function(int bf) { basis_function_ = bf; };

    //! get basis function curvature for visualization of single basis functions
    float basis_elevation() const { return basis_elevation_; }

    //! set basis function curvature for visualization of single basis functions
    void set_basis_elevation(float be) { basis_elevation_ = be; };

    //! draw the mesh
    void draw(const mat4 &projection_matrix, const mat4 &modelview_matrix,
              const std::string draw_mode);

    //! update all opengl buffers for efficient core profile rendering
    void update_opengl_buffers();

    void update_weights(Energy energy);

    int edge_node_to_idx(pmp::Face &f, pmp::Halfedge &he);

    int vertex_node_to_idx(pmp::Face &f, pmp::Halfedge &he);
    //! use color map to visualize scalar fields
    void use_cold_warm_texture();

    //! setup checkerboard texture
    void use_checkerboard_texture();

    void use_coffee_checkerboard_texture();

    void set_energy(Energy e) { energy_ = e; };

    //! load texture from file
    //! \param filename the location and name of the texture
    //! \param format internal format (GL_RGB, GL_RGBA, GL_SRGB8, etc.)
    //! \param min_filter interpolation filter for minification
    //! \param mag_filter interpolation filter for magnification
    //! \param wrap texture coordinates wrap preference
    bool load_texture(const char *filename, GLint format = GL_RGB,
                      GLint min_filter = GL_LINEAR_MIPMAP_LINEAR,
                      GLint mag_filter = GL_LINEAR,
                      GLint wrap = GL_CLAMP_TO_EDGE);

    //! Load mat-cap texture from file. The mat-cap will be used
    //! whenever the drawing mode is "Texture". This also means
    //! that you cannot have texture and mat-cap at the same time.
    //! \param filename the location and name of the texture
    //! \sa See src/apps/mview.cpp for an example usage.
    bool load_matcap(const char *filename);

private:
    //! OpenGL buffers
    GLuint vertex_array_object_;
    GLuint vertex_buffer_;
    GLuint normal_buffer_;
    GLuint tex_coord_buffer_;
    GLuint coeff_buffer_;
    GLuint edge_buffer_;
    GLuint feature_buffer_;

    //! buffer sizes
    GLsizei n_vertices_;
    GLsizei n_edges_;
    GLsizei n_triangles_;
    GLsizei n_features_;
    bool have_texcoords_;

    //! shaders
    MyShader phong_shader_;
    MyShader phong_tess_shader_;
    MyShader matcap_shader_;

    //! material properties
    vec3 front_color_, back_color_;
    float ambient_, diffuse_, specular_, shininess_, alpha_;
    bool srgb_;
    float crease_angle_;

    //! 1D texture for scalar field rendering
    GLuint texture_;
    TextureMode texture_mode_;
    Energy energy_;
    int basis_function_;
    float basis_elevation_;
};

} // namespace pmp
