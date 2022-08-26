//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================
#pragma once
//=============================================================================

#include <pmp/visualization/GL.h>
#include <pmp/visualization/Shader.h>
#include <pmp/MatVec.h>
#include <pmp/Types.h>

#include <OpenVolumeMesh/Attribs/NormalAttrib.hh>

#include <cfloat>

#include "VolumeMesh.h"

//=============================================================================

#define OVM OpenVolumeMesh

using pmp::ivec2;
using pmp::ivec3;
using pmp::mat3;
using pmp::mat4;
using pmp::vec2;
using pmp::vec3;

//=============================================================================

class VolumeMeshGL : public VolumeMesh
{
public:
    //! constructor
    VolumeMeshGL();

    //! destructor
    ~VolumeMeshGL();

    //! get front color
    const vec3& front_color() const { return front_color_; }

    //! set front color
    void set_front_color(const vec3& color) { front_color_ = color; }

    //! get back color
    const vec3& back_color() const { return back_color_; }

    //! set back color
    void set_back_color(const vec3& color) { back_color_ = color; }

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

    //! get embedding value for reducing cell size
    float embedding() const { return embedding_; }

    // set embedding value for reducing cell size
    void set_embedding(float e) { embedding_ = e; }

    //! draw the mesh
    void draw(const mat4& projection_matrix, const mat4& modelview_matrix,
              const std::string& draw_mode);

    //! update all opengl buffers for efficient core profile rendering
    void update_opengl_buffers();

    //! load a gradient texture from blue to red
    void use_gradient_texture();

    //! load a checkerboard texture
    //! not very useful for Volume Mesh
    void use_checkerboard_texture();

    //! get the bounding sphere that encloses the mesh
    std::pair<vec3, float> get_bounding_sphere();

private: //-------------------triangulation of face----------------------------
    struct Triangulation
    {
        explicit Triangulation(pmp::Scalar a = FLT_MAX, int s = -1)
            : area(a), split(s)
        {
        }
        pmp::Scalar area;
        int split;
    };

    // table to hold triangulation data
    std::vector<Triangulation> triangulation_;

    // valence of currently triangulated polygon
    unsigned int polygon_valence_;

    // reserve n*n array for computing triangulation
    void init_triangulation(unsigned int n)
    {
        triangulation_.clear();
        triangulation_.resize(n * n);
        polygon_valence_ = n;
    }

    // access triangulation array
    Triangulation& triangulation(int start, int end)
    {
        return triangulation_[polygon_valence_ * start + end];
    }

    // compute squared area of triangle. used for triangulate().
    inline pmp::Scalar area(const vec3& p0, const vec3& p1,
                            const vec3& p2) const
    {
        return sqrnorm(cross(p1 - p0, p2 - p0));
    }

    // triangulate a polygon such that the sum of squared triangle areas is minimized.
    // this prevents overlapping/folding triangles for non-convex polygons.
    void triangulate(const std::vector<vec3>& points,
                     std::vector<ivec3>& triangles);

private:
    //! OpenGL buffers
    GLuint vertex_array_object_;
    GLuint vertex_buffer_;
    GLuint triangle_buffer_;
    GLuint embedding_buffer_;
    GLuint normal_buffer_;
    GLuint tex_coord_buffer_;
    GLuint edge_buffer_;
    GLuint feature_buffer_;

    //! OpenGl texture
    GLuint texture_;

    //! buffer sizes
    GLsizei n_vertices_;
    GLsizei n_edges_;
    GLsizei n_triangles_;

    //! shaders
    pmp::Shader phong_shader_;

    //! material properties
    vec3 front_color_, back_color_;
    float ambient_, diffuse_, specular_, shininess_, alpha_;

    //! vertex embedding
    float embedding_;

    //! normal attribute for all vertices
    OVM::NormalAttrib<
        OVM::GeometryKernel<OVM::Geometry::Vec3d, OVM::TopologyKernel>>
        normal_attrib_;

    //! texture coordinates for each vertex
    OpenVolumeMesh::VertexPropertyT<vec2> tex_coord_ =
        request_vertex_property<vec2>("texture_coordinates");
};