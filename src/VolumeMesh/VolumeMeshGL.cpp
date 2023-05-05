//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMeshGL.h"
#include "PhongShader.h"

//=============================================================================

VolumeMeshGL::VolumeMeshGL() : normal_attrib_(*this)
{
    // initialize GL buffers to zero
    vertex_array_object_ = 0;
    vertex_buffer_ = 0;
    triangle_buffer_ = 0;
    embedding_buffer_ = 0;
    normal_buffer_ = 0;
    edge_buffer_ = 0;

    // initialize buffer sizes
    n_vertices_ = 0;
    n_edges_ = 0;
    n_triangles_ = 0;

    // material parameters
    front_color_ = vec3(0.5, 0.0, 0.0);
    back_color_ = vec3(0.6, 0.6, 0.6);
    ambient_ = 0.1;
    diffuse_ = 0.8;
    specular_ = 0.6;
    shininess_ = 100.0;
    alpha_ = 1.0;

    texture_ = 0;
}

//-----------------------------------------------------------------------------

VolumeMeshGL::~VolumeMeshGL()
{
    // delete OpenGL buffers
    glDeleteBuffers(1, &vertex_buffer_);
    glDeleteBuffers(1, &embedding_buffer_);
    glDeleteBuffers(1, &normal_buffer_);
    glDeleteBuffers(1, &tex_coord_buffer_);
    glDeleteBuffers(1, &edge_buffer_);
    glDeleteBuffers(1, &feature_buffer_);
    glDeleteVertexArrays(1, &vertex_array_object_);
}

//-----------------------------------------------------------------------------

void VolumeMeshGL::draw(const mat4& projection_matrix,
                        const mat4& modelview_matrix,
                        const std::string& draw_mode)
{
    if (!vertex_array_object_)
    {
        update_opengl_buffers();
    }

    if (!phong_shader_.is_valid())
    {
        if (!phong_shader_.source(phong_vshader, phong_fshader))
        {
            exit(1);
        }
    }

    if (n_vertices() == 0)
    {
        return;
    }

    glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);
    //glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT);

    mat4 mv_matrix = modelview_matrix;
    mat4 mvp_matrix = projection_matrix * modelview_matrix;
    mat3 n_matrix = inverse(transpose(linear_part(mv_matrix)));

    phong_shader_.use();
    phong_shader_.set_uniform("modelview_projection_matrix", mvp_matrix);
    phong_shader_.set_uniform("modelview_matrix", mv_matrix);
    phong_shader_.set_uniform("normal_matrix", n_matrix);
    phong_shader_.set_uniform("point_size", 5.0f);
    phong_shader_.set_uniform("light1", vec3(1.0, 1.0, 1.0));
    phong_shader_.set_uniform("light2", vec3(-1.0, 1.0, 1.0));
    phong_shader_.set_uniform("front_color", front_color_);
    phong_shader_.set_uniform("back_color", back_color_);
    phong_shader_.set_uniform("ambient", ambient_);
    phong_shader_.set_uniform("diffuse", diffuse_);
    phong_shader_.set_uniform("specular", specular_);
    phong_shader_.set_uniform("shininess", shininess_);
    phong_shader_.set_uniform("alpha", alpha_);
    phong_shader_.set_uniform("embedded", false);
    phong_shader_.set_uniform("embedding_size", embedding_);
    phong_shader_.set_uniform("use_lighting", true);
    phong_shader_.set_uniform("use_texture", false);

    glBindVertexArray(vertex_array_object_);

    if (draw_mode == "Points")
    {
        glEnable(GL_PROGRAM_POINT_SIZE);
        glDrawArrays(GL_POINTS, 0, n_vertices_);
    }

    else if (draw_mode == "Line")
    {
        glDepthFunc(GL_LEQUAL);
        glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);
        phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
        phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
        phong_shader_.set_uniform("use_lighting", false);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
        glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
        glDepthFunc(GL_LESS);
    }

    else if (draw_mode == "Hidden Line")
    {
        if (n_faces())
        {
            // draw faces
            glDepthRange(0.001, 1.0);
            glDrawArrays(GL_TRIANGLES, 0, n_vertices_);
            glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

            // overlay edges
            glDepthRange(0.0, 1.0);
            glDepthFunc(GL_LEQUAL);
            phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("use_lighting", false);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
            glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
            glDepthFunc(GL_LESS);
        }
    }

    else if (draw_mode == "Embedded")
    {
        if (n_faces())
        {
            // draw faces
            glDisable(GL_CULL_FACE);
            phong_shader_.set_uniform("embedded", true);
            glDrawArrays(GL_TRIANGLES, 0, n_vertices_);
            phong_shader_.set_uniform("embedded", false);
        }
    }

    else if (draw_mode == "Smooth Shading")
    {
        if (n_faces())
        {
            glDrawArrays(GL_TRIANGLES, 0, n_vertices_);
        }
    }

    else if (draw_mode == "Texture")
    {
        if (n_faces())
        {
            phong_shader_.set_uniform("front_color", vec3(0.3, 0.3, 0.3));
            phong_shader_.set_uniform("back_color", vec3(0.9, 0.9, 0.9));
            phong_shader_.set_uniform("use_texture", true);
            glBindTexture(GL_TEXTURE_2D, texture_);
            glDrawArrays(GL_TRIANGLES, 0, n_vertices_);
        }
    }

    glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);

    glBindVertexArray(0);
    glCheckError();
}

//-----------------------------------------------------------------------------

void VolumeMeshGL::update_opengl_buffers()
{
    if (!vertex_array_object_)
    {
        glGenVertexArrays(1, &vertex_array_object_);
        glBindVertexArray(vertex_array_object_);
        glGenBuffers(1, &vertex_buffer_);
        glGenBuffers(1, &triangle_buffer_);
        glGenBuffers(1, &embedding_buffer_);
        glGenBuffers(1, &normal_buffer_);
        glGenBuffers(1, &edge_buffer_);
        glGenBuffers(1, &tex_coord_buffer_);
    }

    glBindVertexArray(vertex_array_object_);

    auto vertex_indices = request_vertex_property<float>();

    std::vector<vec3> position_array;
    std::vector<vec3> normal_array;
    std::vector<vec3> embedding_array;
    std::vector<vec2> tex_array;
    std::vector<ivec3> triangles;

    OVM::Geometry::Vec3f v;
    OVM::Geometry::Vec3f normal;
    OVM::Geometry::Vec3f embedding;

    size_t vidx(0);

    if (n_faces())
    {
        position_array.reserve(3 * n_faces());
        normal_array.reserve(3 * n_faces());
        embedding_array.reserve(3 * n_faces());
        tex_array.reserve(3 * n_faces());

        std::vector<vec3> corner_positions;
        std::vector<vec3> corner_normals;
        std::vector<vec3> corner_embeddings;
        std::vector<vec2> corner_tex_coords;
        std::vector<OVM::VertexHandle> corner_vertices;

        normal_attrib_.update_face_normals();

        for (auto c_it = c_iter(); c_it.valid(); ++c_it)
        {
            auto c_bary = barycenter(*c_it);

            for (auto chf_it = chf_iter(*c_it); chf_it.valid(); ++chf_it)
            {
                corner_positions.clear();
                corner_normals.clear();
                corner_embeddings.clear();
                corner_tex_coords.clear();
                corner_vertices.clear();

                auto hf = (*chf_it);

                //auto f_vertices = halfface_vertices(*chf_it);
                normal = normal_attrib_[hf];

                // get the needed values from each vertex of the halfface
                for (auto hfv_it = hfv_iter(hf); hfv_it.valid(); ++hfv_it)
                {
                    v = vertex(*hfv_it);
                    embedding = c_bary - v;
                    corner_positions.emplace_back(vec3(v[0], v[1], v[2]));
                    corner_normals.emplace_back(
                        vec3(normal[0], normal[1], normal[2]));
                    corner_embeddings.emplace_back(
                        vec3(embedding[0], embedding[1], embedding[2]));
                    // this should be accessed from a vertex property
                    corner_tex_coords.emplace_back(tex_coord_[*hfv_it]);
                    corner_vertices.emplace_back(hfv_it->idx());
                }
                assert(corner_vertices.size() >= 3);

                // creates triangulation of face and returns the indices of the
                // triangle points
                triangulate(corner_positions, triangles);

                // for each triangle save the corresponding values
                for (auto& t : triangles)
                {
                    int i0 = t[0];
                    int i1 = t[1];
                    int i2 = t[2];

                    position_array.push_back(corner_positions[i0]);
                    position_array.push_back(corner_positions[i1]);
                    position_array.push_back(corner_positions[i2]);

                    normal_array.push_back(corner_normals[i0]);
                    normal_array.push_back(corner_normals[i1]);
                    normal_array.push_back(corner_normals[i2]);

                    embedding_array.push_back(corner_embeddings[i0]);
                    embedding_array.push_back(corner_embeddings[i1]);
                    embedding_array.push_back(corner_embeddings[i2]);

                    tex_array.push_back(corner_tex_coords[i0]);
                    tex_array.push_back(corner_tex_coords[i1]);
                    tex_array.push_back(corner_tex_coords[i2]);

                    vertex_indices[corner_vertices[i0]] = vidx++;
                    vertex_indices[corner_vertices[i1]] = vidx++;
                    vertex_indices[corner_vertices[i2]] = vidx++;
                }
            }
        }
    }
    else if (n_vertices())
    {
        position_array.reserve(n_vertices());
        normal_attrib_.update_vertex_normals();

        for (auto vert_iter = vertices_begin(); vert_iter != vertices_end();
             ++vert_iter)
        {
            v = vertex(*vert_iter);
            normal = normal_attrib_[*vert_iter];

            position_array.emplace_back(vec3(v[0], v[1], v[2]));
            normal_array.emplace_back(vec3(normal[0], normal[1], normal[2]));
        }
    }

    // write all arrays to the OpenGL buffer

    if (!position_array.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_);
        glBufferData(GL_ARRAY_BUFFER, position_array.size() * 3 * sizeof(float),
                     position_array.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(0);
        n_vertices_ = position_array.size();
    }
    else
    {
        n_vertices_ = 0;
    }

    if (!normal_array.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, normal_buffer_);
        glBufferData(GL_ARRAY_BUFFER, normal_array.size() * 3 * sizeof(float),
                     normal_array.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(1);
    }

    if (!embedding_array.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, embedding_buffer_);
        glBufferData(GL_ARRAY_BUFFER,
                     embedding_array.size() * 3 * sizeof(float),
                     embedding_array.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(2);
    }

    if (!tex_array.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, tex_coord_buffer_);
        glBufferData(GL_ARRAY_BUFFER, tex_array.size() * 2 * sizeof(float),
                     tex_array.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(3);
    }

    if (n_edges())
    {
        std::vector<unsigned int> edge_array;
        edge_array.reserve(n_edges());
        for (auto edge_iter = edges_begin(); edge_iter != edges_end();
             ++edge_iter)
        {
            edge_array.emplace_back(
                vertex_indices[edge_vertices(*edge_iter)[0]]);
            edge_array.emplace_back(
                vertex_indices[edge_vertices(*edge_iter)[1]]);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     edge_array.size() * sizeof(unsigned int),
                     edge_array.data(), GL_STATIC_DRAW);
        n_edges_ = edge_array.size();
    }
    else
    {
        n_edges_ = 0;
    }

    glBindVertexArray(0);
}

//-----------------------------------------------------------------------------

void VolumeMeshGL::use_gradient_texture()
{
    // delete old texture
    glDeleteTextures(1, &texture_);

    // generate checkerboard-like image
    const unsigned int res = 512;
    auto* tex = new GLubyte[res * res * 3];
    GLubyte* tp = tex;
    for (unsigned int x = 0; x < res; ++x)
    {
        for (unsigned int y = 0; y < res; ++y)
        {
            *(tp++) = std::floor(255 * ((double)x / (double)res));
            *(tp++) = 0;
            *(tp++) = std::floor(255 * (1.0 - (double)x / (double)res));
        }
    }

    // generate texture
    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, res, res, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // clean up
    delete[] tex;
}

//-----------------------------------------------------------------------------

void VolumeMeshGL::use_checkerboard_texture()
{
    // delete old texture
    glDeleteTextures(1, &texture_);

    // generate checkerboard-like image
    const unsigned int res = 512;
    auto* tex = new GLubyte[res * res * 3];
    GLubyte* tp = tex;
    for (unsigned int x = 0; x < res; ++x)
    {
        for (unsigned int y = 0; y < res; ++y)
        {
            if (((x & 0x10) == 0) ^ ((y & 0x10) == 0))
            {
                *(tp++) = 42;
                *(tp++) = 157;
                *(tp++) = 223;
            }
            else
            {
                *(tp++) = 255;
                *(tp++) = 255;
                *(tp++) = 255;
            }
        }
    }

    // generate texture
    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, res, res, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    // clean up
    delete[] tex;
}

//-----------------------------------------------------------------------------

std::pair<vec3, float> VolumeMeshGL::get_bounding_sphere()
{
    auto bounding_sphere = VolumeMesh::get_bounding_sphere();

    OVM::Vec3d center = bounding_sphere.first;
    float radius = bounding_sphere.second;

    return {vec3(center[0], center[1], center[2]), radius};
}

//-----------------------------------------------------------------------------

void VolumeMeshGL::triangulate(const std::vector<vec3>& points,
                               std::vector<ivec3>& triangles)
{
    const int n = points.size();

    triangles.clear();
    triangles.reserve(n - 2);

    // triangle? nothing to do
    if (n == 3)
    {
        triangles.emplace_back(0, 1, 2);
        return;
    }

    // quad? simply compare to two options
    else if (n == 4)
    {
        if (area(points[0], points[1], points[2]) +
                area(points[0], points[2], points[3]) <
            area(points[0], points[1], points[3]) +
                area(points[1], points[2], points[3]))
        {
            triangles.emplace_back(0, 1, 2);
            triangles.emplace_back(0, 2, 3);
        }
        else
        {
            triangles.emplace_back(0, 1, 3);
            triangles.emplace_back(1, 2, 3);
        }
        return;
    }

    // n-gon with n>4? compute triangulation by dynamic programming
    init_triangulation(n);
    int i, j, m, k, imin;
    pmp::Scalar w, wmin;

    // initialize 2-gons
    for (i = 0; i < n - 1; ++i)
    {
        triangulation(i, i + 1) = Triangulation(0.0, -1);
    }

    // n-gons with n>2
    for (j = 2; j < n; ++j)
    {
        // for all n-gons [i,i+j]
        for (i = 0; i < n - j; ++i)
        {
            k = i + j;

            wmin = FLT_MAX;
            imin = -1;

            // find best split i < m < i+j
            for (m = i + 1; m < k; ++m)
            {
                w = triangulation(i, m).area +
                    area(points[i], points[m], points[k]) +
                    triangulation(m, k).area;

                if (w < wmin)
                {
                    wmin = w;
                    imin = m;
                }
            }

            triangulation(i, k) = Triangulation(wmin, imin);
        }
    }

    // build triangles from triangulation table
    std::vector<ivec2> todo;
    todo.reserve(n);
    todo.emplace_back(0, n - 1);
    while (!todo.empty())
    {
        ivec2 tri = todo.back();
        todo.pop_back();
        int start = tri[0];
        int end = tri[1];
        if (end - start < 2)
            continue;
        int split = triangulation(start, end).split;

        triangles.emplace_back(start, split, end);

        todo.emplace_back(start, split);
        todo.emplace_back(split, end);
    }
}