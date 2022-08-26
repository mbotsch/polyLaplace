// Copyright 2011-2020 the Polygon Mesh Processing Library developers.
// Distributed under a MIT-style license, see LICENSE.txt for details.

#include "MySurfaceMeshGL.h"
#include <pmp/visualization/PhongShader.h>
#include <pmp/visualization/MatCapShader.h>
#include <pmp/visualization/ColdWarmTexture.h>
#include <pmp/algorithms/SurfaceNormals.h>
#include <Misha_old/PolygonBasis.h>
#include "PhongTessShader.h"
#include "../../external/pmp-library/external/stb_image/stb_image.h"

namespace pmp {

MySurfaceMeshGL::MySurfaceMeshGL()
{
    // initialize GL buffers to zero
    vertex_array_object_ = 0;
    vertex_buffer_ = 0;
    normal_buffer_ = 0;
    tex_coord_buffer_ = 0;
    coeff_buffer_ = 0;
    edge_buffer_ = 0;
    feature_buffer_ = 0;

    // initialize buffer sizes
    n_vertices_ = 0;
    n_edges_ = 0;
    n_triangles_ = 0;
    n_features_ = 0;
    have_texcoords_ = false;

    // material parameters
    front_color_ = vec3(0.6, 0.6, 0.6);
    back_color_ = vec3(0.5, 0.0, 0.0);
    ambient_ = 0.1;
    diffuse_ = 0.8;
    specular_ = 0.6;
    shininess_ = 100.0;
    alpha_ = 1.0;
    srgb_ = false;
    crease_angle_ = 180.0;

    // initialize texture
    texture_ = 0;
    texture_mode_ = OtherTexture;
}

MySurfaceMeshGL::~MySurfaceMeshGL()
{
    // delete OpenGL buffers
    glDeleteBuffers(1, &vertex_buffer_);
    glDeleteBuffers(1, &normal_buffer_);
    glDeleteBuffers(1, &tex_coord_buffer_);
    glDeleteBuffers(1, &coeff_buffer_);
    glDeleteBuffers(1, &edge_buffer_);
    glDeleteBuffers(1, &feature_buffer_);
    glDeleteVertexArrays(1, &vertex_array_object_);
    glDeleteTextures(1, &texture_);
}

bool MySurfaceMeshGL::load_texture(const char *filename, GLint format,
                                   GLint min_filter, GLint mag_filter,
                                   GLint wrap)
{
#ifdef __EMSCRIPTEN__
    // emscripen/WebGL does not like mapmapping for SRGB textures
    if ((min_filter == GL_NEAREST_MIPMAP_NEAREST ||
         min_filter == GL_NEAREST_MIPMAP_LINEAR ||
         min_filter == GL_LINEAR_MIPMAP_NEAREST ||
         min_filter == GL_LINEAR_MIPMAP_LINEAR) &&
        (format == GL_SRGB8))
        min_filter = GL_LINEAR;
#endif

    // choose number of components (RGB or RGBA) based on format
    int loadComponents;
    GLint loadFormat;
    switch (format)
    {
        case GL_RGB:
        case GL_SRGB8:
            loadComponents = 3;
            loadFormat = GL_RGB;
            break;

        case GL_RGBA:
        case GL_SRGB8_ALPHA8:
            loadComponents = 4;
            loadFormat = GL_RGBA;
            break;

        default:
            loadComponents = 3;
            loadFormat = GL_RGB;
    }

    // load with stb_image
    int width, height, n;
    stbi_set_flip_vertically_on_load(true);
    unsigned char *img =
        stbi_load(filename, &width, &height, &n, loadComponents);
    if (!img)
        return false;

    // delete old texture
    glDeleteTextures(1, &texture_);

    // setup new texture
    glGenTextures(1, &texture_);
    glBindTexture(GL_TEXTURE_2D, texture_);

    // upload texture data
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, loadFormat,
                 GL_UNSIGNED_BYTE, img);

    // compute mipmaps
    if (min_filter == GL_LINEAR_MIPMAP_LINEAR)
    {
        glGenerateMipmap(GL_TEXTURE_2D);
    }

    // set texture parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, min_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, mag_filter);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, wrap);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, wrap);

    // use SRGB rendering?
    srgb_ = (format == GL_SRGB8);

    // free memory
    stbi_image_free(img);

    texture_mode_ = OtherTexture;
    return true;
}

bool MySurfaceMeshGL::load_matcap(const char *filename)
{
    if (!load_texture(filename, GL_RGBA, GL_LINEAR, GL_LINEAR,
                      GL_CLAMP_TO_EDGE))
        return false;

    texture_mode_ = MatCapTexture;
    return true;
}

void MySurfaceMeshGL::use_cold_warm_texture()
{
    if (texture_mode_ != ColdWarmTexture)
    {
        // delete old texture
        glDeleteTextures(1, &texture_);

        glGenTextures(1, &texture_);
        glBindTexture(GL_TEXTURE_2D, texture_);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 1, 0, GL_RGB,
                     GL_UNSIGNED_BYTE, cold_warm_texture);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_R, GL_MIRRORED_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_MIRRORED_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_MIRRORED_REPEAT);
        srgb_ = false;
        texture_mode_ = ColdWarmTexture;
    }
}

void MySurfaceMeshGL::use_checkerboard_texture()
{
    if (texture_mode_ != CheckerboardTexture)
    {
        // delete old texture
        glDeleteTextures(1, &texture_);

        // generate checkerboard-like image
        const unsigned int res = 1024;
        auto *tex = new GLubyte[res * res * 3];
        GLubyte *tp = tex;
        for (unsigned int x = 0; x < res; ++x)
        {
            for (unsigned int y = 0; y < res; ++y)
            {
                                if (y>512)
                                {
                                    *(tp++) = 255-(y/8);
                                    *(tp++) = 255-(y/8);
                                    *(tp++) = 255;
                                }
                                else
                                {
                                    *(tp++) = 255;
                                    *(tp++) = 255-(1024-y)/8;
                                    *(tp++) = 255-(1024-y)/8;
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

        srgb_ = false;
        texture_mode_ = CheckerboardTexture;
    }
}
void MySurfaceMeshGL::use_coffee_checkerboard_texture()
{

    if (texture_mode_ != CoffeeTexture)
    {
        // delete old texture
        glDeleteTextures(1, &texture_);

        // generate checkerboard-like image
        const unsigned int res = 2048;
        auto *tex = new GLubyte[res * res * 3];
        GLubyte *tp = tex;
        for (unsigned int x = 0; x < res; ++x)
        {
            for (unsigned int y = 0; y < res; ++y)
            {
                if (((x & 0x20) == 0) ^ ((y & 0x20) == 0))
                {
                    *(tp++) = 132;
                    *(tp++) = 172;
                    *(tp++) = 245;
                }
                else
                {
                    *(tp++) = 255 - (int)(y / 16.0);
                    *(tp++) = 255 - (int)(y /32.0);
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

        srgb_ = false;
        texture_mode_ = CoffeeTexture;
    }
}
void MySurfaceMeshGL::set_crease_angle(Scalar ca)
{
    if (ca != crease_angle_)
    {
        crease_angle_ = std::max(Scalar(0), std::min(Scalar(180), ca));
        update_opengl_buffers();
    }
}

int MySurfaceMeshGL::edge_node_to_idx(pmp::Face &f, pmp::Halfedge &he)
{
    int i = 0;
    pmp::Vertex v = from_vertex(he);
    for (auto hhe : halfedges(f))
    {
        pmp::Vertex vv = from_vertex(hhe);
        if ((i == 2 * valence(f) - 4) && (v.idx() == vv.idx()))
        {
            return 2 * valence(f) - 2;
        }
        else if (i == 2 * valence(f) - 2)
        {
            return 2 * valence(f) - 3;
        }
        else if (v.idx() == vv.idx())
        {
            return i + 1;
        }
        i += 2;
    }
}

int MySurfaceMeshGL::vertex_node_to_idx(pmp::Face &f, pmp::Halfedge &he)
{
    int i = 0;
    pmp::Vertex v = from_vertex(he);
    for (auto hhe : halfedges(f))
    {
        pmp::Vertex vv = from_vertex(hhe);
        if ((i == 2 * valence(f) - 2) && (v.idx() == vv.idx()))
        {
            return 2 * valence(f) - 1;
        }
        else if (v.idx() == vv.idx())
        {
            return i;
        }
        i += 2;
    }
}

void MySurfaceMeshGL::update_weights(Energy energy)
{
    if (energy != energy_)
    {
        if (energy == LINEAR_UNIFORM)
        {
            auto vweight = get_halfedge_property<Scalar>("v:weight");
            if (!vweight)
            {
                vweight = add_halfedge_property<Scalar>("v:weight");
                for (auto f : faces())
                {
                    int n = valence(f);
                    for (auto h : halfedges(f))
                    {
                        vweight[h] = (Scalar)1.0 / (Scalar)n;
                    }
                }
            }
        }
        else if (energy == CROSS_EDGE_DIFF)
        {
        }
    }
}

void MySurfaceMeshGL::update_opengl_buffers()
{
    // are buffers already initialized?
    if (!vertex_array_object_)
    {
        glGenVertexArrays(1, &vertex_array_object_);
        glBindVertexArray(vertex_array_object_);
        glGenBuffers(1, &vertex_buffer_);
        glGenBuffers(1, &normal_buffer_);
        glGenBuffers(1, &tex_coord_buffer_);
        glGenBuffers(1, &coeff_buffer_);
        glGenBuffers(1, &edge_buffer_);
        glGenBuffers(1, &feature_buffer_);
    }

    // activate VAO
    glBindVertexArray(vertex_array_object_);

    // get vertex properties
    auto vpos = get_vertex_property<Point>("v:point");
    auto vtex = get_vertex_property<TexCoord>("v:tex");
    auto htex = get_halfedge_property<TexCoord>("h:tex");
    auto vweight = get_halfedge_property<Scalar>("v:weight");

    if (!vweight)
    {
        //remove if computed by astrid
        vweight = add_halfedge_property<Scalar>("v:weight");
        for (auto f : faces())
        {
            int n = valence(f);
            for (auto h : halfedges(f))
            {
                vweight[h] = (Scalar)1.0 / (Scalar)n;
            }
        }
    }

    // index array for remapping vertex indices during duplication
    auto vertex_indices = add_vertex_property<size_t>("v:index");

    // produce arrays of points, normals, and texcoords
    // (duplicate vertices to allow for flat shading)
    std::vector<vec3> positionArray;
    std::vector<vec3> normalArray;
    std::vector<vec2> texArray;
    std::vector<float> coeffArray;

    // we have a mesh: fill arrays by looping over faces
    if (n_faces())
    {
        // reserve memory
        positionArray.reserve(3 * n_faces());
        normalArray.reserve(3 * n_faces());
        if (htex || vtex)
            texArray.reserve(3 * n_faces());
        coeffArray.reserve(3 * n_faces());

        // precompute normals for easy cases
        FaceProperty<Normal> fnormals;
        VertexProperty<Normal> vnormals;
        if (crease_angle_ < 1)
        {
            fnormals = add_face_property<Normal>("gl:fnormal");
            for (auto f : faces())
                fnormals[f] = SurfaceNormals::compute_face_normal(*this, f);
        }
        else if (crease_angle_ > 170)
        {
            vnormals = add_vertex_property<Normal>("gl:vnormal");
            for (auto v : vertices())
                vnormals[v] = SurfaceNormals::compute_vertex_normal(*this, v);
        }

        // convert from degrees to radians
        const Scalar creaseAngle = crease_angle_ / 180.0 * M_PI;

        size_t vidx(0);
        auto tmp_normal = add_vertex_property<vec3>("v:tmp_normal");
        auto tmp_tex = add_vertex_property<vec2>("v:tmp_tex");

        // loop over all faces
        for (auto f : faces())
        {
            Point c_pos(0.0);
            Normal c_normal(0.0);
            TexCoord c_tex(0.0);

            Vertex v;
            Normal n;
            TexCoord tex;

            std::vector<::Point<double, 3>> vertices;
            vertices.reserve(valence(f));

            for (auto h : halfedges(f))
            {
                v = from_vertex(h);

                c_pos += vweight[h] * vpos[v];

                if (fnormals)
                {
                    n = fnormals[f];
                }
                else if (vnormals)
                {
                    n = vnormals[v];
                }
                else
                {
                    n = SurfaceNormals::compute_corner_normal(*this, h,
                                                              creaseAngle);
                }
                c_normal += vweight[h] * n;
                tmp_normal[v] = (vec3)n;

                if (htex)
                {
                    tex = (dvec2)htex[h];
                }
                else if (vtex)
                {
                    tex = (dvec2)vtex[v];
                }
                else
                {
                    tex = dvec2(0.0, 0.0);
                }
                c_tex += vweight[h] * tex;
                tmp_tex[v] = (vec2)tex;
                vertices.push_back(
                    ::Point<double, 3>(vpos[v][0], vpos[v][1], vpos[v][2]));
            }

            EmbeddedRefinedPolygon<3> erPolygon;
            Eigen::MatrixXd P;
            erPolygon.center = ::Point<double, 3>(c_pos[0], c_pos[1], c_pos[2]);
            erPolygon.vertices = vertices;
            if (energy_ == CROSS_EDGE_DIFF)
            {
                typename PolygonElements<3, 2>::PoUProlongationSystem
                    pouProlongationSystem(erPolygon, 4);
                P = pouProlongationSystem.prolongation();
            }
            else if (energy_ == LINEAR_UNIFORM)
            {
                typename PolygonElements<3, 2>::PoUProlongationSystem
                    pouProlongationSystem(erPolygon, 1);
                P = pouProlongationSystem.prolongation();
            }

            Vertex from_v;
            Vertex to_v;

            int val = valence(f);
            int v_idx = 0;

            //TODO update center vertex and edge vertex based on the prolongation matrix

            for (auto h : halfedges(f))
            {
                from_v = from_vertex(h);
                to_v = to_vertex(h);

                positionArray.push_back((vec3)vpos[from_v]);
                positionArray.push_back((vec3)vpos[to_v]);
                positionArray.push_back((vec3)c_pos);
                positionArray.push_back(
                    (vec3)(0.5 * vpos[from_v] + 0.5 * vpos[to_v]));
                positionArray.push_back((vec3)(0.5 * vpos[to_v] + 0.5 * c_pos));
                positionArray.push_back(
                    (vec3)(0.5 * c_pos + 0.5 * vpos[from_v]));

                normalArray.push_back(tmp_normal[from_v]);
                normalArray.push_back(tmp_normal[to_v]);
                normalArray.push_back((vec3)c_normal);
                normalArray.push_back(0.5 * tmp_normal[from_v] +
                                      0.5 * tmp_normal[to_v]);
                normalArray.push_back(0.5 * tmp_normal[to_v] +
                                      0.5 * (vec3)c_normal);
                normalArray.push_back(0.5 * (vec3)c_normal +
                                      0.5 * tmp_normal[from_v]);

                if (htex || vtex || true)
                {
                    texArray.push_back(tmp_tex[from_v]);
                    texArray.push_back(tmp_tex[to_v]);
                    texArray.push_back((vec2)c_tex);
                    texArray.push_back(0.5 * tmp_tex[from_v] +
                                       0.5 * tmp_tex[to_v]);
                    texArray.push_back(0.5 * tmp_tex[to_v] + 0.5 * (vec2)c_tex);
                    texArray.push_back(0.5 * (vec2)c_tex +
                                       0.5 * tmp_tex[from_v]);
                }

                uint bf = basis_function_ % (2 * val);

                int start_idx = vertex_node_to_idx(f, h);
                int e_idx = edge_node_to_idx(f, h);
                auto nh = next_halfedge(h);
                int vv_idx = vertex_node_to_idx(f, nh);
                coeffArray.push_back(P.coeffRef(start_idx, bf));
                coeffArray.push_back(P.coeffRef(vv_idx, bf));
                coeffArray.push_back(P.coeffRef(3 * val, bf));
                coeffArray.push_back(P.coeffRef(e_idx, bf));
                coeffArray.push_back(
                    P.coeffRef(2 * val + (v_idx + 1) % val, bf));
                coeffArray.push_back(P.coeffRef(2 * val + v_idx, bf));
                ++v_idx;

                vertex_indices[from_vertex(h)] = vidx++;
                vertex_indices[to_vertex(h)] = vidx++;
                vidx += 4;
            }
        }

        // clean up
        if (vnormals)
            remove_vertex_property(vnormals);
        if (fnormals)
            remove_face_property(fnormals);
        remove_vertex_property(tmp_normal);
        remove_vertex_property(tmp_tex);
    }

    // we have a point cloud
    else if (n_vertices())
    {
        auto position = vertex_property<Point>("v:point");
        if (position)
        {
            positionArray.reserve(n_vertices());
            for (auto v : vertices())
                positionArray.push_back((vec3)position[v]);
        }

        auto normals = get_vertex_property<Point>("v:normal");
        if (normals)
        {
            normalArray.reserve(n_vertices());
            for (auto v : vertices())
                normalArray.push_back((vec3)normals[v]);
        }
    }

    // upload vertices
    if (!positionArray.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer_);
        glBufferData(GL_ARRAY_BUFFER, positionArray.size() * 3 * sizeof(float),
                     positionArray.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(0);
        n_vertices_ = positionArray.size();
    }
    else
        n_vertices_ = 0;

    // upload normals
    if (!normalArray.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, normal_buffer_);
        glBufferData(GL_ARRAY_BUFFER, normalArray.size() * 3 * sizeof(float),
                     normalArray.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(1);
    }

    // upload texture coordinates
    if (!texArray.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, tex_coord_buffer_);
        glBufferData(GL_ARRAY_BUFFER, texArray.size() * 2 * sizeof(float),
                     texArray.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(2);
        have_texcoords_ = true;
    }
    else
        have_texcoords_ = false;

    if (!coeffArray.empty())
    {
        glBindBuffer(GL_ARRAY_BUFFER, coeff_buffer_);
        glBufferData(GL_ARRAY_BUFFER, coeffArray.size() * sizeof(float),
                     coeffArray.data(), GL_STATIC_DRAW);
        glVertexAttribPointer(3, 1, GL_FLOAT, GL_FALSE, 0, nullptr);
        glEnableVertexAttribArray(3);
    }

    // edge indices
    if (n_edges())
    {
        std::vector<unsigned int> edgeArray;
        edgeArray.reserve(n_edges());
        for (auto e : edges())
        {
            edgeArray.push_back(vertex_indices[vertex(e, 0)]);
            edgeArray.push_back(vertex_indices[vertex(e, 1)]);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     edgeArray.size() * sizeof(unsigned int), edgeArray.data(),
                     GL_STATIC_DRAW);
        n_edges_ = edgeArray.size();
    }
    else
        n_edges_ = 0;

    // feature edges
    auto efeature = get_edge_property<bool>("e:feature");
    if (efeature)
    {
        std::vector<unsigned int> features;

        for (auto e : edges())
        {
            if (efeature[e])
            {
                features.push_back(vertex_indices[vertex(e, 0)]);
                features.push_back(vertex_indices[vertex(e, 1)]);
            }
        }

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, feature_buffer_);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                     features.size() * sizeof(unsigned int), features.data(),
                     GL_STATIC_DRAW);
        n_features_ = features.size();
    }
    else
        n_features_ = 0;

    // unbind vertex arry
    glBindVertexArray(0);

    // remove vertex index property again
    remove_vertex_property(vertex_indices);
}

void MySurfaceMeshGL::draw(const mat4 &projection_matrix,
                           const mat4 &modelview_matrix,
                           const std::string draw_mode)
{
    // did we generate buffers already?
    if (!vertex_array_object_)
    {
        update_opengl_buffers();
    }

    // load shader?
    if (!phong_shader_.is_valid())
    {
        if (!phong_shader_.source(phong_vshader, phong_fshader))
            exit(1);
    }

    // load shader?
    if (!phong_tess_shader_.is_valid())
    {
        if (!phong_tess_shader_.source(phong_tess_vshader, phong_tess_fshader,
                                       phong_tess_gshader, phong_tess_tcshader,
                                       phong_tess_teshader))
            exit(1);
    }

    // load shader?
    if (!matcap_shader_.is_valid())
    {
        if (!matcap_shader_.source(matcap_vshader, matcap_fshader))
            exit(1);
    }

    // we need some texture, otherwise WebGL complains
    if (!texture_)
    {
        use_cold_warm_texture();
    }

    // empty mesh?
    if (is_empty())
        return;

    // allow for transparent objects
    glEnable(GL_SAMPLE_ALPHA_TO_COVERAGE);

    // setup matrices
    mat4 p_matrix = projection_matrix;
    mat4 mv_matrix = modelview_matrix;
    mat4 mvp_matrix = projection_matrix * modelview_matrix;
    mat3 n_matrix = inverse(transpose(linear_part(mv_matrix)));

    // setup shader
    phong_tess_shader_.use();
    phong_tess_shader_.set_uniform("modelview_matrix", mv_matrix);
    phong_tess_shader_.set_uniform("normal_matrix", n_matrix);
    phong_tess_shader_.set_uniform("point_size", 5.0f);
    phong_tess_shader_.set_uniform("tess_level", 16.0f);
    phong_tess_shader_.set_uniform("projection_matrix", p_matrix);
    phong_tess_shader_.set_uniform("use_basis_function_texture", false);
    phong_tess_shader_.set_uniform("use_basis_function_elevation", false);
    phong_tess_shader_.set_uniform("basis_function_elevation",
                                   basis_elevation_);
    phong_tess_shader_.set_uniform("draw_edges", false);
    phong_tess_shader_.set_uniform("light1", vec3(1.0, 1.0, 1.0));
    phong_tess_shader_.set_uniform("light2", vec3(-1.0, 1.0, 1.0));
    phong_tess_shader_.set_uniform("front_color", front_color_);
    phong_tess_shader_.set_uniform("back_color", back_color_);
    phong_tess_shader_.set_uniform("ambient", ambient_);
    phong_tess_shader_.set_uniform("diffuse", diffuse_);
    phong_tess_shader_.set_uniform("specular", specular_);
    phong_tess_shader_.set_uniform("shininess", shininess_);
    phong_tess_shader_.set_uniform("alpha", alpha_);
    phong_tess_shader_.set_uniform("use_lighting", true);
    phong_tess_shader_.set_uniform("use_texture", false);
    phong_tess_shader_.set_uniform("use_srgb", false);
    phong_tess_shader_.set_uniform("show_texture_layout", false);

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
    phong_shader_.set_uniform("use_lighting", true);
    phong_shader_.set_uniform("use_texture", false);
    phong_shader_.set_uniform("use_srgb", false);
    phong_shader_.set_uniform("show_texture_layout", false);

    glBindVertexArray(vertex_array_object_);

    if (draw_mode == "Points")
    {
        phong_shader_.use();
#ifndef __EMSCRIPTEN__
        glEnable(GL_PROGRAM_POINT_SIZE);
#endif
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
        glDrawElements(GL_POINTS, 2 * n_edges_, GL_UNSIGNED_INT, nullptr);
    }
    else if (draw_mode == "Hidden Line")
    {
        if (n_faces())
        {
            phong_tess_shader_.use();
            // draw faces
            glDepthRange(0.01, 1.0);
            glPatchParameteri(GL_PATCH_VERTICES, 6);
            glDrawArrays(GL_PATCHES, 0, n_vertices_);
            glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

            // overlay edges
            glDepthRange(0.0, 1.0);
            glDepthFunc(GL_LEQUAL);
            phong_shader_.use();
            phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("use_lighting", false);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
            glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
            glDepthFunc(GL_LESS);
        }
    }
    else if (draw_mode == "Smooth Shading")
    {
        if (n_faces())
        {
            phong_tess_shader_.use();
            glPatchParameteri(GL_PATCH_VERTICES, 6);
            glDrawArrays(GL_PATCHES, 0, n_vertices_);
        }
    }
    else if (draw_mode == "Texture")
    {
        if (n_faces())
        {
            if (texture_mode_ == MatCapTexture)
            {
                matcap_shader_.use();
                matcap_shader_.set_uniform("modelview_projection_matrix",
                                           mvp_matrix);
                matcap_shader_.set_uniform("normal_matrix", n_matrix);
                matcap_shader_.set_uniform("alpha", alpha_);
                glBindTexture(GL_TEXTURE_2D, texture_);
                glDrawArrays(GL_PATCHES, 0, n_vertices_);
            }
            else
            {

                phong_tess_shader_.use();
                phong_tess_shader_.set_uniform("front_color",
                                               vec3(0.9, 0.9, 0.9));
                phong_tess_shader_.set_uniform("back_color",
                                               vec3(0.3, 0.3, 0.3));
                phong_tess_shader_.set_uniform("use_texture", true);
                phong_tess_shader_.set_uniform("use_srgb", srgb_);
                glPatchParameteri(GL_PATCH_VERTICES, 6);
                glBindTexture(GL_TEXTURE_2D, texture_);
                glDrawArrays(GL_PATCHES, 0, n_vertices_);
            }
        }
    }
    else if (draw_mode == "Texture Line")
    {
        if (n_faces())
        {

            glDepthRange(0.01, 1.0);
            if (texture_mode_ == MatCapTexture)
            {
                matcap_shader_.use();
                matcap_shader_.set_uniform("modelview_projection_matrix",
                                           mvp_matrix);
                matcap_shader_.set_uniform("normal_matrix", n_matrix);
                matcap_shader_.set_uniform("alpha", alpha_);
                glBindTexture(GL_TEXTURE_2D, texture_);
                glDrawArrays(GL_PATCHES, 0, n_vertices_);
            }
            else
            {
                phong_tess_shader_.use();
                phong_tess_shader_.set_uniform("front_color",
                                               vec3(0.9, 0.9, 0.9));
                phong_tess_shader_.set_uniform("back_color",
                                               vec3(0.3, 0.3, 0.3));
                phong_tess_shader_.set_uniform("use_texture", true);
                phong_tess_shader_.set_uniform("use_srgb", srgb_);
                glBindTexture(GL_TEXTURE_2D, texture_);
                glDrawArrays(GL_PATCHES, 0, n_vertices_);
            }

            // overlay edges
            glDepthRange(0.0, 1.0);
            glDepthFunc(GL_LEQUAL);
            phong_shader_.use();
            phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("use_lighting", false);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
            glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
            glDepthFunc(GL_LESS);
        }
    }
    else if (draw_mode == "Basis Functions")
    {
        phong_tess_shader_.use();
        phong_tess_shader_.set_uniform("use_texture", true);
        phong_tess_shader_.set_uniform("use_basis_function_texture", true);
        phong_tess_shader_.set_uniform("use_basis_function_elevation", true);
        glPatchParameteri(GL_PATCH_VERTICES, 6);
        glBindTexture(GL_TEXTURE_2D, texture_);
        glDrawArrays(GL_PATCHES, 0, n_vertices_);

        if (n_faces())
        {
            // draw faces
            glDepthRange(0.01, 1.0);
            glPatchParameteri(GL_PATCH_VERTICES, 6);
            glDrawArrays(GL_PATCHES, 0, n_vertices_);
            glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

            // overlay edges
            glDepthRange(0.0, 1.0);
            glDepthFunc(GL_LEQUAL);
            phong_shader_.use();
            phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("use_lighting", false);
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
            glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
            glDepthFunc(GL_LESS);
        }
    }
    else if (draw_mode == "Texture Layout")
    {
        std::cout << "Texture Layout currently unavailable" << std::endl;
        if (n_faces() && have_texcoords_ && false)
        {
            phong_shader_.set_uniform("show_texture_layout", true);
            phong_shader_.set_uniform("use_lighting", false);

            // draw faces
            phong_shader_.set_uniform("front_color", vec3(0.8, 0.8, 0.8));
            phong_shader_.set_uniform("back_color", vec3(0.9, 0.0, 0.0));
            glDepthRange(0.01, 1.0);
            glDrawArrays(GL_TRIANGLES, 0, n_vertices_);

            // overlay edges
            glDepthRange(0.0, 1.0);
            glDepthFunc(GL_LEQUAL);
            phong_shader_.set_uniform("front_color", vec3(0.1, 0.1, 0.1));
            phong_shader_.set_uniform("back_color", vec3(0.1, 0.1, 0.1));
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, edge_buffer_);
            glDrawElements(GL_LINES, n_edges_, GL_UNSIGNED_INT, nullptr);
            glDepthFunc(GL_LESS);
        }
    }

    // draw feature edges
    if (n_features_)
    {
        phong_shader_.set_uniform("front_color", vec3(0, 1, 0));
        phong_shader_.set_uniform("back_color", vec3(0, 1, 0));
        phong_shader_.set_uniform("use_lighting", false);
        glDepthRange(0.0, 1.0);
        glDepthFunc(GL_LEQUAL);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, feature_buffer_);
        glDrawElements(GL_LINES, n_features_, GL_UNSIGNED_INT, nullptr);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        glDepthFunc(GL_LESS);
    }

    // disable transparency (doesn't work well with imgui)
    glDisable(GL_SAMPLE_ALPHA_TO_COVERAGE);

    glBindVertexArray(0);
    glCheckError();
}

} // namespace pmp