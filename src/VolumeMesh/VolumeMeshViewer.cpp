//=============================================================================
// Copyright 2021 Hendrik Meyer, Astrid Bunge, Mario Botsch
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "VolumeMeshViewer.h"
#include <imgui.h>
#include <iostream>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

typedef OpenVolumeMesh::Geometry::Vec3f Vec3f;
typedef OpenVolumeMesh::VertexHandle VHandle;

//=============================================================================

VolumeMeshViewer::VolumeMeshViewer(const char* title, int width, int height,
                                   bool showgui)
    : pmp::TrackballViewer(title, width, height, showgui)
{
    // setup draw modes
    clear_draw_modes();
    add_draw_mode("Points");
    add_draw_mode("Line");
    add_draw_mode("Hidden Line");
    add_draw_mode("Texture");
    add_draw_mode("Embedded");
    add_draw_mode("Smooth Shading");

    set_draw_mode("Smooth Shading");

    embedding_ = 0.1f;
    alpha_ = 1.0f;

    // add help items
    add_help_item("Backspace", "Reload mesh", 3);

    //file_manager_.setVerbosityLevel(0);
}

//-----------------------------------------------------------------------------

VolumeMeshViewer::~VolumeMeshViewer() = default;

//-----------------------------------------------------------------------------

bool VolumeMeshViewer::load_mesh(const char* filename)
{
    if (mesh_.read(filename))
    {
        update_mesh();
        int Dof = 0;
        for (auto v : mesh_.vertices())
        {
            if (!mesh_.is_boundary(v))
            {
                Dof++;
            }
        }
        std::cout << "linear Dof: " << Dof << "\n";
        for (auto e : mesh_.edges())
        {
            if (!mesh_.is_boundary(e))
            {
                Dof++;
            }
        }
        std::cout << "quadratic Dof: " << Dof << "\n";
        auto bounding_sphere = mesh_.get_bounding_sphere();
        set_scene(bounding_sphere.first, bounding_sphere.second);

        mesh_.use_gradient_texture();

        set_draw_mode("Smooth Shading");

        filename_ = filename;
        return true;
    }

    std::cerr << "Failed to read mesh from " << filename << " !" << std::endl;
    return false;
}


void VolumeMeshViewer::update_mesh()
{
    mesh_.update_bounding_sphere();
    mesh_.update_opengl_buffers();
}


void VolumeMeshViewer::process_imgui()
{
    if (ImGui::CollapsingHeader("Mesh Info", ImGuiTreeNodeFlags_DefaultOpen))
    {
        // output mesh statistics
        ImGui::BulletText("%d vertices", (int)mesh_.n_logical_vertices());
        ImGui::BulletText("%d edges", (int)mesh_.n_logical_edges());
        ImGui::BulletText("%d faces", (int)mesh_.n_logical_faces());
        ImGui::BulletText("%d cells", (int)mesh_.n_logical_cells());
    }

    if (ImGui::CollapsingHeader("Controls", ImGuiTreeNodeFlags_DefaultOpen))
    {
        //control embedding
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Cell embedding", &embedding_, 0.0, 1.0);
        ImGui::PopItemWidth();
        if (embedding_ != mesh_.embedding())
        {
            if (embedding_ > 0)
            {
                mesh_.set_embedding(embedding_);
            }
            else
            {
                mesh_.set_embedding(1e-4);
            }
        }

        //control near clipping plane
        ImGui::PushItemWidth(100);
        float near = near_clipping_;
        ImGui::SliderFloat("Near clipping", &near, 0.0, 1.0);
        ImGui::PopItemWidth();
        if (near != near_clipping_)
        {
            near_clipping_ = near;

            pmp::vec4 mc(center_, 1.0);
            pmp::vec4 ec = modelview_matrix_ * mc;
            float z = -ec[2];
            near = std::max(0.001f * radius_,
                            z - radius_ + 2 * near_clipping_ * radius_);

            std::vector<std::vector<HFHandle>> cells;
            std::vector<CHandle> cell_del;
            std::vector<bool> cell_tag, edge_tag, vertex_tag;
            for (auto c : mesh_.cells())
            {
                for (auto c_v : mesh_.cell_vertices(c))
                {
                    pmp::vec4 m_vertex((float)mesh_.vertex(c_v)[0],
                                       (float)mesh_.vertex(c_v)[1],
                                       (float)mesh_.vertex(c_v)[2], 1.0);
                    pmp::vec4 e_vertex = modelview_matrix_ * m_vertex;
                    if (-e_vertex[2] < near)
                    {
                        mesh_.delete_vertex(c_v);
                    }
                }
            }
            mesh_.collect_garbage();
            mesh_.update_opengl_buffers();
        }

        // control opacity
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Opacity", &alpha_, 0.0, 1.0);
        ImGui::PopItemWidth();
        if (alpha_ != mesh_.alpha())
        {
            mesh_.set_alpha(alpha_);
        }
    }
}

//-----------------------------------------------------------------------------

void VolumeMeshViewer::draw(const std::string& draw_mode)
{
    mesh_.draw(projection_matrix_, modelview_matrix_, draw_mode);
}

//-----------------------------------------------------------------------------

void VolumeMeshViewer::keyboard(int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS && action != GLFW_REPEAT)
        return;

    switch (key)
    {
        case GLFW_KEY_BACKSPACE: // reload model
        {
            load_mesh(filename_.c_str());
            break;
        }

        case GLFW_KEY_W: // write mesh to output.ovm
        {
            file_manager_.writeFile("output.ovm", mesh_);
            break;
        }

        default: {
            TrackballViewer::keyboard(key, scancode, action, mods);
            break;
        }
    }
}

//-----------------------------------------------------------------------------

void VolumeMeshViewer::mouse(int button, int action, int mods)
{
    TrackballViewer::mouse(button, action, mods);
}

//-----------------------------------------------------------------------------

void VolumeMeshViewer::display()
{
    // clear buffers
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // adjust clipping planes to tightly fit bounding sphere
    pmp::vec4 mc(center_, 1.0);
    pmp::vec4 ec = modelview_matrix_ * mc;
    float z = -ec[2];
    fovy_ = 45.0f;
    near_ = std::max(0.001f * radius_, z - radius_);
    // std::max(0.001f * radius_, z - radius_ + 2 * near_clipping_ * radius_);
    far_ = std::max(0.002f * radius_, z + radius_);

    // update projection matrix
    projection_matrix_ = pmp::perspective_matrix(
        fovy_, (float)width() / (float)height(), near_, far_);

    // draw the scene in current draw mode
    if (draw_mode_ < draw_mode_names_.size())
        draw(draw_mode_names_[draw_mode_]);
    else
        draw("");
}

//-----------------------------------------------------------------------------

VHandle VolumeMeshViewer::pick_vertex(double x, double y)
{
    VHandle vmin;

    vec3 p;
    float d, dmin(FLT_MAX);

    // get 3D position
    if (TrackballViewer::pick(x, y, p))
    {
        Vec3f point(p[0], p[1], p[2]);
        for (auto v_it = mesh_.v_iter(); v_it.valid(); ++v_it)
        {
            // compare 3D position to each vertex
            d = (mesh_.vertex(*v_it) - point).length();
            if (d < dmin)
            {
                dmin = d;
                vmin = *v_it;
            }
        }
    }
    return vmin;
}
