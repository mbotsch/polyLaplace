//=============================================================================

#include "SurfaceViewer.h"
#include "Surface/diffgeo.h"
#include "Surface/Curvature.h"
#include "Surface/SpectralProcessing.h"
#include <pmp/algorithms/Triangulation.h>
#include <pmp/algorithms/Subdivision.h>
#include <pmp/algorithms/Normals.h>
#include <imgui.h>
#include "Surface/GeodesicsInHeat.h"
#include "Surface/[dGBD20]Laplace.h"
#include "Surface/Poisson_System.h"
#include <pmp/Timer.h>
#include <random>

using namespace pmp;

//=============================================================================

enum InsertedPoint
{
    Centroid_ = 0,
    AreaMinimizer = 2,
};

enum LaplaceMethods {
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
    Harmonic = 4,
};

void Viewer::keyboard(int key, int scancode, int action, int mods)
{
    if (action != GLFW_PRESS) // only react on key press events
        return;

    switch (key)
    {
        default: {
            MeshViewer::keyboard(key, scancode, action, mods);
            break;
        }
    }
}

//----------------------------------------------------------------------------

void Viewer::process_imgui()
{
    // add standard mesh info stuff
    pmp::MeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::CollapsingHeader("Settings", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Hyperparameter Alexa & Wardetzky Laplace",
                           &poly_laplace_lambda_, 0.01, 3.0);
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Hyperparameter deGoes Laplace",
                           &deGoes_laplace_lambda_, 0.01, 3.0);
        ImGui::PopItemWidth();

        ImGui::Spacing();

        ImGui::Text("Choose your Laplacian");

        ImGui::Spacing();

        static int laplace = 0;
        ImGui::RadioButton("Polysimple Laplace", &laplace, 0);
        ImGui::RadioButton("Alexa & Wardetzky Laplace", &laplace, 1);
        ImGui::RadioButton("Diamond", &laplace, 2);
        ImGui::RadioButton("deGoes Laplace", &laplace, 3);
        ImGui::RadioButton("Harmonic", &laplace, 4);

        ImGui::Spacing();

        ImGui::Text("Choose your minimizing Point ");

        ImGui::Spacing();

        static int min_point = 2;
        ImGui::RadioButton("Centroid", &min_point, 0);
        ImGui::RadioButton("Area Minimizer", &min_point, 2);

        static int ts = 0;
        ImGui::Text("Choose your diffusion time step ");

        ImGui::RadioButton("Mean edge length", &ts, 0);
        ImGui::RadioButton("Max edge length", &ts, 1);
        ImGui::RadioButton("Max diagonal length", &ts, 2);

        laplace_matrix = laplace;
        min_point_ = min_point;
        time_step_ = DiffusionStep(ts);
    }
    ImGui::Spacing();
    ImGui::Spacing();

    // turn mesh into non-triangles
    if (ImGui::CollapsingHeader("Polygons!"))
    {
        // Catmull-Clark subdivision
        if (ImGui::Button("Catmull-Clark"))
        {
            Subdivision(mesh_).catmull_clark();
            update_mesh();
        }
        if (ImGui::Button("insert virtual points"))
        {
            insert_points(min_point_);
        }
        // Catmull-Clark subdivision
        if (ImGui::Button("Centroids"))
        {
            Centroid();
        }
        // dualize the mesh
        if (ImGui::Button("Dualize mesh"))
        {
            dualize();
        }

        if (ImGui::Button("Relocate to mesh center"))
        {
            Point p = centroid(mesh_);
            std::cout << "Center Point: " << p << std::endl;
            for (auto v : mesh_.vertices())
            {
                mesh_.position(v) -= p;
            }
            std::cout << "new Center Point: " << centroid(mesh_) << std::endl;

            update_mesh();
        }
        // close holes by polygons
        if (ImGui::Button("Close holes"))
        {
            close_holes();
        }

        if (ImGui::Button("Triangulate mesh (min area)"))
        {
            Triangulation tesselator(mesh_);
            tesselator.triangulate(Triangulation::Objective::MIN_AREA);
            update_mesh();
        }

        if (ImGui::Button("Triangulate mesh (max angle)"))
        {
            Triangulation tesselator(mesh_);
            tesselator.triangulate(Triangulation::Objective::MAX_ANGLE);
            update_mesh();
        }
        if (ImGui::Button("Kugelize"))
        {
            for (auto v : mesh_.vertices())
                mesh_.position(v) = normalize(mesh_.position(v));
            update_mesh();
        }
        if (ImGui::Button("Add noise"))
        {
            for (auto v : mesh_.vertices())
            {
                Point n = Normals::compute_vertex_normal(mesh_, v);
                Scalar r = 2.0 * static_cast<float>(rand()) /
                               static_cast<float>(RAND_MAX) -
                           1.0;
                Scalar h = r * 0.01 * radius_;
                mesh_.position(v) += h * n;
            }
            update_mesh();
        }
        static bool normal_, fixed_;

        ImGui::Checkbox("noise in normal direction?", &normal_);
        ImGui::Checkbox("Constant normal noise?", &fixed_);
        if (ImGui::Button("Noise"))
        {
            // create random generator
            // upper and lower bounds are proportional to bounding box and inverse of smoothness
            auto l = mesh_.bounds().size();
            double upper_bound = l / 1000.0;
            double lower_bound = -upper_bound;
            std::uniform_real_distribution<double> unif(lower_bound,
                                                        upper_bound);
            std::default_random_engine re;

            re.seed(42); // fixed seed

            double rand =0.0;
            if (fixed_)
            {
                rand = unif(re);
            }

            for (auto v : mesh_.vertices())
            {
                auto n = pmp::Normals::compute_vertex_normal(mesh_, v);
                if (normal_)
                {
                    if (fixed_)
                    {
                        mesh_.position(v) -= n * rand;
                    }
                    else
                    {
                        mesh_.position(v) -= n * unif(re);
                    }
                }
                else
                {
                    mesh_.position(v) +=
                        pmp::Point(unif(re), unif(re), unif(re));
                }
            }

            update_mesh();
        }
        if (ImGui::Button("collapse Edge"))
        {
            Edge origin_edge;
            auto min = DBL_MAX;
            for(auto e: mesh_.edges()){
                Vertex v1 = mesh_.vertex(e,1);
                Vertex v0 = mesh_.vertex(e,0);
                Point p = 0.5*(mesh_.position(v0)+mesh_.position(v1));
                Point origin(0.0,0.0,0.0);
                if(norm(origin-p)<min){
                    origin_edge = e;
                    min = norm(origin-p);
                }
            }
            Vertex v1 = mesh_.vertex(origin_edge,1);
            Vertex v0 = mesh_.vertex(origin_edge,0);
            Point dir = (mesh_.position(v1)-mesh_.position(v0));
            std::cout << mesh_.position(v0) << std::endl;
            mesh_.position(v0) -= 0.1*dir;
            std::cout << mesh_.position(v0) << std::endl;
            update_mesh();
        }
        if (ImGui::Button("Add tangential noise")) {
            // create random generator
            // upper and lower bounds are proportional to bounding box and inverse of smoothness
            auto l = mesh_.bounds().size();
            double upper_bound = l / 1000.0;
            double lower_bound = -upper_bound;
            std::uniform_real_distribution<double> unif(lower_bound,
                                                        upper_bound);
            std::default_random_engine re;

            re.seed(42); // fixed seed

            for (auto e: mesh_.edges()) {
                auto n = mesh_.position(mesh_.vertex(e, 0)) - mesh_.position(mesh_.vertex(e, 1));
                n.normalize();
                mesh_.position(mesh_.vertex(e, 0)) += n * unif(re) ;

            }
            update_mesh();
        }
        if (ImGui::Button("Check for non-planarity")) {
            auto vpoint = mesh_.get_vertex_property<Point>("v:point");
            int ctr = 0;
            double max_sum = 0.0;
            for (auto f: mesh_.faces()) {
                Point p0, p1, p2;
                for(auto h :mesh_.halfedges(f)) {
                    p0 = vpoint[mesh_.from_vertex(h)];
                    p1 = vpoint[mesh_.to_vertex(h)];
                    p2 = vpoint[mesh_.to_vertex(mesh_.next_halfedge(h))];
                    break;
                }
                Point n = normalize(cross(p0-p1,p2-p1));
                for (auto v: mesh_.vertices(f)) {
                    Point p = vpoint[v];
                    if (abs(dot(n, p - p1)) > 0.0001) {
                        ctr += 1;
                        break;
                    }
                }

                double max = -1000000;
                for(auto h :mesh_.halfedges(f)) {
                    p0 = vpoint[mesh_.from_vertex(h)];
                    p1 = vpoint[mesh_.to_vertex(h)];
                    p2 = vpoint[mesh_.to_vertex(mesh_.next_halfedge(h))];

                    n = normalize(cross(p0-p1,p2-p1));
                    for (auto v: mesh_.vertices(f)) {
                        Point p = vpoint[v];
                        if (abs(dot(n, p - p1)) > max) {
                            max = abs(dot(n, p - p1));
                        }
                    }
                }
                max_sum += max;
            }
            std::cout << "Nr. non planar faces: " << ctr << std::endl;
            std::cout << "Mean maximum plane distance " << max_sum/mesh_.n_faces() << std::endl;
        }
    }

    ImGui::Spacing();
    ImGui::Spacing();

    // discrete harmonic parameterization
    if (ImGui::CollapsingHeader("Spectral Processing",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        if (ImGui::Button("SH Reproduction"))
        {
            bool lumped = true;
            if (laplace_matrix == Diamond)
                lumped = false;
            rmse_sh(mesh_, laplace_matrix, min_point_, lumped);
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Poisson System",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::PushItemWidth(100);
        static int function = 2;
        ImGui::RadioButton("Franke 2D (planar)", &function, 2);
        ImGui::RadioButton("Spherical Harmonics", &function, 0);
        ImGui::PushItemWidth(100);
        static int l = 4;
        static int m = 2;
        ImGui::SliderInt("SH l (degree):",&l, 0, 5);
        ImGui::SliderInt("SH m (band):",&m, -l, l);
        ImGui::PopItemWidth();
        if (ImGui::Button("Solve!"))
        {
            solve_poisson_system(mesh_, laplace_matrix, min_point_,function, l, m);
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();

    // curvature visualization
    if (ImGui::CollapsingHeader("Curvature", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static bool curvature_sphere_ = false;
        ImGui::Checkbox("Compare to unit sphere curvatures",
                        &curvature_sphere_);

        if (ImGui::Button("Mean Curvature"))
        {
            Curvature analyzer(mesh_, curvature_sphere_);

                analyzer.visualize_curvature(laplace_matrix, min_point_, true);
            mesh_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
    if (ImGui::CollapsingHeader("Geodesics in Heat",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {

        static int geodesics = 2;
        ImGui::RadioButton("Compare distances to arc lengths", &geodesics,0);
        ImGui::RadioButton("Compare to euclidean distances", &geodesics,1);
        ImGui::RadioButton("No comparison", &geodesics,2);

        if(geodesics == 0){
            compare_sphere = true;
            compare_cube = false;
        }else if (geodesics == 1){
            compare_sphere = false;
            compare_cube = true;
        }else if (geodesics == 2){
            compare_sphere = false;
            compare_cube = false;
        }
        if (ImGui::Button("Compute Distances Vertex 0"))
        {

            GeodesicsInHeat heat(mesh_, laplace_matrix, min_point_,
                                 compare_sphere, compare_cube, time_step_);
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance(0, dist, geodist);

            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");

            mesh_.use_checkerboard_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
}

//----------------------------------------------------------------------------

void Viewer::draw(const std::string &draw_mode)
{
    // normal mesh draw
    mesh_.draw(projection_matrix_, modelview_matrix_, draw_mode);
}

//----------------------------------------------------------------------------

void Viewer::dualize()
{
    SurfaceMeshGL dual;

    auto fvertex = mesh_.add_face_property<Vertex>("f:vertex");
    for (auto f : mesh_.faces())
    {
        fvertex[f] = dual.add_vertex(centroid(mesh_, f));
    }

    for (auto v : mesh_.vertices())
    {
        if (!mesh_.is_boundary(v))
        {
            std::vector<Vertex> vertices;
            for (auto f : mesh_.faces(v))
                vertices.push_back(fvertex[f]);
            dual.add_face(vertices);
        }
    }

    mesh_ = dual;
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::update_mesh()
{
    // re-compute face and vertex normals
    mesh_.update_opengl_buffers();
}

//----------------------------------------------------------------------------

void Viewer::close_holes()
{
    bool finished = false;
    std::vector<Face> holes;
    while (!finished)
    {
        finished = true;

        // loop through all vertices
        for (auto v : mesh_.vertices())
        {
            // if we find a boundary vertex...
            if (mesh_.is_boundary(v))
            {
                // trace boundary loop
                std::vector<Vertex> vertices;
                vertices.push_back(v);
                for (Halfedge h = mesh_.halfedge(v); mesh_.to_vertex(h) != v;
                     h = mesh_.next_halfedge(h))
                {
                    vertices.push_back(mesh_.to_vertex(h));
                }

                // add boudary loop as polygonal face
                Face f = mesh_.add_face(vertices);
                holes.push_back(f);
                // start over
                finished = false;
                break;
            }
        }
    }
    update_mesh();
}

//----------------------------------------------------------------------------
void Viewer::insert_points(unsigned int minpoint)
{

    Eigen::MatrixXd poly;
    Eigen::VectorXd w;

    for (Face f : mesh_.faces())
    {
        const int n = (int)mesh_.valence(f);
        poly.resize(n, 3);
        int i = 0;
        for (Vertex v : mesh_.vertices(f))
        {
            for (int h = 0; h < 3; h++)
            {
                poly.row(i)(h) = mesh_.position(v)[h];
            }
            i++;
        }

        // compute weights for the polygon
        Eigen::Vector3d p;
        if (minpoint == Centroid_)
        {
            auto val = poly.rows();
            w = Eigen::MatrixXd::Ones(val, 1);
            w /= (double)val;
        }
        else
        {
            find_area_minimizer_weights(poly, w);
        }
        Eigen::Vector3d point = poly.transpose() * w;
        Vertex ver = mesh_.add_vertex(Point(point(0), point(1), point(2)));
        //        std::cout << "norm of point: " << norm(Point(point(0), point(1), point(2)))  << std::endl;
        mesh_.split(f, ver);
    }
    mesh_.garbage_collection();
    update_mesh();
}
//----------------------------------------------------------------------------

void Viewer::Centroid()
{
    for (Face f : mesh_.faces())
    {
        Point c = centroid(mesh_, f);
        Vertex hidx = mesh_.add_vertex(c);
        mesh_.split(f, hidx);
    }
    mesh_.garbage_collection();
    update_mesh();
}

void Viewer::mouse(int button, int action, int mods)
{
    if (action == GLFW_PRESS && button == GLFW_MOUSE_BUTTON_MIDDLE &&
        mods == GLFW_MOD_SHIFT)
    {
        double x, y;
        cursor_pos(x, y);
        Vertex v = pick_vertex((int)x, (int)y);
        std::cout << "Vertex Idx: " << v.idx() << std::endl;
        if (mesh_.is_valid(v))
        {
            GeodesicsInHeat heat(mesh_, laplace_matrix, min_point_,
                                 compare_sphere, compare_cube, time_step_);
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance((int)v.idx(), dist, geodist);
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }else
    {
        MeshViewer::mouse(button, action, mods);
    }
}

void Viewer::load_mesh(const char *filename)
{
    MeshViewer::load_mesh(filename);
    set_draw_mode("Hidden Line");
}

//=============================================================================
