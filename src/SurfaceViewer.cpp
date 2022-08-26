//=============================================================================

#include "SurfaceViewer.h"
#include "Surface/diffgeo.h"
#include "Surface/Curvature.h"
#include "Surface/SpectralProcessing.h"
#include <pmp/algorithms/SurfaceTriangulation.h>
#include <pmp/algorithms/SurfaceSubdivision.h>
#include <pmp/algorithms/SurfaceNormals.h>
#include <imgui.h>
#include "Surface/GeodesicsInHeat.h"
#include "Surface/DisneyLaplace.h"
#include "Surface/Poisson_System.h"
#include "Surface/SECLaplace.h"
#include "Surface/AQAPoly_Laplacian.h"
#include "Surface/SmoothSubdivBasis.h"
#include <pmp/Timer.h>
#include <random>
#include <Eigen/CholmodSupport>

using namespace pmp;

//=============================================================================

enum InsertedPoint
{
    Centroid_ = 0,
    AbsAreaMinimizer = 1,
    AreaMinimizer = 2,
    Triangle_Circumcenter = 3
};

enum LaplaceMethods
{

    SandwichLaplace = 0,
    AlexaLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    IntrinsicDelaunay = 4,
    Disney = 5,
    SEC = 6,
    AQAPoly_Laplace = 7
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
        //        ImGui::Checkbox("Clamp cotan", &clamp_cotan_);

        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Hyperparameter Alexa Laplace",
                           &poly_laplace_lambda_, 0.01, 3.0);
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Hyperparameter Disney Laplace",
                           &disney_laplace_lambda_, 0.01, 3.0);
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("Hyperparameter SEC Laplace", &sec_laplace_lambda_,
                           0.01, 3.0);
        ImGui::SliderInt("SEC Laplace lvl", &sec_laplace_lvl, 1, 8);
        ImGui::PopItemWidth();

        ImGui::Spacing();

        ImGui::Text("Choose your Laplacian");

        ImGui::Spacing();

        static int laplace = 0;
        ImGui::RadioButton("Sandwich Laplace", &laplace, 0);
        ImGui::RadioButton("Alexa Laplace", &laplace, 1);
        ImGui::RadioButton("Cotan Laplace", &laplace, 2);
        ImGui::RadioButton("Diamond", &laplace, 3);
        ImGui::RadioButton("Intrinsic Delaunay", &laplace, 4);
        ImGui::RadioButton("Disney Laplace", &laplace, 5);
        ImGui::RadioButton("New CC Laplace", &laplace, 6);
        ImGui::RadioButton("AQuadAp Laplace", &laplace, 7);
        ImGui::RadioButton("quadratic Triangle Laplace", &laplace, 8);
        ImGui::RadioButton("Harmonic", &laplace, 10);

        ImGui::Text("Degree of Basis functions");

        ImGui::Spacing();
        static int degree = 2;
        ImGui::RadioButton("Linear", &degree, 1);
        ImGui::RadioButton("Quadratic", &degree, 2);
        ImGui::RadioButton("Cubic", &degree, 3);
        ImGui::RadioButton("Quartic", &degree, 4);

        degree_ = degree;
        ImGui::Spacing();
        ImGui::Text("DoF to keep");

        ImGui::Spacing();
        static int coarsening = Edges;
        ImGui::RadioButton("Vertices", &coarsening, Vertices);
        ImGui::RadioButton("Edges", &coarsening, Edges);
        ImGui::RadioButton("Refined mesh", &coarsening, Refined_mesh);

        coarseningType_ = (CoarseDimension)coarsening;
        ImGui::Spacing();

        ImGui::Text("Choose your minimizing Point ");

        ImGui::Spacing();

        static int min_point = 2;
        if (ImGui::Button("Faces without point recreation"))
        {
            problem_faces(mesh_, min_point_);
            update_mesh();
        }
        ImGui::RadioButton("Centroid", &min_point, 0);
        ImGui::RadioButton("Area Minimizer", &min_point, 2);

        static int ts = 0;
        ImGui::Text("Choose your diffusion timestep ");

        ImGui::RadioButton("Mean edge length", &ts, 0);
        ImGui::RadioButton("Max edge length", &ts, 1);

        laplace_matrix = laplace;
        min_point_ = min_point;
        if (ts == 0)
        {
            time_step_ = true;
        }
        else
        {
            time_step_ = false;
        }
    }
//    if (ImGui::Button("Print simpel Prolongation"))
//    {
//        Eigen::SparseMatrix<double> P;
//        setup_simple_2D_prolongation(mesh_, P);
//    }
//    if (ImGui::Button("Test Quadratic Reproduction"))
//    {
//        test_quadratic_reproduction(mesh_, coarseningType_, degree_);
//    }
//    if (ImGui::Button("Test Multigrid"))
//    {
//        solve_AQAPoly_Poisson_mg(mesh_, coarseningType_, degree_, false);
//    }
//    if (ImGui::Button("Mean edge"))
//    {
//        double mean = 0.0;
//        for (auto e : mesh_.edges())
//        {
//            mean += mesh_.edge_length(e);
//        }
//        std::cout << "mean edge length: " << mean / mesh_.n_edges()
//                  << std::endl;
//        std::cout << "Quad DOF: " << mesh_.n_edges() + mesh_.n_vertices()
//                  << std::endl;
//        std::cout << "Linear DOF: " << mesh_.n_vertices() << std::endl;
//    }
//    if (ImGui::Button("Linear Precision"))
//    {
//        solve_laplace_equation(mesh_, laplace_matrix, min_point_, degree_,
//                               coarseningType_);
//    }
//    if (ImGui::Button("Mirror Mesh"))
//    {
//        mirror_mesh();
//    }
    const int nv = basis_values.cols() - 1;
    static int nr_basis = 0;
    static int basis = 1;
    ImGui::SliderInt("Subdiv Laplace lvl", &subdiv_lvl, 1, 8);

    ImGui::RadioButton("CC_Subdiv_CC_P", &basis, 0);
    ImGui::RadioButton("Lin_Subdiv_CC_P", &basis, 1);
    ImGui::RadioButton("Lin_Subdiv_Lin_P", &basis, 2);
    ImGui::RadioButton("SEC", &basis, 3);

    ImGui::SliderInt("Basis vis", &nr_basis, 0, nv);
    if(ImGui::Button("Non Dirichlet Curvature")){
        curvature_non_dirichlet(mesh_, basis);
        mesh_.use_cold_warm_texture();
        update_mesh();
        set_draw_mode("Texture");
     }
    if(ImGui::Button("Non Dirichlet Poisson")){
        solve_poisson_non_dirichlet(mesh_,1, basis);
    }
    if(ImGui::Button("Non Dirichlet Linear Precision")){
        solve_poisson_non_dirichlet(mesh_,0,basis);
        mesh_.use_cold_warm_texture();
        update_mesh();
//        set_draw_mode("Basis Shading Grid");
        set_draw_mode("Solution Error");

    }
    if (ImGui::Button("CC Basis"))
    {
        Eigen::SparseMatrix<double> S,M;
        setup_smooth_basis_matrix(mesh_, basis_values, basis);
        if(basis != 0)
        {
            for (int i = 0; i < subdiv_lvl; i++)
            {
                linear_interpolation_catmull_clark(mesh_);
            }
        }else{
            SurfaceSubdivision divider = SurfaceSubdivision(mesh_);

            for (int i = 0; i < subdiv_lvl; i++)
            {
                divider.catmull_clark();
            }
        }
        update_mesh();
    }
    if (ImGui::Button("visualise CC Basis"))
    {

        // sort basis values
//        std::cout << basis_values << std::endl;
        std::vector<Scalar> values;
        auto basis_i = mesh_.add_vertex_property<Scalar>("v:basis");
        values.reserve(mesh_.n_vertices());
        for (auto v : mesh_.vertices())
        {
            values.push_back(basis_values(v.idx(), nr_basis));
            basis_i[v] = basis_values(v.idx(), nr_basis);
        }
        std::sort(values.begin(), values.end());
        unsigned int n = values.size() - 1;

        // generate 1D texture coordiantes
        auto tex = mesh_.vertex_property<TexCoord>("v:tex");
        double sum = 0.0;
        for (auto v : mesh_.vertices())
        {
//            tex[v] = TexCoord((basis_i[v] - kmin) / (kmax - kmin), 0.0);
                tex[v] = TexCoord(basis_i[v] , 0.0);
                if(v.idx() < basis_values.cols()){
                    std::cout << v.idx() << ": "<< tex[v] << std::endl;
                    sum += tex[v][0];
                }
        }
        std::cout << "Sum function values: " << sum << std::endl;
        std::cout << "---------------------------------"  << std::endl;
        Scalar bb_size = mesh_.bounds().size();
        for (auto v : mesh_.vertices())
        {
            tex[v] = TexCoord(basis_i[v] / bb_size, 0.0);
        }

//         remove per-halfedge texture coordinates
        auto htex = mesh_.get_halfedge_property<TexCoord>("h:tex");
        if (htex)
            mesh_.remove_halfedge_property(htex);
        mesh_.use_cold_warm_texture();
//        mesh_.use_checkerboard_texture();
        update_mesh();
        set_draw_mode("Basis Shading Grid");


        mesh_.remove_vertex_property(basis_i);
    }
    ImGui::Spacing();
    ImGui::Spacing();

    if (ImGui::Button("Condition Nr")){
        std::vector<double> condnr;
        std::cout << "Coarsening: " << coarseningType_ << " Degree: " << degree_ << std::endl;
        AQAPoly_condition_nr(mesh_,  coarseningType_, degree_, condnr);

    }

    // turn mesh into non-triangles
    if (ImGui::CollapsingHeader("Polygons!"))
    {
        // Catmull-Clark subdivision
        if (ImGui::Button("Catmull-Clark"))
        {
            SurfaceSubdivision(mesh_).catmull_clark();
            update_mesh();
        }
        if (ImGui::Button("Modified Butterfly ##2"))
        {
            Eigen::SparseMatrix<double> P_cc;
            setup_mod_butterfly_P_matrix(mesh_, P_cc, 1);
            update_mesh();
        }

        if (ImGui::Button("Interpolating Catmull-Clark"))
        {
            linear_interpolation_catmull_clark(mesh_);
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
            SurfaceTriangulation tesselator(mesh_);
            tesselator.triangulate(SurfaceTriangulation::Objective::MIN_AREA);
            update_mesh();
        }

        if (ImGui::Button("Triangulate mesh (max angle)"))
        {
            SurfaceTriangulation tesselator(mesh_);
            tesselator.triangulate(SurfaceTriangulation::Objective::MAX_ANGLE);
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
                Point n = SurfaceNormals::compute_vertex_normal(mesh_, v);
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

            double rand;
            if (fixed_)
            {
                rand = unif(re);
            }

            for (auto v : mesh_.vertices())
            {
                auto n = pmp::SurfaceNormals::compute_vertex_normal(mesh_, v);
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
            rmse_sh(mesh_, laplace_matrix, min_point_, lumped, degree_,
                    coarseningType_);
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
        //        ImGui::RadioButton("Franke Halfsphere", &function, 3);

        ImGui::PopItemWidth();

        if (ImGui::Button("Solve!"))
        {
            if (laplace_matrix == AQAPoly_Laplace && function == 2)
            {
                solve_AQAPoly_Poisson_mg(mesh_, coarseningType_, degree_, true);
                //                solve_poisson_system(mesh_, laplace_matrix, min_point_,
                //                     function, degree_, coarseningType_, 2, 0);
            }
            else if (laplace_matrix == AQAPoly_Laplace)
            {
                //                solve_poisson_system(mesh_, laplace_matrix, min_point_,
                //                                     function, 2, Edges, 2, 0);
                solve_AQAPoly_Poisson_mg(mesh_, coarseningType_, degree_, true);
            }
            else
            {
                solve_poisson_system(mesh_, laplace_matrix, min_point_,
                                     function, 1, Vertices, 2, 0);
            }
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();

    // implicit smoothing
    if (ImGui::CollapsingHeader("Smoothing", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static float timestep = 0.1;
        float lb = 0.001;
        float ub = 1.0;
        ImGui::PushItemWidth(100);
        ImGui::SliderFloat("TimeStep", &timestep, lb, ub);
        ImGui::PopItemWidth();

        if (ImGui::Button("Implicit Smoothing"))
        {
            close_holes();

            Scalar dt = timestep;
            smooth_.implicit_smoothing_misha(dt, laplace_matrix, min_point_);
            update_mesh();
            BoundingBox bb = mesh_.bounds();
            set_scene((vec3)bb.center(), 0.5 * bb.size());
            open_holes();
        }
        if (ImGui::Button("20 * Implicit Smoothing"))
        {
            close_holes();

            Scalar dt = timestep;
            for (int i = 0; i < 20; i++)
            {
                smooth_.implicit_smoothing_misha(dt, laplace_matrix,
                                                 min_point_);
                update_mesh();
            }
            BoundingBox bb = mesh_.bounds();
            set_scene((vec3)bb.center(), 0.5 * bb.size());
            open_holes();
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
            if (laplace_matrix == AQAPoly_Laplace)
            {
                analyzer.visualize_curvature(laplace_matrix, min_point_, true,
                                             degree_, coarseningType_);
            }
            else
            {
                analyzer.visualize_curvature(laplace_matrix, min_point_, true,
                                             1);
            }
            mesh_.use_cold_warm_texture();
            update_mesh();
            set_draw_mode("Texture");
        }
    }
    if (ImGui::CollapsingHeader("Geodesics in Heat",
                                ImGuiTreeNodeFlags_DefaultOpen))
    {

        static bool geodesic_sphere_ = false;
        static bool geodesic_cube_ = false;
        ImGui::Checkbox("Compare distances to arc lengths", &geodesic_sphere_);
        ImGui::Checkbox("Compare to euclidean distances", &geodesic_cube_);
        compare_sphere = geodesic_sphere_;
        compare_cube = geodesic_cube_;
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
    holes_ = holes;
    update_mesh();
}

//----------------------------------------------------------------------------
void Viewer::insert_points(unsigned int minpoint)
{

    Eigen::MatrixXd poly;
    Eigen::VectorXd w;

    for (Face f : mesh_.faces())
    {
        const int n = mesh_.valence(f);
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
            int val = poly.rows();
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

void Viewer::mirror_mesh()
{

    double max_x = std::numeric_limits<double>::min();
    double min_x = std::numeric_limits<double>::max();

    for (auto v : mesh_.vertices())
    {
        if (mesh_.position(v)[0] > max_x)
        {
            max_x = mesh_.position(v)[0];
        }
        if (mesh_.position(v)[0] < min_x)
        {
            min_x = mesh_.position(v)[0];
        }
    }
    double d = max_x - min_x;
    //    std::cout << "min: " << min_x << std::endl;
    //    std::cout << "max: " << max_x << std::endl;
    //    std::cout << "D: " << d << std::endl;

    for (Face f : mesh_.faces())
    {
        std::vector<pmp::Vertex> face;
        for (Vertex v : mesh_.vertices(f))
        {
            pmp::Point vp = mesh_.position(v);
            //            std::cout << "Vertex position: " << vp << std::endl;
            if (abs(vp[0] - max_x) < 0.00001)
            {
                face.emplace_back(v);
            }
            else
            {
                bool new_vertex = true;
                for (auto vv : mesh_.vertices())
                {
                    if (norm(mesh_.position(vv) -
                             pmp::Point(max_x + (max_x - vp[0]), vp[1],
                                        vp[2])) < 0.0001)
                    {
                        face.emplace_back(vv);
                        new_vertex = false;
                        break;
                    }
                }
                if (new_vertex)
                {
                    pmp::Point np =
                        pmp::Point(max_x + (max_x - vp[0]), vp[1], vp[2]);
                    pmp::Vertex nv = mesh_.add_vertex(np);
                    face.emplace_back(nv);
                }
            }
        }
        std::reverse(face.begin(), face.end());
        mesh_.add_face(face);
    }

    mesh_.garbage_collection();
    update_mesh();
}

//----------------------------------------------------------------------------

void Viewer::Centroid()
{

    for (Face f : mesh_.faces())
    {
        Vertex hidx = mesh_.add_vertex(centroid(mesh_, f));
        mesh_.split(f, hidx);
    }
    mesh_.garbage_collection();
    update_mesh();
}

void Viewer::open_holes()
{
    for (Face f : holes_)
    {
        mesh_.delete_face(f);
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
        Vertex v = pick_vertex(x, y);
        std::cout << "Vertex Idx: " << v.idx() << std::endl;
        if (mesh_.is_valid(v))
        {
            GeodesicsInHeat heat(mesh_, laplace_matrix, min_point_,
                                 compare_sphere, compare_cube, time_step_);
            Eigen::VectorXd dist, geodist;

            heat.compute_geodesics();
            heat.getDistance(v.idx(), dist, geodist);
            update_mesh();
            mesh_.use_checkerboard_texture();
            set_draw_mode("Texture");
        }
    }
    else
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
