//=============================================================================

#include <imgui.h>
#include <VolumeSubdivision.h>
#include "VolumeViewer.h"
#include "Volume/Eigenmodes.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "VolumeMeshIO.h"
#include "Volume/HarmonicBasis.h"

//=============================================================================

enum VolumePoints
{
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints
{
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

enum LaplaceMethods
{
    Diamond = 0,
    Harmonic = 1,
    PolySimpleLaplace = 2,
};

bool VolumeViewer::load_mesh(const char* filename)
{
    if (VolumeMeshViewer::load_mesh(filename))
    {
        filename_ = filename;
        return true;
    }
    return false;
}

void VolumeViewer::keyboard(int key, int code, int action, int mod)
{
    if (action != GLFW_PRESS) // only react on key press events
        return;

    switch (key)
    {
        case GLFW_KEY_B:
        {
            auto mesh = static_cast<VolumeMesh>(mesh_);
            std::string name = filename_;
            int last_slash = name.find_last_of('/');
            int last_dot = name.find_last_of('.');
            name = name.substr(last_slash + 1, last_dot - last_slash - 1);

            file_manager_.writeFile(name + ".ovm", mesh);
            VolumeMeshIO output(name + ".HYBRID");
            output.write(mesh);
            break;
        }

        default:
        {
            VolumeMeshViewer::keyboard(key, code, action, mod);
            break;
        }
    }
}
//----------------------------------------------------------------------------

void VolumeViewer::process_imgui()
{
    VolumeMeshViewer::process_imgui();

    ImGui::Spacing();
    ImGui::Spacing();

    static int laplace = 0;
    ImGui::RadioButton("Diamond", &laplace, 0);
    ImGui::RadioButton("Harmonic", &laplace, 1);
    ImGui::RadioButton("Polysimple Laplace", &laplace, 2);
    ImGui::Spacing();
    laplace_matrix = laplace;
    ImGui::Spacing();
    ImGui::Spacing();

    static int face_point = 0;
    ImGui::RadioButton("Area minimizer", &face_point, Quadratic_Areas_);
    ImGui::RadioButton("Face Centroid", &face_point, Face_Centroid);
    face_point_ = face_point;
    ImGui::Spacing();
    ImGui::Spacing();

    static int cell_point = 0;
    ImGui::RadioButton("Volume minimizer", &cell_point, Quadratic_Volume_);
    ImGui::RadioButton("Cell Centroid", &cell_point, Cell_Centroid_);
    cell_point_ = cell_point;

    if (ImGui::Button("RMSE Eigenvalues Sphere"))
    {
        Eigen::VectorXd evalues;
        solve_eigenvalue_problem(mesh_, evalues, laplace_matrix, face_point_,
                                 cell_point_);
    }
    if (ImGui::Button("RMSE Franke Poisson System"))
    {
        if (laplace == Harmonic)
        {
            solve_3D_Franke_harmonic(filename_);
        }
        else
        {
            solve_franke_poisson(mesh_, laplace_matrix, face_point_,
                                 cell_point_);
        }
    }
    ImGui::Spacing();
    ImGui::Spacing();
    if (ImGui::CollapsingHeader("Polyhedra!"))
    {
        if (ImGui::Button("Virtual Points"))
        {
            VolumeSubdivision(mesh_).tetrahedra(face_point_, cell_point_);
            update_mesh();
        }
        if (ImGui::Button("Irregular Pyramids"))
        {
            VolumeSubdivision(mesh_).irregular_mesh(5);
            update_mesh();
        }
        if (ImGui::Button("Full Truncation"))
        {
            VolumeSubdivision(mesh_).full_truncation();
            update_mesh();
        }
        if (ImGui::Button("Quads"))
        {
            VolumeSubdivision(mesh_).quads();
            update_mesh();
        }
        if (ImGui::Button("Linear Subdivision"))
        {
            VolumeSubdivision(mesh_).linear_subdivision();
            update_mesh();
        }
    }
}

void VolumeViewer::mouse(int button, int action, int mods)
{
    VolumeMeshViewer::mouse(button, action, mods);
}
