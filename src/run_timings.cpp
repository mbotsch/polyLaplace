#include <iostream>
#include <Eigen/Dense>
#include <pmp/Timer.h>
#include <pmp/io/io.h>
#include <unsupported/Eigen/SparseExtra>
#include "Surface/LaplaceConstruction.h"
#include "VolumeMesh.h"
#include "Volume/LaplaceConstruction_3D.h"
#include "Volume/Franke_PoissonSystem_3D.h"
#include "Volume/HarmonicBasis.h"

enum LaplaceMethods2D
{
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    Diamond = 2,
    deGoesLaplace = 3,
    Harmonic = 4
};
enum LaplaceMethods3D
{
    Diamond3D = 0,
    Harmonic3D = 1,
    PolySimpleLaplace3D = 2,

};
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

using namespace std;

double StiffnessConstructionTimings2D(SurfaceMesh& mesh, int laplace)
{
    Eigen::SparseMatrix<double> L;
    pmp::Timer t;
    t.start();

    for (int i = 0; i < 20; i++)
    {
        setup_stiffness_matrices(mesh, L, laplace, Quadratic_Areas_);
    }
    t.stop();
    std::cout << "Time for Laplace Construction : " << t.elapsed() / 20.0
              << "ms" << std::endl;
    return t.elapsed() / 20.0;
}

double StiffnessConstructionTimings3D(VolumeMesh& mesh, int laplace)
{
    Eigen::SparseMatrix<double> L;
    pmp::Timer t;
    t.start();

    for (int i = 0; i < 20; i++)
    {
        setup_3D_stiffness_matrix(mesh, L, laplace, Quadratic_Areas_,
                                  Quadratic_Volume_);
    }
    t.stop();
    std::cout << "Time for Laplace Construction : " << t.elapsed() / 20.0
              << "ms" << std::endl;
    return t.elapsed() / 20.0;
}

void takeLaplaceTiming2D(SurfaceMesh& mesh, int laplace,
                         std::ofstream& file_timings)
{
    auto points = mesh.vertex_property<Point>("v:point");
    Eigen::SparseMatrix<double> L, M, A;

    pmp::Timer t_stiffness, t_solve;

    t_stiffness.start();
    setup_stiffness_matrices(mesh, L, laplace, Quadratic_Areas_);
    t_stiffness.stop();
    std::cout << "Time for Laplace Construction : " << t_stiffness.elapsed()
              << "ms" << std::endl;
    file_timings << t_stiffness.elapsed() << ",";

    setup_mass_matrices(mesh, M, laplace, Quadratic_Areas_);
    Eigen::MatrixXd B(mesh.n_vertices(), 3);

    for (auto v : mesh.vertices())
    {
        B(v.idx(), 0) = points[v][0];
        B(v.idx(), 1) = points[v][1];
        B(v.idx(), 2) = points[v][2];
    }

    A = M - 0.1 * L;
    Eigen::MatrixXd M_B = M * B;

    static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

    pmp::Timer t;

    t_solve.start();
    for (int i = 0; i < 10; i++)
    {
        solver.compute(A);
    }
    t_solve.stop();
    Eigen::SparseMatrix<double> C = solver.matrixL();
    std::cout << "Nonzeros Cholesky: " << C.nonZeros() << std::endl;
    std::cout << "Ratio cholesky to system : "
              << C.nonZeros() / (float)A.nonZeros() << endl;
    std::cout << "Time for compute A : " << t.elapsed() / 10.0 << "ms"
              << std::endl;

    file_timings << L.nonZeros() << ",";
    t_solve.start();
    for (int i = 0; i < 10; i++)
    {
        Eigen::MatrixXd X = solver.solve(M_B);
    }
    t_solve.stop();
    file_timings << t_solve.elapsed() / 10.0 << ",\n";
    std::cout << "Time for solving CMCf : " << t_solve.elapsed() / 10.0 << "ms"
              << std::endl;
}

void takeLaplaceTiming3D(const std::string& meshname, int Laplace,
                         std::ofstream& file_timings)
{
    Eigen::SparseMatrix<double> L, M, A;

    pmp::Timer t_stiffness, t_solve;
    VolumeMesh mesh;
    mesh.read(meshname);
    Eigen::VectorXd b(mesh.n_vertices());

    auto c_prop = mesh.request_cell_property<VolumeMesh::PointT>("cell points");
    auto f_prop = mesh.request_face_property<VolumeMesh::PointT>("face points");
    if (Laplace == Harmonic3D)
    {
        PolyhedralMesh mesh2(meshname);
        t_stiffness.start();
        setup_3D_harmonic_stiffness_matrix(mesh2, L);
        t_stiffness.stop();
        std::cout << "Time for Laplace Construction : " << t_stiffness.elapsed()
                  << "ms" << std::endl;
        file_timings << t_stiffness.elapsed() << ",";
        setup_3D_harmonic_matrices(mesh2, L, M);
        L *= -1.0;
    }
    else
    {
        t_stiffness.start();
        setup_3D_stiffness_matrix(mesh, L, Laplace, Quadratic_Areas_,
                                  Quadratic_Volume_);
        t_stiffness.stop();
        std::cout << "Time for Laplace Construction : " << t_stiffness.elapsed()
                  << "ms" << std::endl;
        file_timings << t_stiffness.elapsed() << ",";
        setup_3D_mass_matrix(mesh, M, Laplace, Quadratic_Areas_,
                             Quadratic_Volume_);
    }

    for (auto v : mesh.vertices())
    {
        b(v.idx()) = laplace_franke3d(mesh.vertex(v));
    }

    A = M - 0.1 * L;
    Eigen::MatrixXd M_B = M * b;

    static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;

    t_solve.start();
    for (int i = 0; i < 10; i++)
    {
        solver.compute(A);
    }
    t_solve.stop();
    Eigen::SparseMatrix<double> C = solver.matrixL();
    std::cout << "Nonzeros Cholesky: " << C.nonZeros() << std::endl;
    std::cout << "Ratio cholesky to system : "
              << C.nonZeros() / (float)A.nonZeros() << endl;
    std::cout << "Time for compute A : " << t_solve.elapsed() / 10.0 << "ms"
              << std::endl;

    file_timings << L.nonZeros() << ",";
    t_solve.start();
    for (int i = 0; i < 10; i++)
    {
        Eigen::MatrixXd X = solver.solve(M_B);
        if (solver.info() != Eigen::Success)
        {
            std::cerr << "FUCK PROBLEM" << std::endl;
        }
    }
    t_solve.stop();
    file_timings << t_solve.elapsed() / 10.0 << ",\n";
    std::cout << "Time for solving CMCf : " << t_solve.elapsed() / 10.0 << "ms"
              << std::endl;
}

void taketimings2D(std::string meshname, std::string path)
{
    std::string filename_ = "./timings_" + meshname + ".csv";
    SurfaceMesh mesh;
    pmp::read(mesh, path);
    std::ofstream file(filename_);

    file << "Method,Build,Nonzeros,Solve,\n";

    file << "Alexa,";
    takeLaplaceTiming2D(mesh, AlexaWardetzkyLaplace, file);

    file << "deGoes,";
    takeLaplaceTiming2D(mesh, deGoesLaplace, file);
    file << "Diamond,";
    takeLaplaceTiming2D(mesh, Diamond, file);

    file << "Polysimple,";
    takeLaplaceTiming2D(mesh, PolySimpleLaplace, file);

    file << "Harmonic,";
    takeLaplaceTiming2D(mesh, Harmonic, file);
    file.close();
}

void taketimings3D(std::string meshname, std::string path)
{
    std::string filename_ = "./timings_" + meshname + ".csv";
    std::ofstream file(filename_);

    file << "Method,Build,Nonzeros,Solve,\n";
    file << "Diamond,";
    takeLaplaceTiming3D(path, Diamond3D, file);

    file << "Polysimple,";
    takeLaplaceTiming3D(path, PolySimpleLaplace3D, file);

    file << "Harmonic,";
    takeLaplaceTiming3D(path, Harmonic3D, file);
    file.close();
}

int main()
{
    std::string path;
    //    Surface timings
    path = "../data/surface_meshes/grid/quad_5.obj";
    taketimings2D("2DQuad", path);
    path = "../data/surface_meshes/grid/clean/voronoi_5.obj";
    taketimings2D("2DVoronoi", path);

    //Volume timings
    path = "../data/volume_meshes/cubes/cube_hexahedron_4.ovm";
    taketimings3D("3DHexahedra", path);
    path = "../data/volume_meshes/cubes/cube_voronoi_4.ovm";
    taketimings3D("3DVoronoi", path);
}
