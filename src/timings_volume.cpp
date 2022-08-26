#include <iostream>
#include <Eigen/Dense>

#include <pmp/Timer.h>
#include <unsupported/Eigen/SparseExtra>
#include "Surface/PolyLaplace.h"
#include "Surface/LaplaceConstruction.h"
#include "Volume/LaplaceConstruction_3D.h"
#include "Volume/diffgeo_3D.h"
#include <Eigen/CholmodSupport>

enum LaplaceMethods {
    Diamond = 0,
    Primal_Laplace = 1,
    Dual_Laplace = 2,
    DDFV_CeVe = 3,
    Sandwich = 4,
    DDFV_CeVeFE = 5,
    AQAPoly = 6
};


enum VolumePoints {
    Quadratic_Volume_ = 0,
    Cell_Centroid_ = 1
};

enum AreaPoints {
    Quadratic_Areas_ = 0,
    Face_Centroid = 1
};

using namespace std;


void testOperator_3d(VolumeMesh &mesh, int laplace, int face_point, int cell_point) {
    Eigen::SparseMatrix<double> S, M;

    setup_3D_stiffness_matrix(mesh, S, laplace, face_point, cell_point);

    if (laplace == Diamond) {
        std::cout << "Diamond Laplace :" << std::endl;

    } else if (laplace == Primal_Laplace) {
        std::cout << "Primal Laplace :" << std::endl;
    } else if (laplace == Dual_Laplace) {
        std::cout << "Dual Laplace :" << std::endl;
    }

    std::cout << "Nonzero Entries Laplace "
              << round(((double) sparsity(S) / (S.cols() * S.rows()) * 100.0) * 1000) / 1000
              << "%" << std::endl;
}

//======================Volume Timings==================================================
void takeTimingLGS_3D(VolumeMesh &mesh, vector<double> &compute, vector<double> &solve, int Laplace, int face_point,
                      int cell_point, bool cholmod) {
//
    Eigen::SparseMatrix<double> L, M;
    int degree = 2;
    setup_3D_stiffness_matrix(mesh, S, Laplace, face_point, cell_point,degree);
    setup_3D_mass_matrix(mesh, M, Laplace, face_point, cell_point,degree);

    int nv = mesh.n_vertices();
    int nf = mesh.n_faces();
    int nc = mesh.n_cells();

    Eigen::MatrixXd B(mesh.n_vertices(), 3);
    if (Laplace == DDFV_CeVeFE) {
        B.resize(mesh.n_vertices() + mesh.n_cells() + mesh.n_faces() + mesh.n_edges(), 3);
    }
    if (Laplace == DDFV_CeVe) {
        B.resize(mesh.n_vertices() + mesh.n_cells(), 3);
    }
    for (auto v : mesh.vertices()) {
        B(v.idx(), 0) = mesh.vertex(v)[0];
        B(v.idx(), 1) = mesh.vertex(v)[1];
        B(v.idx(), 2) = mesh.vertex(v)[2];
    }

    if (Laplace == DDFV_CeVeFE) {
        for (auto f: mesh.faces()) {
            VolumeMesh::PointT f_bary = mesh.barycenter(f);
            B(nv + f.idx(), 0) = f_bary[0];
            B(nv + f.idx(), 1) = f_bary[1];
            B(nv + f.idx(), 2) = f_bary[2];
        }
        for (auto c: mesh.cells()) {
            VolumeMesh::PointT c_bary = mesh.barycenter(c);
            B(nv + nf + c.idx(), 0) = c_bary[0];
            B(nv + nf + c.idx(), 1) = c_bary[1];
            B(nv + nf + c.idx(), 2) = c_bary[2];
        }
        for (auto e: mesh.edges()) {
            VolumeMesh::PointT e_bary = mesh.barycenter(e);
            B(nv + nf + nc + e.idx(), 0) = e_bary[0];
            B(nv + nf + nc + e.idx(), 1) = e_bary[1];
            B(nv + nf + nc + e.idx(), 2) = e_bary[2];
        }
    }
    if (Laplace == DDFV_CeVe) {
        for (auto c: mesh.cells()) {
            VolumeMesh::PointT c_bary = mesh.barycenter(c);
            B(nv + c.idx(), 0) = c_bary[0];
            B(nv + c.idx(), 1) = c_bary[1];
            B(nv + c.idx(), 2) = c_bary[2];
        }
    }

    Eigen::SparseMatrix<double> A = M - 0.1 * L;
    Eigen::MatrixXd M_B = M * B;

//    //------------Cholmod----------------------
    if (cholmod) {
//    Eigen::CholmodSimplicialLLT<Eigen::SparseMatrix<double> > cd;
        Eigen::CholmodSupernodalLLT<Eigen::SparseMatrix<double> > cd;
        pmp::Timer t;

        t.start();
        for (int i = 0; i < 20.0; i++) {
            cd.compute(A);
        }
        t.stop();
        std::cout << "Factorization of system matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
        compute.emplace_back(t.elapsed() / 20.0);
        t.start();
        for (int i = 0; i < 20.0; i++) {
            Eigen::MatrixXd X = cd.solve(M_B);
        }
        t.stop();
        std::cout << "Time for solving CMCf : " << t.elapsed() / 20.0 << "ms" << std::endl;
        solve.emplace_back(t.elapsed() / 20.0);
        //-----------------------------------------
    } else {

        static Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
        pmp::Timer t;

        t.start();
        for (int i = 0; i < 20; i++) {
            solver.compute(A);
        }
        t.stop();
        std::cout << "Factorization of system matrix : " << t.elapsed() / 20.0 << "ms" << std::endl;
        compute.emplace_back(t.elapsed() / 20.0);

        t.start();
        for (int i = 0; i < 20; i++) {
            Eigen::MatrixXd X = solver.solve(M_B);
        }
        t.stop();
        std::cout << "Time for solving CMCf : " << t.elapsed() / 20.0 << "ms" << std::endl;
        solve.emplace_back(t.elapsed() / 20.0);

    }


}


double stiffness_sparsity(VolumeMesh &mesh, unsigned int laplace, unsigned int face_point = 0, unsigned int cell_point = 0) {
    Eigen::SparseMatrix<double> S;
    setup_3D_stiffness_matrix(mesh, S,laplace, face_point, cell_point);

//    double val = round(((double) sparsity(S) / (S.cols() * S.rows()) * 100.0) * 1000) / 1000;
    double dim = S.rows()*S.cols();
    double nonzeros= S.nonZeros();
    std::cout << "Entries: " << dim << " non Zeros: " << nonzeros << std::endl;
    return nonzeros/dim;
}


void compute_sparsity() {
    std::string filename_;
    filename_ = "./volume_sparsity.csv";
    std::ofstream file_sparse(filename_);
    std::vector<double> compute_ceve, compute_diamond, compute_sandwich, compute_cevefe;
    file_sparse << "Laplace,Hexas,Pyramids,Truncated,Tetrahedra\n";

    VolumeMesh mesh;

    mesh.read("../data/Volume_data/cubes/cube_hexahedron_5.ovm");
    compute_ceve.emplace_back(stiffness_sparsity(mesh, DDFV_CeVe, Quadratic_Volume_,Quadratic_Areas_));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, Quadratic_Volume_,Quadratic_Areas_));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, Sandwich, Quadratic_Volume_,Quadratic_Areas_));
    compute_cevefe.emplace_back(stiffness_sparsity(mesh, DDFV_CeVeFE, Quadratic_Volume_,Quadratic_Areas_));

    mesh.read("../data/Volume_data/cubes/cube_pyramids_5.ovm");
    compute_ceve.emplace_back(stiffness_sparsity(mesh, DDFV_CeVe, Quadratic_Volume_,Quadratic_Areas_));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, Quadratic_Volume_,Quadratic_Areas_));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, Sandwich, Quadratic_Volume_,Quadratic_Areas_));
    compute_cevefe.emplace_back(stiffness_sparsity(mesh, DDFV_CeVeFE, Quadratic_Volume_,Quadratic_Areas_));


    mesh.read("../data/Volume_data/cubes/cube_truncated_5.ovm");
    compute_ceve.emplace_back(stiffness_sparsity(mesh, DDFV_CeVe, Quadratic_Volume_,Quadratic_Areas_));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, Quadratic_Volume_,Quadratic_Areas_));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, Sandwich, Quadratic_Volume_,Quadratic_Areas_));
    compute_cevefe.emplace_back(stiffness_sparsity(mesh, DDFV_CeVeFE, Quadratic_Volume_,Quadratic_Areas_));


    mesh.read("../data/Tetgen/cubes/cube3.mesh");
    compute_ceve.emplace_back(stiffness_sparsity(mesh, DDFV_CeVe, Quadratic_Volume_,Quadratic_Areas_));
    compute_diamond.emplace_back(stiffness_sparsity(mesh, Diamond, Quadratic_Volume_,Quadratic_Areas_));
    compute_sandwich.emplace_back(stiffness_sparsity(mesh, Sandwich, Quadratic_Volume_,Quadratic_Areas_));
    compute_cevefe.emplace_back(stiffness_sparsity(mesh, DDFV_CeVeFE, Quadratic_Volume_,Quadratic_Areas_));



    file_sparse << "Diamond,";
    for (double i : compute_diamond) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "Sandwich,";
    for (double i : compute_sandwich) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "CeVe,";
    for (double i : compute_ceve) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse << "CeVeFE,";
    for (double i : compute_cevefe) {
        file_sparse << i << ",";
    }
    file_sparse << "\n";

    file_sparse.close();


}

void
test_volume_Meshes(vector<double> &compute, vector<double> &solve, int laplace, int minpoint = 0, int cellpoint = 0,
                   bool cholmod = true) {
    using namespace std;

    VolumeMesh mesh;
//    if (laplace != Primal_Laplace && laplace != Dual_Laplace) {
        std::cout << " ----------------------------------------" << std::endl;
        std::cout << ": Hexaeder " << endl;
        std::cout << " ----------------------------------------" << std::endl;
        mesh.read("../data/Volume_data/cubes/cube_hexahedron_5.ovm");
        std::cout << "mesh n vertices: " << mesh.n_vertices() << std::endl;
        takeTimingLGS_3D(mesh, compute, solve, laplace, minpoint, cellpoint, cholmod);
//        testOperator_3d(mesh, laplace, minpoint, cellpoint);
        std::cout << " ----------------------------------------" << std::endl;
        std::cout << ": Pyramids " << endl;
        std::cout << " ----------------------------------------" << std::endl;
        mesh.read("../data/Volume_data/cubes/cube_pyramids_5.ovm");
        std::cout << "mesh n vertices: " << mesh.n_vertices() << std::endl;
        takeTimingLGS_3D(mesh, compute, solve, laplace, minpoint, cellpoint, cholmod);
//        testOperator_3d(mesh, laplace, minpoint, cellpoint);
        std::cout << " ----------------------------------------" << std::endl;
        std::cout << ": Truncated " << endl;
        std::cout << " ----------------------------------------" << std::endl;
        mesh.read("../data/Volume_data/cubes/cube_truncated_5.ovm");
        std::cout << "mesh n vertices: " << mesh.n_vertices() << std::endl;
        takeTimingLGS_3D(mesh, compute, solve, laplace, minpoint, cellpoint, cholmod);
//        testOperator_3d(mesh, laplace, minpoint, cellpoint);
//    }
    std::cout << " ----------------------------------------" << std::endl;
    std::cout << ":Tetraeder: " << endl;
    std::cout << " ----------------------------------------" << std::endl;
    mesh.read("../data/Tetgen/cubes/cube3.mesh");
    std::cout << "mesh n vertices: " << mesh.n_vertices() << std::endl;
    takeTimingLGS_3D(mesh, compute, solve, laplace, minpoint, cellpoint, cholmod);
//    testOperator_3d(mesh, laplace, minpoint, cellpoint);

}

int main() {
//
//    std::string compute_eigen, compute_cholmod, solve_eigen, solve_cholmod;
//    compute_cholmod = "./volume_timings_compute_cholmod_supernodal.csv";
////    compute_eigen = "./volume_timings_compute_cholmod_llt.csv";
//    compute_eigen = "./volume_timings_compute_eigen_llt.csv";
//    solve_cholmod = "./volume_timings_solve_cholmod_supernodal.csv";
////    compute_cholmod = "./volume_timings_solve_cholmod_llt.csv";
//    solve_eigen = "./volume_timings_solve_eigen_llt.csv";
//
//    std::ofstream file_compute_eigen(compute_eigen);
//    std::ofstream file_compute_cholmod(compute_cholmod);
//    std::ofstream file_solve_eigen(solve_eigen);
//    std::ofstream file_solve_cholmod(solve_cholmod);
//
//    file_compute_eigen << " ,Hexas,Pyramids,Truncated,Tetraheder,\n";
//    file_compute_cholmod << " ,Hexas,Pyramids,Truncated,Tetraheder,\n";
//    file_solve_eigen << " ,Hexas,Pyramids,Truncated,Tetraheder,\n";
//    file_solve_cholmod << " ,Hexas,Pyramids,Truncated,Tetraheder,\n";
//    {
//
//        std::vector<double> compute_primal, compute_diamond, compute_dual, compute_sandwich, compute_CeVeFE, compute_CeVe;
//        std::vector<double> solve_primal, solve_diamond, solve_dual, solve_sandwich, solve_CeVeFE, solve_CeVe;
//
//        bool cholmod = true;
//        std::cout << "Diamond:" << std::endl;
//        test_volume_Meshes(compute_diamond, solve_diamond, Diamond, Quadratic_Areas_, Quadratic_Volume_, cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "CeVe:" << std::endl;
//        test_volume_Meshes(compute_CeVe, solve_CeVe, DDFV_CeVe, Quadratic_Areas_, Quadratic_Volume_, cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "CeVeFE:" << std::endl;
//        test_volume_Meshes(compute_CeVeFE, solve_CeVeFE, DDFV_CeVeFE, Quadratic_Areas_, Quadratic_Volume_, cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "Sandwich:" << std::endl;
//        test_volume_Meshes(compute_sandwich, solve_sandwich, Sandwich, Quadratic_Areas_, Quadratic_Volume_, cholmod);
//        cout << " ==================================" << std::endl;
//        file_compute_cholmod << "Diamond,";
//        for (int i = 0; i < compute_diamond.size(); i++) {
//            file_compute_cholmod << compute_diamond[i] << ",";
//        }
//        file_compute_cholmod << "\n";
//        file_compute_cholmod << "Sandwich,";
//        for (int i = 0; i < compute_sandwich.size(); i++) {
//            file_compute_cholmod << compute_sandwich[i] << ",";
//        }
//        file_compute_cholmod << "\n";
//        file_compute_cholmod << "CeVeFE,";
//        for (int i = 0; i < compute_CeVeFE.size(); i++) {
//            file_compute_cholmod << compute_CeVeFE[i] << ",";
//        }
//        file_compute_cholmod << "\n";
//        file_compute_cholmod << "CeVe,";
//        for (int i = 0; i < compute_CeVe.size(); i++) {
//            file_compute_cholmod << compute_CeVe[i] << ",";
//        }
//        file_compute_cholmod << "\n";
//
//
//        file_solve_cholmod << "Diamond,";
//        for (int i = 0; i < solve_diamond.size(); i++) {
//            file_solve_cholmod << solve_diamond[i] << ",";
//        }
//        file_solve_cholmod << "\n";
//        file_solve_cholmod << "Sandwich,";
//        for (int i = 0; i < solve_sandwich.size(); i++) {
//            file_solve_cholmod << solve_sandwich[i] << ",";
//        }
//        file_solve_cholmod << "\n";
//        file_solve_cholmod << "CeVeFE,";
//        for (int i = 0; i < solve_CeVeFE.size(); i++) {
//            file_solve_cholmod << solve_CeVeFE[i] << ",";
//        }
//        file_solve_cholmod << "\n";
//        file_solve_cholmod << "CeVe,";
//        for (int i = 0; i < solve_CeVe.size(); i++) {
//            file_solve_cholmod << solve_CeVe[i] << ",";
//        }
//        file_solve_cholmod << "\n";
//        file_solve_cholmod.close();
//        file_compute_cholmod.close();
//    }
//    //////////////////////////////////////////////////////////
//    {
//        std::vector<double> compute_primal, compute_diamond, compute_dual, compute_sandwich, compute_CeVeFE, compute_CeVe;
//        std::vector<double> solve_primal, solve_diamond, solve_dual, solve_sandwich, solve_CeVeFE, solve_CeVe;
//
//        bool cholmod = true;
//        std::cout << "Diamond:" << std::endl;
//        test_volume_Meshes(compute_diamond, solve_diamond, Diamond, Quadratic_Areas_, Quadratic_Volume_, !cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "CeVe:" << std::endl;
//        test_volume_Meshes(compute_CeVe, solve_CeVe, DDFV_CeVe, Quadratic_Areas_, Quadratic_Volume_, !cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "CeVeFE:" << std::endl;
//        test_volume_Meshes(compute_CeVeFE, solve_CeVeFE, DDFV_CeVeFE, Quadratic_Areas_, Quadratic_Volume_, !cholmod);
//        cout << " ==================================" << std::endl;
//        std::cout << "Sandwich:" << std::endl;
//        test_volume_Meshes(compute_sandwich, solve_sandwich, Sandwich, Quadratic_Areas_, Quadratic_Volume_, !cholmod);
//        cout << " ==================================" << std::endl;
//        file_compute_eigen << "Diamond,";
//        for (int i = 0; i < compute_diamond.size(); i++) {
//            file_compute_eigen << compute_diamond[i] << ",";
//        }
//        file_compute_eigen << "\n";
//        file_compute_eigen << "Sandwich,";
//        for (int i = 0; i < compute_sandwich.size(); i++) {
//            file_compute_eigen << compute_sandwich[i] << ",";
//        }
//        file_compute_eigen << "\n";
//        file_compute_eigen << "CeVeFE,";
//        for (int i = 0; i < compute_CeVeFE.size(); i++) {
//            file_compute_eigen << compute_CeVeFE[i] << ",";
//        }
//        file_compute_eigen << "\n";
//        file_compute_eigen << "CeVe,";
//        for (int i = 0; i < compute_CeVe.size(); i++) {
//            file_compute_eigen << compute_CeVe[i] << ",";
//        }
//        file_compute_eigen << "\n";
//
//        //////////////////////////////////////////////////////////
//
//        file_solve_eigen << "Diamond,";
//        for (int i = 0; i < solve_diamond.size(); i++) {
//            file_solve_eigen << solve_diamond[i] << ",";
//        }
//        file_solve_eigen << "\n";
//        file_solve_eigen << "Sandwich,";
//        for (int i = 0; i < solve_sandwich.size(); i++) {
//            file_solve_eigen << solve_sandwich[i] << ",";
//        }
//        file_solve_eigen << "\n";
//        file_solve_eigen << "CeVeFE,";
//        for (int i = 0; i < solve_CeVeFE.size(); i++) {
//            file_solve_eigen << solve_CeVeFE[i] << ",";
//        }
//        file_solve_eigen << "\n";
//        file_solve_eigen << "CeVe,";
//        for (int i = 0; i < solve_CeVe.size(); i++) {
//            file_solve_eigen << solve_CeVe[i] << ",";
//        }
//        file_solve_eigen << "\n";
//        file_solve_eigen.close();
//        file_compute_eigen.close();
//    }
    compute_sparsity();
}
