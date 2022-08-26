
#include "LaplaceConstruction.h"
#include "[AW11]Laplace.h"
#include "DiamondLaplace_2D.h"
#include "[dGBD20]Laplace.h"
#include "unsupported/Eigen/SparseExtra"
#include "PolySimpleLaplace.h"
#include "diffgeo.h"

//=============================================================================

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

using namespace std;

enum LaplaceMethods
{
    PolySimpleLaplace = 0,
    AlexaWardetzkyLaplace = 1,
    CotanLaplace = 2,
    Diamond = 3,
    deGoesLaplace = 4
};


void setup_stiffness_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &S,
                              int Laplace, int minpoint)
{

    if (Laplace == AlexaWardetzkyLaplace)
    {
        setup_poly_Laplace_matrix(mesh, S);
        S *= 0.5;
    }
    else if (Laplace == CotanLaplace)
    {
        setup_triangle_Laplace_matrix(mesh, S);
    }
    else if (Laplace == PolySimpleLaplace)
    {
        setup_stiffness_matrix(mesh, S, minpoint);
    }
    else if (Laplace == Diamond)
    {
        setup_face_point_properties(mesh, minpoint);
        Eigen::SparseMatrix<double> G, D, Gra, Div, P;
        setup_prolongation_matrix(mesh, P);
        compute_primal_points(mesh, minpoint);
        setup_diamond_gradient_divergence_intrinsic(mesh, G, D);
        Gra = G * P;
        Div = P.transpose() * D;
        S = Div * Gra;
    }
    else if (Laplace == deGoesLaplace) {
        setup_disney_laplace_operator(mesh, S);
    }
}

//----------------------------------------------------------------------------------

void setup_mass_matrices(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &M,
                         int Laplace, int minpoint,bool lumped)
{

    if (Laplace == AlexaWardetzkyLaplace)
    {
        setup_poly_mass_matrix(mesh, M);
    }
    else if (Laplace == CotanLaplace)
    {
        setup_triangle_mass_matrix(mesh, M);
    }
    else if (Laplace == PolySimpleLaplace)
    {
        setup_mass_matrix(mesh, M, minpoint);
        if (lumped)
        {
            lump_matrix(M);
        }
    }
    else if (Laplace == Diamond)
    {
        Eigen::SparseMatrix<double> M_, P;
        setup_prolongation_matrix(mesh, P);
        setup_diamond_mass_matrix(mesh, M_);
        M = P.transpose() * M_ * P;
    }
    else if (Laplace == deGoesLaplace)
    {
        setup_disney_mass_matrix(mesh, M);
    }
    double area = 0.0;
    for (int k = 0; k < M.outerSize(); ++k)
    {
        for (SparseMatrix::InnerIterator it(M, k); it; ++it)
        {
            area += it.value();
        }
    }
    std::cout << "Surface area: " << area << std::endl;
}
//----------------------------------------------------------------------------------

void setup_gradient(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &G, int Laplace, int minpoint ){
    if (Laplace == AlexaWardetzkyLaplace)
    {
        setup_poly_gradient_operator(mesh, G);
        G *= 0.5;
    }
    else if (Laplace == deGoesLaplace)
    {
        setup_disney_gradient_operator(mesh, G);
    }
    else if (Laplace == PolySimpleLaplace)
    {
        G.resize(3 * mesh.n_faces(), mesh.n_vertices());
        setup_sandwich_gradient_matrix(mesh, G, minpoint);
    }
}
//----------------------------------------------------------------------------------

void setup_divergence(SurfaceMesh &mesh, Eigen::SparseMatrix<double> &D, int Laplace, int minpoint ){
    if (Laplace == AlexaWardetzkyLaplace)
    {
        setup_poly_divergence_operator(mesh, D);
    }
    else if (Laplace == deGoesLaplace)
    {
        setup_disney_divergence_operator(mesh, D);
    }
    else
    {
        D.resize(mesh.n_vertices(), 3 * mesh.n_faces());
        setup_sandwich_divergence_matrix(mesh, D, minpoint);
    }
}