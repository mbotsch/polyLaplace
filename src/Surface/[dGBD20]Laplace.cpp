//=============================================================================
// Copyright 2023 Astrid Bunge, Mario Botsch.
// Distributed under MIT license, see file LICENSE for details.
//=============================================================================

#include "[dGBD20]Laplace.h"
#include "[AW11]Laplace.h"
#include "diffgeo.h"

//=============================================================================

using SurfaceMesh = pmp::SurfaceMesh;
using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

//=============================================================================

// recommended in paper
float deGoes_laplace_lambda_ = 1.0;

//=============================================================================

void setup_deGoes_gradient_operator(pmp::SurfaceMesh& mesh,
                                    Eigen::SparseMatrix<double>& G)
{
    // ------------------------------------ Gradient in Paper -----------------------------------------
    G.resize(3 * (int)mesh.n_faces(), (int)mesh.n_vertices());

    std::vector<Eigen::Triplet<double>> triplets;
    for (auto f : mesh.faces())
    {
        const int n = (int)mesh.valence(f);
        Eigen::VectorXi F(n);

        int i = 0;
        for (auto v : mesh.vertices(f))
        {
            F(i) = (int)v.idx();
            ++i;
        }

        Eigen::MatrixXd E, B, Avg, EtAvg;
        Avg = Eigen::MatrixXd::Zero(n, n);
        for (i = 0; i < (int)n; i++)
        {
            if (i == (int)n - 1)
            {
                Avg(i, 0) = 0.5;
            }
            else
            {
                Avg(i, i + 1) = 0.5;
            }
            Avg(i, i) = 0.5;
        }
        setup_E_and_B_perFace(mesh, f, E, B);
        EtAvg.resize(3, n);
        EtAvg = E.transpose() * Avg;

        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
        Eigen::Vector3d nf = af;
        nf.normalize();

        Eigen::MatrixXd d(n, n);
        Eigen::MatrixXd G_f;
        d.setZero();
        G_f.resize(3, n);
        for (i = 0; i < n; ++i)
        {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }

        double area = -1.0 / af.norm();
        int j = 0;
        for (pmp::Vertex v : mesh.vertices(f))
        {
            Eigen::Vector3d x = EtAvg.col(j);
            Eigen::Vector3d gradient = nf.cross(x);
            G_f.col(j) = area * gradient;
            triplets.emplace_back(3 * f.idx(), v.idx(), area * gradient(0));
            triplets.emplace_back(3 * f.idx() + 1, v.idx(), area * gradient(1));
            triplets.emplace_back(3 * f.idx() + 2, v.idx(), area * gradient(2));
            j++;
        }
    }

    G.setFromTriplets(triplets.begin(), triplets.end());
}

void setup_deGoes_divergence_operator(pmp::SurfaceMesh& mesh,
                                      Eigen::SparseMatrix<double>& D)
{
    //---------------------------- Divergence as Gradient Transpose * area ------------------------------

    Eigen::SparseMatrix<double> G, A;
    setup_deGoes_gradient_operator(mesh, G);

    setup_deGoes_face_area_matrix(mesh, A);
    D.resize((int)mesh.n_vertices(), 3 * (int)mesh.n_faces());
    D = -G.transpose() * A;

    //---------------------------- Divergence in Paper------------------------------

    //    std::vector<Eigen::Triplet<double>> triplets;
    //    int colCnt = 0;
    //
    //    for (auto f : mesh.faces()) {
    //        const unsigned int n = mesh.valence(f);
    //        Eigen::VectorXi F(n);
    //        Eigen::MatrixXd E, B;
    //        setup_E_and_B_perFace(mesh, f, E, B);
    //        int i = 0;
    //        for (auto v : mesh.vertices(f)) {
    //            F(i) = v.idx();
    //            ++i;
    //        }
    //        // compute vector area
    //        const Eigen::Matrix3d Af = E.transpose() * B;
    //        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    //
    //        const double area = af.norm();
    //        af /= area;
    //
    //        // compute d
    //        Eigen::MatrixXd d(n, n);
    //        Eigen::MatrixXd M, Lf;
    //        setup_disney_inner_product_per_face(mesh,f,M);
    //        d.setZero();
    //
    //        for (unsigned int i = 0; i < n; ++i) {
    //            d(i, i) = -1;
    //            d(i, (i + 1) % n) = 1;
    //        }
    //
    //        Lf = d.transpose() * M;
    //
    //        // add local laplacian to global matrix entries
    //        for (unsigned int k = 0; k < n; ++k) {
    //            for (unsigned int l = 0; l < n; ++l) {
    //                triplets.push_back(
    //                        Eigen::Triplet<double>(F(k), colCnt + l, - Lf(k, l)));
    //            }
    //        }
    //
    //        colCnt += n;
    //    }
    //
    //    D.resize(mesh.n_vertices(), colCnt);
    //    D.setFromTriplets(triplets.begin(), triplets.end());
}

void setup_deGoes_sharp_operator_per_face(pmp::SurfaceMesh& mesh, pmp::Face f,
                                          Eigen::MatrixXd& U)
{
    const int n = (int)mesh.valence(f);
    U.resize(3, n);
    Eigen::MatrixXd E, B, Bt;
    setup_E_and_B_perFace(mesh, f, E, B);
    Bt.resize(3, n);
    Bt = B.transpose();
    pmp::Point c = centroid(mesh, f);
    Eigen::Vector3d cf;
    cf << c[0], c[1], c[2];
    // compute vector area
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    Eigen::Vector3d nf = af;
    nf.normalize();

    double area = 1.0 / af.norm();
    for (int i = 0; i < n; i++)
    {
        Eigen::Vector3d x = Bt.col(i) - cf;
        Eigen::Vector3d cross_prod = nf.cross(x);
        U.col(i) = area * cross_prod;
    }
}

void setup_deGoes_discrete_flat_operator_per_face(pmp::SurfaceMesh& mesh,
                                                  pmp::Face f,
                                                  Eigen::MatrixXd& V)
{
    const int n = (int)mesh.valence(f);
    V.resize(n, 3);
    Eigen::MatrixXd E, B, NNT;
    setup_E_and_B_perFace(mesh, f, E, B);

    // compute vector area
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    Eigen::Vector3d nf = af;
    nf.normalize();
    NNT = Eigen::MatrixXd::Identity(3, 3);
    NNT -= nf * nf.transpose();
    V = E * NNT;
}

void setup_deGoes_projection_operator_per_face(pmp::SurfaceMesh& mesh,
                                               pmp::Face f, Eigen::MatrixXd& P)
{
    const unsigned int n = mesh.valence(f);
    P.resize(n, n);
    P.setIdentity();
    Eigen::MatrixXd U, V;
    setup_deGoes_discrete_flat_operator_per_face(mesh, f, V);
    setup_deGoes_sharp_operator_per_face(mesh, f, U);

    P -= V * U;
}

void setup_deGoes_inner_product_per_face(pmp::SurfaceMesh& mesh, pmp::Face f,
                                         Eigen::MatrixXd& M)
{
    const unsigned int n = mesh.valence(f);

    double lambda = deGoes_laplace_lambda_;
    M.resize(n, n);
    Eigen::MatrixXd E, B, U, P;
    setup_E_and_B_perFace(mesh, f, E, B);
    setup_deGoes_sharp_operator_per_face(mesh, f, U);
    setup_deGoes_projection_operator_per_face(mesh, f, P);
    const Eigen::Matrix3d Af = E.transpose() * B;
    Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
    double area = af.norm();

    M = area * U.transpose() * U + lambda * P.transpose() * P;
}

void setup_deGoes_laplace_operator(pmp::SurfaceMesh& mesh,
                                   Eigen::SparseMatrix<double>& L)
{
    L.resize((int)mesh.n_vertices(), (int)mesh.n_vertices());

    std::vector<Eigen::Triplet<double>> triplets;
    for (auto f : mesh.faces())
    {
        const int n = (int)mesh.valence(f);

        Eigen::VectorXi F(n);
        int i = 0;
        for (auto v : mesh.vertices(f))
        {
            F(i) = (int)v.idx();
            ++i;
        }
        Eigen::MatrixXd M, Lf;
        setup_deGoes_inner_product_per_face(mesh, f, M);
        // compute d
        Eigen::MatrixXd d(n, n);
        d.setZero();
        for (i = 0; i < n; ++i)
        {
            d(i, i) = -1;
            d(i, (i + 1) % n) = 1;
        }
        Lf = d.transpose() * M * d;

        // add local laplacian to global matrix entries
        for (int k = 0; k < n; ++k)
        {
            for (int l = 0; l < n; ++l)
            {
                triplets.emplace_back(F(k), F(l), -Lf(k, l));
            }
        }
    }
    L.setFromTriplets(triplets.begin(), triplets.end());
    //std::cout << "deGoes Laplace hyperparameter: " << deGoes_laplace_lambda_
    //          << std::endl;
}

void setup_deGoes_mass_matrix(pmp::SurfaceMesh& mesh,
                              Eigen::SparseMatrix<double>& M)
{
    M.resize((int)mesh.n_vertices(), (int)mesh.n_vertices());

    std::vector<Triplet> tripletsM;
    double sum = 0.0;
    Eigen::MatrixXd E, B;
    for (auto v : mesh.vertices())
    {
        for (auto f : mesh.faces(v))
        {
            const unsigned int n = mesh.valence(f);
            setup_E_and_B_perFace(mesh, f, E, B);
            //          compute vector area
            const Eigen::Matrix3d Af = E.transpose() * B;
            Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
            const double area = af.norm();
            sum += area / (double)(n);
        }
        tripletsM.emplace_back(v.idx(), v.idx(), sum);
        sum = 0.0;
    }
    // build sparse matrix from triplets
    M.setFromTriplets(tripletsM.begin(), tripletsM.end());
}

void setup_deGoes_face_area_matrix(pmp::SurfaceMesh& mesh,
                                   Eigen::SparseMatrix<double>& A)
{
    std::vector<Eigen::Triplet<double>> area_triplets;
    A.resize(3 * (int)mesh.n_faces(), 3 * (int)mesh.n_faces());
    for (auto f : mesh.faces())
    {
        Eigen::MatrixXd E, B, Gf;
        setup_E_and_B_perFace(mesh, f, E, B);

        // compute vector area
        const Eigen::Matrix3d Af = E.transpose() * B;
        Eigen::Vector3d af(-Af(1, 2), Af(0, 2), -Af(0, 1));
        area_triplets.emplace_back(3 * f.idx(), 3 * f.idx(), af.norm());
        area_triplets.emplace_back(3 * f.idx() + 1, 3 * f.idx() + 1, af.norm());
        area_triplets.emplace_back(3 * f.idx() + 2, 3 * f.idx() + 2, af.norm());
    }
    A.setFromTriplets(area_triplets.begin(), area_triplets.end());
}
